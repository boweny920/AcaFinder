import subprocess
import sys
import pandas as pd
from Bio import SeqIO
import os
import pprint as pp

def make_output_folder(GCF_Name):
    os.makedirs(GCF_Name)
    return os.path.dirname(GCF_Name)

def get_file_directory_path(file_path):
    directory_path=str(file_path).rsplit("/",1)
    return directory_path[0]

def length(list):
    gene_length=int(list[4])-int(list[3]) + 1
    return gene_length

def contig(list):
    contig=list[0]
    return contig

def strand(list):
    gene_strand=list[6]
    return gene_strand

def proteinInfo(list):
    return list[-1]

def proteinID(IDinfo):
    IDinfo = IDinfo.split(";")
    for v in IDinfo:
        if "protein_id=" in v:
            ProID = str(v).strip('protein_id=')
            return ProID

def protein_pos_start(list):
    protein_start=int(list[3])
    return protein_start

def protein_pos_end(list):
    protein_end = int(list[4])
    return protein_end

def proteinInfo_List_process(proteinInfo_list):
    combined_list=[]
    for sublist in proteinInfo_list:
        loci_list_proteins_WithPseudo = []
        for number,value in enumerate(sublist):
            if "pseudo=true" not in value:
                loci_list_proteins_WithPseudo.append([int(number),proteinID(value)])
            elif "pseudo=true" in value:
                loci_list_proteins_WithPseudo.append([int(number),"pseudo"])
        pseudo_gene_number=0
        for x in loci_list_proteins_WithPseudo:pseudo_gene_number=pseudo_gene_number+x.count("pseudo")
        if len(loci_list_proteins_WithPseudo) - pseudo_gene_number >= 2 :
        ##see how many genes left without the pseudo gene
            combined_list.append(loci_list_proteins_WithPseudo)
    return combined_list

def result_check_list(list):
    gene_strand=strand(list)
    gene_length=length(list)
    if "pseudo=true" not in list[-1]:ProID=proteinID(proteinInfo(list))
    elif "pseudo=true" in list[-1]:ProID="pseudo"
    # output result in format ProteinID; contig;gene strand; gene length；gene start;gene end, in a list
    check_result = [ProID, contig(list), gene_strand, gene_length, list[3], list[4]]
    return check_result

def loci_select_before_diamond(file_dic,all_protein_lenbp=600,intergenic_dist=250):
    ###################################################################
    ###   Get the regions that contain at least two genes, each gene must
    ###   be on the same strand, each gene must be less than 200 aa or 600 nucleotides
    key=0
    lst=[]
    loci_list=[]
    lst_result_check=[]
    loci_list_result_check=[]
    while key in file_dic.keys():
        ##First filter gene length:
        if length(file_dic[key]) < all_protein_lenbp:
            #This is working!!
            if key-1 in file_dic.keys() and strand(file_dic[key]) == strand(file_dic[key-1]) and contig(file_dic[key]) == contig(file_dic[key-1]) and length(file_dic[key-1]) < all_protein_lenbp and protein_pos_start(file_dic[key]) - protein_pos_end(file_dic[key-1]) < intergenic_dist:
                #print(contig(file_dic[key]))
                lst.append(proteinInfo(file_dic[key]))
                # This is for result checking
                lst_result_check.append(result_check_list(file_dic[key]))
                key=key+1
            elif key+1 in file_dic.keys() and strand(file_dic[key])==strand(file_dic[key+1]) and contig(file_dic[key]) == contig(file_dic[key+1]) and length(file_dic[key+1]) < all_protein_lenbp and protein_pos_start(file_dic[key+1]) - protein_pos_end(file_dic[key]) < intergenic_dist:
                if len(lst) >= 2:
                    #if previous lst has 2 or more genes, append the lst to the loci_list
                    loci_list.append(lst)
                    loci_list_result_check.append(lst_result_check)
                lst = []
                lst_result_check = []
                lst.append(proteinInfo(file_dic[key]))
                lst.append(proteinInfo(file_dic[key+1]))
                #This is for result checking
                lst_result_check.append(result_check_list(file_dic[key]))
                lst_result_check.append(result_check_list(file_dic[key+1]))

                key=key+2
            else:
                key = key + 1
                if len(lst) >= 2:
                    # print(lst)
                    loci_list.append(lst)
                    loci_list_result_check.append(lst_result_check)
                lst = []
                lst_result_check = []
        else:
            key=key+1
            if len(lst) >= 2:
                # print(lst)
                loci_list.append(lst)
                # This is for result checking
                loci_list_result_check.append(lst_result_check)
            lst=[]
            lst_result_check = []
    else:
        if len(lst) >= 2:
            # print(lst)
            loci_list.append(lst)
            # This is for result checking
            loci_list_result_check.append(lst_result_check)
        lst = []
        lst_result_check = []

    loci_list_w_pseudo=proteinInfo_List_process(loci_list)
    # pp.pprint(loci_list_w_pseudo)
    # pp.pprint(loci_list_result_check)
    return loci_list_w_pseudo,loci_list_result_check

def is_non_zero_file(fpath):
    #Check if file is empty or does not exsit
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def make_and_check_output_directory(dpath):
    if os.path.isdir(dpath) is False:
        os.makedirs(dpath)
    return str(dpath)

def run_diamond(database_file,outputdir,query,evalue_cutoff,coverage_cut_off,threads): #add coverage cutoff
    ##Use this as dmond_out
    subprocess.Popen(['diamond','makedb','--in',database_file,'-d',os.path.join(outputdir,"diamond_db")]).wait()
    subprocess.Popen(['diamond','blastp','-d',os.path.join(outputdir,"diamond_db"),'-f','6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen' ,'qstart' ,'qend' ,'sstart' ,'send' ,'evalue' ,'bitscore', 'qlen', 'slen','--quiet','--more-sensitive','-q',query,'-p',threads,'-o',os.path.join(outputdir,"diamond_blastp_result.txt"),'-e',evalue_cutoff]).wait()
    ## Output format 6 columns: Default: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
    subprocess.Popen(["awk -F '\t' '($8-$7)/$(NF-1) > %s {print}' %s > %s"%(str(coverage_cut_off),os.path.join(outputdir,"diamond_blastp_result.txt"),os.path.join(outputdir,"diamond_blastp_result.coverageParsed.txt"))],shell=True).wait()
    return os.path.join(outputdir,"diamond_blastp_result.coverageParsed.txt")

def parse_diamond_get_proteinID(dmond_out):
    ##Use with run_diamond; Use on NON-EMPTY diamond output file
    ##Return also a dict of Protein ID to Acr type
    dic_acr={}
    protein_ID_list=[]
    Pro_ID = subprocess.Popen("awk -F '\t' '{print $1,$2}' %s" % (dmond_out), shell=True, stdout=subprocess.PIPE)
    for line in Pro_ID.stdout:
        line = line.rstrip().decode('utf-8').split()
        proID = str(line[1])
        protein_ID_list.append(proID)
        dic_acr.setdefault(proID,str(line[0]))
    return protein_ID_list,dic_acr

def faa_file_wrote(one_acr_aca_locus,faa_file,newfile_name_dirctory):
    ##The one_acr_aca_locus is a list of proteinID.
    ##The function takes the IDs and the faa file with protein sequences
    ##And make a fasta newfile containing ProteinIDs with corresponding amino acid sequences
    record_dict = SeqIO.to_dict(SeqIO.parse(faa_file, "fasta"))
    with open(newfile_name_dirctory,"w") as newfile:
        for proID in one_acr_aca_locus:
            SeqIO.write(record_dict[proID],newfile,"fasta")

def faa_file_wrote_potential_AcrAca(one_acr_aca_locus,faa_file,newfile_name_dirctory):
    record_dict = SeqIO.to_dict(SeqIO.parse(faa_file, "fasta"))
    with open(newfile_name_dirctory, "w") as newfile:
        for gene_order,proID in one_acr_aca_locus:
            if proID in record_dict.keys():
                SeqIO.write(record_dict[proID], newfile, "fasta")

def aca_select_process(one_acr_aca_locus,faa_file,newfile_name_dirctory,Acr_protein_NP_list,Acr_aca_inBetweenGenes):
    ##This function takes the proteins that are not Acrs and write them into a new faa file
    ##The one_acr_aca_locus is a list of [gene order, proteinID].
    ##The function takes the IDs and the faa file with protein sequences
    ##And make a fasta newfile containing ProteinIDs with gene orders with corresponding amino acid sequences
    record_dict = SeqIO.to_dict(SeqIO.parse(faa_file, "fasta"))
    unhmmer_aca = []
    with open(newfile_name_dirctory, "w") as newfile:
        acr_location=[]
        for order, proID in one_acr_aca_locus:
            if proID in Acr_protein_NP_list: acr_location.append(int(order)+1)
            # "+1" is to get rid of 0, which makes the stops the function when occured
        for order1, proID1 in one_acr_aca_locus:
            order1=order1+1
            # "+1" is to get rid of 0, which makes the stops the function when occured
            #if proID1 not in Acr_protein_NP_list and any(v for v in acr_location if abs(v-int(order1))<=4):
            if any(v for v in acr_location if abs(v - int(order1)) <= int(Acr_aca_inBetweenGenes)):
                if proID1 in record_dict.keys():
                    SeqIO.write(record_dict[proID1], newfile, "fasta")
                    unhmmer_aca.append(proID1)
                ##Check the gene distance between Acr and Aca, it should be less than 3 genes
    return unhmmer_aca

def run_hmmscan(faafile,evalue_cut_off,hmmfile,outdir,threads):
    hmm_outfile=faafile+".hmmout"
    subprocess.Popen(['hmmscan','--domtblout',hmm_outfile,'-o',os.path.join(outdir,'log.hmm'),'--noali',"--cpu",threads,'-E', evalue_cut_off,hmmfile,faafile]).wait()
    return hmm_outfile

def parse_hmmOutfile(hmm_outfile,hmm_coverage_cutoff): # redo the coverage calculation part!
    parsed_hmm_outfile = str(hmm_outfile) + ".Coverage_parsed"
    #Parse the hmm output in terms of subject coverage: subjust aligned length/ Total subject length
    subprocess.Popen(["grep -v '#' %s | awk '($17-$16)/$3 > %s {print}'  > %s" % (hmm_outfile, str(hmm_coverage_cutoff),parsed_hmm_outfile)],shell=True).wait()
    return parsed_hmm_outfile

def potential_new_aca_faa_filemake(parsed_hmm_outfile,newAcrAca_faaFile,unhmmer_aca):
    ## Make a fasta file of potential Aca amino acid sequences
    ## Make a list of Aca, for checking results
    Aca_Pro_lst_dic={}
    record_dict = SeqIO.to_dict(SeqIO.parse(newAcrAca_faaFile, "fasta"))
    newfile_name= parsed_hmm_outfile+".new_found_ACA.faa"
    with open(newfile_name,"w") as new:
        acaout=subprocess.Popen(["awk '{print $1,$2,$4}' %s |sort -u"%parsed_hmm_outfile],shell=True, stdout=subprocess.PIPE)
        for line in acaout.stdout:
            HTH,pfamID,ID=line.rstrip().decode('utf-8').split()
            if ID in unhmmer_aca:
                # If Aca is not in unhmmer, then it should not be considered as a potential Aca(here is to make sure that condiction "gene distance < 3" is applied)
                if ID not in Aca_Pro_lst_dic.keys():
                    Aca_Pro_lst_dic.setdefault(ID,[HTH+"="+pfamID])
                    SeqIO.write(record_dict[ID], new, "fasta")
                elif ID in Aca_Pro_lst_dic.keys():
                    Aca_Pro_lst_dic[ID].append(HTH+"="+pfamID)

                #This is to consider the < than 3 genes condition, which had been implemented before in "faa_file_wrote_diamond_special" function
    return newfile_name,Aca_Pro_lst_dic

def final_result_check_output_generation(Acr_homolog_candidate_list_resultCheck,newfile_name,Aca_Pro_lst_dic,publishedAcaHMM_hits_dic):
    #generation of final result checking file
    #Table will be tsv format
    #### GCF_number    ProteinID   Contig  Strand  Length  Start   End     Species     Acr     Aca  AcaHMM_HIT ####
    Acr_Aca_loci_by_GBA = []
    for i in Acr_homolog_candidate_list_resultCheck:
        if any(v for v in i if v[0] in Aca_Pro_lst_dic):
            for v in i:
                # Check which ProteinID is the Aca homolog
                if v[0] in Aca_Pro_lst_dic:
                    v.append("Aca_protein|"+";".join(Aca_Pro_lst_dic[v[0]]))
                else:
                    v.append("NA")
                # Check which ProteinID has AcaHMM hits
                if v[0] in publishedAcaHMM_hits_dic:
                    v.append(";".join(publishedAcaHMM_hits_dic[v[0]]))
                else:
                    v.append("NA")
            Acr_Aca_loci_by_GBA.append(i)

    with open(newfile_name, "w") as newfile:
        for pro in Acr_Aca_loci_by_GBA:
            for value in pro:
                # for info in value: newfile.write("%s\t" % str(info))
                newfile.write("\t".join([str(v) for v in value]))
                newfile.write("\n")

def pos_rearrange(pos_list):
    if pos_list[0] < pos_list[1]: return pos_list
    else:
        tmp=pos_list[0]
        pos_list[0]=pos_list[1]
        pos_list[1]=tmp
        return pos_list

def distance_cal(pos1,pos2):
    #map function, very useful
    pos1=pos_rearrange([int(v) for v in pos1.split("-")])
    pos2=pos_rearrange([int(v) for v in pos2.split("-")])
    # Handy function to calculate distances
    if pos2[0] >pos1[1]: return pos2[0]-pos1[1]
    if pos1[0]<pos2[0]<pos1[1] or pos1[0]<pos2[1]<pos1[1] : return 0
    if pos2[1]<pos1[0]:return pos1[0]-pos2[1]
    else: return 0

def find_complete_CRISPR_Cas_and_SelfTargeting(fna,outputfile,threads,general_folder):
    subprocess.Popen(["cctyper",fna,outputfile,"--no_plot","-t",threads,"--prodigal","meta"]).wait()
    CRISPR_CasTable = subprocess.run(["find", outputfile, "-name", "CRISPR_Cas.tab"], stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip()
    if CRISPR_CasTable == "":
        print("No complete CRISPR-Cas found input")
        return None
    else:
        print("Complete CRISPR-Cas found and can be found in %s"%CRISPR_CasTable)
        CasTable = subprocess.run(["find", outputfile, "-name", "cas_operons.tab"],
                                  stdout=subprocess.PIPE).stdout.decode(
            'utf-8').rstrip()
        Crispr_Table = subprocess.run(["find", outputfile, "-name", "crisprs_near_cas.tab"],
                                      stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip()
        df_CasTable = pd.read_csv(CasTable, sep="\t")
        df_Crispr_Table = pd.read_csv(Crispr_Table, sep="\t")
        fna_dic = SeqIO.to_dict(SeqIO.parse(fna, "fasta"))

        df_CRISPRcas = pd.read_csv(CRISPR_CasTable, sep="\t")
        CC_list = []
        fna_blastdb = os.path.join(outputfile, os.path.basename(fna) + ".blastDB")
        subprocess.Popen(["makeblastdb", "-dbtype", "nucl", "-in", fna, "-out", fna_blastdb]).wait()
        df_CCtable = pd.DataFrame()
        CC_contig = []
        CC_location = []
        C_operon_withLocation = []
        CC_type = []
        Cas_operon_withLocation = []
        STSS_region = []
        contig_length = []
        for index, CRISPR_Cas in df_CRISPRcas.iterrows():
            location = "-".join(CRISPR_Cas["Operon_Pos"].lstrip("[").rstrip("]").replace(" ", "").split(","))
            CRCasType = CRISPR_Cas["Prediction"]
            spacer_names = CRISPR_Cas["CRISPRs"].lstrip("[").rstrip("]").replace("'", "").replace(" ", "").split(
                ",")  # Some CC have multiple spacers, thisis to get all the spacers
            # print("CRISPR-Cas spacers: ")
            ## Look for self-targeting regions in the genome
            self_targeting_regions = []
            CC_contig.append(CRISPR_Cas["Contig"])
            CC_type.append(CRCasType)
            CC_location.append(location)
            tmp_list = []
            for spacer in spacer_names:
                tmp_list.append(spacer + "|" + "-".join(
                    [df_Crispr_Table[df_Crispr_Table["CRISPR"] == spacer]["Start"].to_string(index=False),
                     df_Crispr_Table[df_Crispr_Table["CRISPR"] == spacer]["End"].to_string(index=False)]))
            C_operon_withLocation.append(";".join(tmp_list))
            Cas_operon_withLocation.append(CRISPR_Cas["Operon"] + "|" + "-".join(
                [df_CasTable[df_CasTable["Operon"] == CRISPR_Cas["Operon"]]["Start"].to_string(index=False),
                 df_CasTable[df_CasTable["Operon"] == CRISPR_Cas["Operon"]]["End"].to_string(index=False)]))
            contig_length.append(len(fna_dic[CRISPR_Cas["Contig"]].seq))
            for spacer in spacer_names:
                spacer_fna = subprocess.run(["find", outputfile, "-name", spacer + ".fa"],
                                            stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip()
                spacer_fna_blastOutfile = os.path.join(outputfile,
                                                       os.path.basename(fna) + "_VS_" + spacer + ".blastnOUT")
                spacer_contig = CRISPR_Cas["Contig"]
                spacer_location = "-".join(
                    [df_Crispr_Table[df_Crispr_Table["CRISPR"] == spacer]["Start"].to_string(index=False),
                     df_Crispr_Table[df_Crispr_Table["CRISPR"] == spacer]["End"].to_string(index=False)])
                CRISPR_spacer_with_cas_locations=str(min([int(q) for q in spacer_location.split("-") + location.split("-")]))+"-"+str(max([int(p) for p in spacer_location.split("-") + location.split("-")]))

                subprocess.Popen(
                    ["blastn", "-query", spacer_fna, "-db", fna_blastdb, "-out", spacer_fna_blastOutfile,
                     "-num_threads",
                     str(threads), "-outfmt", "6"]).wait()
                df_blastOUT = pd.read_csv(spacer_fna_blastOutfile, header=None, sep="\t")
                for index, row in df_blastOUT.iterrows():
                    blast_location = str(row[8]) + "-" + str(row[9])
                    target_contig = str(row[1])
                    if distance_cal(CRISPR_spacer_with_cas_locations, blast_location) > 5000 and target_contig == spacer_contig:
                        # print(spacer + ":" + spacer_location + "=>" + target_contig + ":" + blast_location)
                        self_targeting_regions.append(
                            spacer + ":" + spacer_location + "=>" + target_contig + ":" + blast_location)
                    elif target_contig != spacer_contig:
                        # print(spacer + ":" + spacer_location + "=>" + target_contig + ":" + blast_location)
                        self_targeting_regions.append(spacer + ":" + spacer_location + "=>" + target_contig + ":" + blast_location)
            if len(self_targeting_regions) > 0:
                CC_list.append(CRISPR_Cas["Contig"] + "|" + CRCasType + "|" + location + "|" + "STSS=" + "+".join(
                    self_targeting_regions)) #This location indicate the Cas operon
                STSS_region.append("+".join(self_targeting_regions))
            else:
                CC_list.append(CRISPR_Cas["Contig"] + "|" + CRCasType + "|" + location + "|" + "No_STSS")
                STSS_region.append("No_STSS")
            # Contig|CasTyper|CasPosition|STSS_info
        df_CCtable["Contig"] = CC_contig
        df_CCtable["CRISPR-Cas Type"] = CC_type
        # df_CCtable["CRISPR-Cas Location"] = CC_location
        df_CCtable["CRISPR operon and location"] = C_operon_withLocation
        df_CCtable["Cas operon and location"] = Cas_operon_withLocation
        df_CCtable["STSS"] = STSS_region
        df_CCtable["Contig Length"] = contig_length
        df_CCtable.to_csv(os.path.join(general_folder, "CRISPR-Cas_found.csv"), index=False)
        return CC_list


def find_prophage(fna,outputdir,threads):
    subprocess.Popen(["VIBRANT_run.py","-i",fna,"-folder",outputdir,"-t",threads,"-no_plot"]).wait()
    phage_combined=subprocess.run(["find %s -name *phages_combined.faa"%outputdir],stdout=subprocess.PIPE,shell=True).stdout.decode('utf-8').rstrip()
    if os.path.isfile(phage_combined) and os.path.getsize(phage_combined) > 0:
        print("All prophage found stored in file %s"%phage_combined)
        phage_pos_start_dic = {}
        phage_pos_end_dic = {}
        for line in subprocess.Popen(["grep", ">", phage_combined], stdout=subprocess.PIPE).stdout:
            line = line.decode('utf-8').rstrip().split()
            contig = line[0].lstrip(">") # this is actually proID
            if "fragment" in contig:
                fragment=[contig.split("fragment_")[-1].split("_")[0]]
            else: fragment=[]
            positions_start_end = [int(s) for s in
                                   [v.lstrip("(").rstrip(")") for v in line if ".." in v][0].split("..")]
            if len(fragment) > 0:
                contig_fragment_key = contig.split("_fragment")[0] + "|" + fragment[0]
                if contig_fragment_key not in phage_pos_start_dic:
                    phage_pos_start_dic.setdefault(contig_fragment_key, [positions_start_end[0]])
                else:
                    phage_pos_start_dic[contig_fragment_key].append(positions_start_end[0])
                if contig_fragment_key not in phage_pos_end_dic:
                    phage_pos_end_dic.setdefault(contig_fragment_key, [positions_start_end[1]])
                else:
                    phage_pos_end_dic[contig_fragment_key].append(positions_start_end[1])
            else:
                contig_key = contig.rsplit("_",1)[0]
                if contig_key not in phage_pos_start_dic:
                    phage_pos_start_dic.setdefault(contig_key, [positions_start_end[0]])
                else:
                    phage_pos_start_dic[contig_key].append(positions_start_end[0])
                if contig_key not in phage_pos_end_dic:
                    phage_pos_end_dic.setdefault(contig_key, [positions_start_end[1]])
                else:
                    phage_pos_end_dic[contig_key].append(positions_start_end[1])
        phage_locations=[]
        prophage_out_table = os.path.join(outputdir, "prophage_locations.csv")
        newfile = open(prophage_out_table, "w")
        newfile.write(",".join(["Contig","Start","End","Contig Length"])+"\n")
        fna_dic=SeqIO.to_dict(SeqIO.parse(fna,"fasta"))
        for key in phage_pos_start_dic:
            pos_end = max(phage_pos_start_dic[key])
            pos_start = min(phage_pos_end_dic[key])
            phage_locations.append(key.split("|")[0] + ":" + str(pos_start) + "-" + str(pos_end)) # Contig:startPos-endPos
            newfile.write(",".join([key.split("|")[0],str(pos_start),str(pos_end),str(len(fna_dic[key.split("|")[0]].seq))])+"\n")
        return phage_locations
    else:
        print("No prophage found in sequence")
        return  None

def prophage_harboring_operonFind(df_outTable,prophage_regions):
    operon_location = df_outTable.head(1)["Start"].to_string(index=False) + "-" + df_outTable.tail(1)["End"].to_string(index=False)
    operon_contig = df_outTable.head(1)["Contig ID"].to_string(index=False)
    prophage_withOperon=[v for v in prophage_regions if distance_cal(v.split(":")[-1],operon_location) == 0 and operon_contig == v.split(":")[0]]
    if len(prophage_withOperon)>0:
        return ";".join(prophage_withOperon) # Contig:startPos-endPos
    else:
        return None

def pfamScan_run(operon_faa_file, pfam_hmm_dir,threads,HTH_alignment_evalue):
    pfamOut_file=operon_faa_file+".pfamScanOut"
    subprocess.Popen(["pfam_scan.pl","-fasta",operon_faa_file,"-dir",pfam_hmm_dir,"-outfile",pfamOut_file,"-e_seq",HTH_alignment_evalue,"-cpu",threads]).wait()
    dic_pfam={}
    for record in SeqIO.parse(operon_faa_file,"fasta"):
        info_list=[]
        for line in subprocess.Popen(["grep", record.id,pfamOut_file],stdout=subprocess.PIPE).stdout:
            line=line.decode('utf-8').rstrip().split()
            if len(line)>1:
                info_list.append("|".join([line[5],line[6],line[12],line[14]])) #hmmacc|hmmname|evalue|clan
        if len(info_list)>0:
            dic_pfam.setdefault(record.id,";".join(info_list))
    return dic_pfam

def Aca_HMM_search(aca_candidate_file,published_acaHMM,threads,hmm_outfile,evalue_cut_off,outdir,coverage_cutoff,gff_file,fna_file):
    fna_dic=SeqIO.to_dict(SeqIO.parse(fna_file,"fasta"))
    subprocess.Popen(
        ['hmmsearch', '--domtblout', hmm_outfile, '-o', os.path.join(outdir, 'log_aca.hmm'), '--noali', "--cpu", threads,
         '-E', evalue_cut_off, published_acaHMM, aca_candidate_file]).wait() # add coverage
    AcaPub_Pro_lst_dic={}
    df=pd.DataFrame()
    contig_list=[]
    wpid_list=[]
    start_list=[]
    end_list=[]
    strand_list=[]
    contig_length_list=[]
    hmm_ID_list = []
    hmm_coverage_list=[]
    evalue_list=[]
    for line in open(hmm_outfile).readlines():
        line = line.split()
        if "#" not in line[0]:
            ID = line[0]
            aca = line[3]
            coverage=(int(line[16]) - int(line[15]) + 1) / int(line[5]) #hmm coverage
            if coverage > coverage_cutoff: #Coverage filter
                info_from_gff=subprocess.run(["grep CDS %s|grep %s"%(gff_file,ID)],shell=True,stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip().split("\t")
                contig_list.append(info_from_gff[0])
                wpid_list.append(ID)
                start_list.append(info_from_gff[3])
                end_list.append(info_from_gff[4])
                strand_list.append(info_from_gff[6])
                contig_length_list.append(len(fna_dic[info_from_gff[0]].seq))
                hmm_ID_list.append(aca)
                hmm_coverage_list.append(str(coverage))
                evalue_list.append(line[6])

                if ID not in AcaPub_Pro_lst_dic:
                    AcaPub_Pro_lst_dic.setdefault(ID, [aca])
                elif ID in AcaPub_Pro_lst_dic:
                    AcaPub_Pro_lst_dic[ID].append(aca)
    if len(AcaPub_Pro_lst_dic) > 0:
        aca_fasta_file = open(os.path.join(outdir, "Aca-like_protein.faa"), "w")
        faa_dic = SeqIO.to_dict(SeqIO.parse(aca_candidate_file, "fasta"))
        for key in AcaPub_Pro_lst_dic.keys():
            SeqIO.write(faa_dic[key],aca_fasta_file,"fasta")
        df["Contig"]=contig_list
        df["Protein ID"]=wpid_list
        df["Start Location"]= start_list
        df["End Location"]= end_list
        df["Strand"] = strand_list
        df["Contig Length (nt)"]=contig_length_list
        df["Aca HMM ID"]=hmm_ID_list
        df["Aca-like Protein Coverage"]=hmm_coverage_list
        df["Aca HMM Evalue"]=evalue_list
        df.to_csv(os.path.join(outdir, 'Aca-like_protein.csv'),index=False)

    return AcaPub_Pro_lst_dic
