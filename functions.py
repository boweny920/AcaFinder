# Need to make a seperate function file for prodigal runs and ncbi annotation runs, with only a few functions, not the whole thing!

import subprocess
import sys
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
    ##modified to fit prodigal output
    return [list[0], list[-1]]

def proteinID(IDinfo_prodiagal):
    ##modified to fit prodigal output
    IDinfo = IDinfo_prodiagal[-1].split(";")
    for v in IDinfo:
        if "ID=" in v:
            ID = str(v).split("_",1)
            ProID=str(IDinfo_prodiagal[0])+"_"+ID[-1]
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

def run_diamond(database_file,outputdir,query,evalue_cutoff,coverage_cut_off): #add coverage cutoff
    ##Use this as dmond_out
    subprocess.Popen(['diamond','makedb','--in',database_file,'-d',os.path.join(outputdir,"diamond_db")]).wait()
    subprocess.Popen(['diamond','blastp','-d',os.path.join(outputdir,"diamond_db"),'-f','6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen' ,'qstart' ,'qend' ,'sstart' ,'send' ,'evalue' ,'bitscore', 'qlen', 'slen','--quiet','--more-sensitive','-q',query,'-p','5','-o',os.path.join(outputdir,"diamond_blastp_result.txt"),'-e',evalue_cutoff]).wait()
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

def run_hmmscan(faafile,evalue_cut_off,hmmfile,outdir):
    hmm_outfile=faafile+".hmmout"
    subprocess.Popen(['hmmscan','--domtblout',hmm_outfile,'-o',os.path.join(outdir,'log.hmm'),'--noali','-E', evalue_cut_off,hmmfile,faafile]).wait()
    return hmm_outfile

def parse_hmmOutfile(hmm_outfile,hmm_coverage_cutoff): # redo the coverage calculation part!
    parsed_hmm_outfile = str(hmm_outfile) + ".Coverage_parsed"
    #Parse the hmm output in terms of subject coverage: subjust aligned length/ Total subject length
    subprocess.Popen(["grep -v '#' %s | awk '($17-$16)/$3 > %s {print}'  > %s" % (hmm_outfile, str(hmm_coverage_cutoff),parsed_hmm_outfile)],shell=True).wait()
    return parsed_hmm_outfile

def potential_new_aca_faa_filemake(parsed_hmm_outfile,newAcrAca_faaFile,unhmmer_aca):
    ## Make a fasta file of potential Aca amino acid sequences
    ## Make a list of Aca, for checking results
    Aca_Pro_lst=[]
    record_dict = SeqIO.to_dict(SeqIO.parse(newAcrAca_faaFile, "fasta"))
    newfile_name= parsed_hmm_outfile+".new_found_ACA.faa"
    with open(newfile_name,"w") as new:
        acaout=subprocess.Popen(["awk '{print $4}' %s |sort -u "%parsed_hmm_outfile],shell=True, stdout=subprocess.PIPE)
        for ID in acaout.stdout:
            ID=ID.strip().decode('utf-8')
            if ID in unhmmer_aca:
                Aca_Pro_lst.append(ID)
                # If Aca is not in unhmmer, then it should not be considered as a potential Aca(here is to make sure that condiction "gene distance < 3" is applied)
                SeqIO.write(record_dict[ID], new, "fasta")
                #This is to consider the < than 3 genes condition, which had been implemented before in "faa_file_wrote_diamond_special" function
    return newfile_name,Aca_Pro_lst

def final_result_check_output_generation(Acr_homolog_candidate_list_resultCheck,newfile_name,Aca_Pro_lst):
    #generation of final result checking file
    #Table will be tsv format
    #### GCF_number    ProteinID   Contig  Strand  Length  Start   End     Species     Acr     Aca ####
    Acr_Aca_loci_by_GBA = []
    for i in Acr_homolog_candidate_list_resultCheck:
        if any(v for v in i if v[0] in Aca_Pro_lst):
            for v in i:
                # Check which ProteinID is the Aca homolog
                if v[0] in Aca_Pro_lst:
                    v.append("Aca_protein")
                else:
                    v.append("NA")
            Acr_Aca_loci_by_GBA.append(i)

    with open(newfile_name, "w") as newfile:
        for pro in Acr_Aca_loci_by_GBA:
            for value in pro:
                for info in value: newfile.write("%s\t" % str(info))
                newfile.write("\n")