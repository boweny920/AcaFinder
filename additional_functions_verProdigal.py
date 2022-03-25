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
########################
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

def faa_file_wrote(one_acr_aca_locus,faa_file,newfile_name_dirctory):
    ##The one_acr_aca_locus is a list of proteinID.
    ##The function takes the IDs and the faa file with protein sequences
    ##And make a fasta newfile containing ProteinIDs with corresponding amino acid sequences
    record_dict = SeqIO.to_dict(SeqIO.parse(faa_file, "fasta"))
    with open(newfile_name_dirctory,"w") as newfile:
        for proID in one_acr_aca_locus:
            SeqIO.write(record_dict[proID],newfile,"fasta")
