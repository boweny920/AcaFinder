from functions import *
import os

class Aca_Find_process:

    def __init__(self,gff_file,faa_file,outputFolder,Acr_alignment_evalue,Acr_alignment_coverage,HTH_alignment_evalue,HTH_alignment_coverage,
                 all_protein_length_in_AcrAca_operon,intergenic_dist_in_AcrAca_operon,KnowAcrFaa,HTHhmm,Acr_Aca_inBetweenGenes):
        self.faa=faa_file
        self.gff=gff_file
        self.outputFolder=outputFolder
        self.Acr_alignment_evalue=Acr_alignment_evalue
        self.Acr_alignment_coverage = Acr_alignment_coverage
        self.HTH_alignment_evalue=HTH_alignment_evalue
        self.HTH_alignment_coverage = HTH_alignment_coverage
        self.all_protein_length_in_AcrAca_operon=all_protein_length_in_AcrAca_operon
        self.intergenic_dist_in_AcrAca_operon=intergenic_dist_in_AcrAca_operon
        self.KnowAcrFaa= KnowAcrFaa
        self.HTHhmm=HTHhmm
        self.Acr_Aca_inBetweenGenes=Acr_Aca_inBetweenGenes

    def run_process(self):
        sub_outputfolder_path=self.outputFolder
        file_list=[]
        file=subprocess.Popen(['grep "CDS" %s'%(self.gff)],shell=True, stdout=subprocess.PIPE)
        for line in file.stdout:
            line = line.decode('utf-8').rstrip().split('\t')
            file_list.append(line)
        file_dic={}
        for i in enumerate(file_list):
            file_dic.setdefault(i[0],i[1])

        loci_list,loci_list_result_check=loci_select_before_diamond(file_dic,self.all_protein_length_in_AcrAca_operon,self.intergenic_dist_in_AcrAca_operon)
        diamond_outfile=run_diamond(self.faa,sub_outputfolder_path,self.KnowAcrFaa,self.Acr_alignment_evalue,self.Acr_alignment_coverage) # include KnownAcr.faa in tool directory

        if is_non_zero_file(diamond_outfile) is not False:
            print("Acr homologs found in annotation genome","...")
            Acr_homolog_candidate_list = []
            # For result checking
            Acr_homolog_candidate_list_resultCheck = []
            protein_NP_list, dic_Acr = parse_diamond_get_proteinID(diamond_outfile) #Need to add coverage too
            ##Write the Acr homologs to a fasta file
            faa_file_wrote(protein_NP_list, self.faa, os.path.join(sub_outputfolder_path, "Acr_homologs.faa"))
            for i in loci_list:
                # Note: v[1] would be the protein ID not the gene order
                if any(v for v in i if v[1] in protein_NP_list): Acr_homolog_candidate_list.append(i)
                ##To see if any sub-list of loci_list has values in protein_NP_list.

            ### For result checking, v[0] is the Protein ID
            for i in loci_list_result_check:
                if any(v for v in i if v[0] in protein_NP_list):
                    for v in i:
                        # Check which ProteinID is the Acr homolog
                        if v[0] in protein_NP_list:
                            v.append(dic_Acr[v[0]])
                        else:
                            v.append("NA")
                    Acr_homolog_candidate_list_resultCheck.append(i)

            if len(Acr_homolog_candidate_list) > 0:
                print("Acr homologs found in identified short gene operons","...")
                # Put the new faa files in a directory called "Acr_homolog_positive_short_gene_operons"
                output_1_DirPath = os.path.join(sub_outputfolder_path,"Acr_homolog_positive_Short_Gene_Operons")
                os.makedirs(output_1_DirPath)
                # Write each locus with Acr homolog into individual fasta files
                for order, locus in enumerate(Acr_homolog_candidate_list):
                    potential_Acr_Aca_filename = os.path.join(output_1_DirPath,"Acr_Homolog_poisitve_SGO_OperonNumber-" + str(order) + ".faa")
                    potential_Aca_filename = os.path.join(output_1_DirPath,"Aca_candidates_within_Acr_Homolog_poisitve_SGO_OperonNumber-" + str(order) + ".faa") # This file contains non-Acr proteins in the SGO
                    unhmmer_aca = aca_select_process(locus, self.faa, potential_Aca_filename, protein_NP_list,self.Acr_Aca_inBetweenGenes)
                    faa_file_wrote_potential_AcrAca(locus, self.faa, potential_Acr_Aca_filename)
                    print("Acr homolog positive short gene operon protein sequences saved in %s" % potential_Acr_Aca_filename)

                    ##run hmmscan on this "potential_Aca_filename"
                    hmm_outfile = run_hmmscan(potential_Acr_Aca_filename,self.HTH_alignment_evalue,self.HTHhmm,self.outputFolder)
                    if is_non_zero_file(hmm_outfile) is not False:
                        ##parse the hmm_outfile based on coverage>0.5
                        parsed_hmm_outfile = parse_hmmOutfile(hmm_outfile,self.HTH_alignment_coverage) # Need to correct coverage, the coverage in this function is not correct
                        if is_non_zero_file(parsed_hmm_outfile) is not False:
                            aca_candidate_file, Aca_Pro_lst = potential_new_aca_faa_filemake(parsed_hmm_outfile, self.faa,
                                                                                             unhmmer_aca)
                            if is_non_zero_file(aca_candidate_file) is not False:
                                print("Potential Aca genes found and saved in file %s" % aca_candidate_file)

                            ####For result checking#####
                            Acr_Aca_loci_checkResult = os.path.join(output_1_DirPath, "GBA_identified_AcrAca_loci_OperonNumber-" + str(
                                order) + ".check_Result")
                            final_result_check_output_generation(Acr_homolog_candidate_list_resultCheck,
                                                                 Acr_Aca_loci_checkResult, Aca_Pro_lst)
                            #####For result checking#####
            elif len(Acr_homolog_candidate_list) == 0:
                print("%s have %s number of identified Acr homologs, but not in the identified short gene operons, "
                      "maybe try and adjust the 'Short Gene Operon Identification' parameters" % (self.faa, len(protein_NP_list)))
        else:
            print("No Acr homologs in annotated protein file '%s' under the current set condition, try loosen it a little :)" %self.faa)
