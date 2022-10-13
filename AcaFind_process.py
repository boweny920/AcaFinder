import os
import pandas as pd
from functions import *
from Bio import SeqIO
import numpy as np

class Aca_Find_process:

    def __init__(self,fna_file,gff_file,faa_file,outputFolder,Acr_alignment_evalue,Acr_alignment_coverage,HTH_alignment_evalue,HTH_alignment_coverage,
                 all_protein_length_in_AcrAca_operon,intergenic_dist_in_AcrAca_operon,KnowAcrFaa,HTHhmm,Acr_Aca_inBetweenGenes,Virus_Flag,threads,phamDir,published_acaHMM,acaHMM_evalue,acaHMM_cov):
        self.fna=fna_file
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
        self.Viral_Flag=Virus_Flag
        self.threads=threads
        self.phamDir=phamDir
        self.published_acaHMM=published_acaHMM
        self.acaHMM_evalue=acaHMM_evalue
        self.acaHMM_cov=acaHMM_cov

    def run_process(self):
        sub_outputfolder_path=self.outputFolder
        file_list=[]
        file=subprocess.Popen(['grep "CDS" %s|grep -v "#"'%(self.gff)],shell=True, stdout=subprocess.PIPE)
        for line in file.stdout:
            line = line.decode('utf-8').rstrip().split('\t')
            file_list.append(line)
        file_dic={}
        for i in enumerate(file_list):
            file_dic.setdefault(i[0],i[1])

        diamond_outfile=run_diamond(self.faa,sub_outputfolder_path,self.KnowAcrFaa,self.Acr_alignment_evalue,self.Acr_alignment_coverage,self.threads) # include KnownAcr.faa in tool directory
        #Search for prophage and CRISPR-Cas
        if self.Viral_Flag is False:
            ##check if any complete CRISPR-Cas found, and check for self-targeting locations if any
            complete_CRISPR_Cas_systems = find_complete_CRISPR_Cas_and_SelfTargeting(self.fna, os.path.join(sub_outputfolder_path, "CRISPR_Cas_Found"), self.threads, sub_outputfolder_path)  # Name of cctyper output directory is "CRISPR_Cas_Found"
            # [Contig|CasTyper|Position|self-targeting regions,...]/None
            prophage_regions = find_prophage(self.fna, sub_outputfolder_path, self.threads)
            # [Contig:startPos-endPos,...]/None
        else:
            complete_CRISPR_Cas_systems = None
            prophage_regions = None
        # Search the entire genome for potetnial Aca hmm hits
        publishedAcaHMM_hits_dic = Aca_HMM_search(self.faa, self.published_acaHMM, self.threads,os.path.join(sub_outputfolder_path, "Aca_HMM_hits.hmmOut"),self.acaHMM_evalue,sub_outputfolder_path,self.acaHMM_cov,self.gff,self.fna,complete_CRISPR_Cas_systems,prophage_regions)

        if is_non_zero_file(diamond_outfile) is not False:
            print("Acr homologs found in annotation genome","...")
            Acr_homolog_candidate_list = []
            # For result checking
            Acr_homolog_candidate_list_resultCheck = []
            protein_NP_list, dic_Acr = parse_diamond_get_proteinID(diamond_outfile) #Need to add coverage too
            ##Write the Acr homologs to a fasta file
            faa_file_wrote(protein_NP_list, self.faa, os.path.join(sub_outputfolder_path, "Acr_homologs.faa"))
            loci_list, loci_list_result_check = loci_select(file_dic, dic_Acr, self.KnowAcrFaa, self.all_protein_length_in_AcrAca_operon, self.intergenic_dist_in_AcrAca_operon)
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
                            v.append(dic_Acr[v[0]].replace(',',';'))
                        else:
                            v.append("NA")
                    Acr_homolog_candidate_list_resultCheck.append(i)

            if len(Acr_homolog_candidate_list) > 0:
                print("Acr homologs found in identified short gene operons","...")
                # Put the new faa files in a directory called "Acr_homolog_positive_short_gene_operons"
                output_1_DirPath = os.path.join(sub_outputfolder_path,"Acr_homolog_positive_Short_Gene_Operons")
                os.makedirs(output_1_DirPath)

                # Write each locus with Acr homolog into individual fasta files
                output_checkResult_tables=[]
                for order, locus in enumerate(Acr_homolog_candidate_list):
                    potential_Acr_Aca_filename = os.path.join(output_1_DirPath,"Acr_Homolog_poisitve_SGO_OperonNumber-" + str(order) + ".faa")
                    potential_Aca_filename = os.path.join(output_1_DirPath,"Aca_candidates_within_Acr_Homolog_poisitve_SGO_OperonNumber-" + str(order) + ".faa") # This file contains non-Acr proteins in the SGO
                    unhmmer_aca = aca_select_process(locus, self.faa, potential_Aca_filename, protein_NP_list,self.Acr_Aca_inBetweenGenes)
                    faa_file_wrote_potential_AcrAca(locus, self.faa, potential_Acr_Aca_filename)
                    print("Acr homolog positive short gene operon protein sequences saved in %s" % potential_Acr_Aca_filename)

                    ##run hmmscan on this "potential_Aca_filename"
                    hmm_outfile = run_hmmscan(potential_Acr_Aca_filename,self.HTH_alignment_evalue,self.HTHhmm,self.outputFolder,self.threads)
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
                                                                 Acr_Aca_loci_checkResult, Aca_Pro_lst,publishedAcaHMM_hits_dic)
                            output_checkResult_tables.append(Acr_Aca_loci_checkResult)
                            #####For result checking#####

                ###CRISPR-Cas PRophage finding, and final output table make###
                if len(output_checkResult_tables) > 0:

                    protein_faa_dic=SeqIO.to_dict(SeqIO.parse(self.faa,"fasta"))
                    fna_dic=SeqIO.to_dict(SeqIO.parse(self.fna,"fasta"))
                    df_allResult=pd.DataFrame(columns=["Operon Number","Protein ID","Contig ID","Strand","Protein Length (nt)","Start","End","Acr Homolog","Potential Aca","AcaHMM HIT","Pfam","Complete CRISPR-Cas and STSS", "Operon in Prophage","Protein Sequence"])
                    for out_table in output_checkResult_tables:
                        if is_non_zero_file(out_table):
                            faa_list = []
                            operon_number_list = []
                            prophage_withOperon_location_list=[]
                            Complete_CRISPR_Cas_inContig_list = []
                            pfam_list=[]
                            contig_length_list=[]
                            df_outTable = pd.read_csv(out_table,header=None,sep="\t")
                            # df_outTable.drop(9, axis=1, inplace=True) #Need to check, need to drop still?
                            df_outTable.columns=["Protein ID","Contig ID","Strand","Protein Length (nt)","Start","End","Acr Homolog","Potential Aca","AcaHMM HIT"]

                            ##check if aca operon in prophage
                            if prophage_regions is not None:
                                prophage_containing_AcaOperon = prophage_harboring_operonFind(df_outTable, prophage_regions) # Contig:startPos-endPos/None
                            else: prophage_containing_AcaOperon = None

                            #Run pfamscan on the operon faa file
                            operon_faa_file = os.path.join(output_1_DirPath, "Acr_Homolog_poisitve_SGO_OperonNumber-" + out_table.split("-")[-1].rstrip(".check_Result") + ".faa")
                            dic_pfam = pfamScan_run(operon_faa_file,self.phamDir,self.threads,self.HTH_alignment_evalue)

                            for index, row in df_outTable.iterrows():
                                if "pseudo" in row["Protein ID"]:
                                    faa_list.append("pseudo")
                                else:
                                    faa_list.append(str(protein_faa_dic[row["Protein ID"]].seq))
                                operon_number_list.append(out_table.rstrip(".check_Result").split("_")[-1])
                                if prophage_containing_AcaOperon is not None: prophage_withOperon_location_list.append(prophage_containing_AcaOperon)
                                else: prophage_withOperon_location_list.append(np.nan)
                                if complete_CRISPR_Cas_systems is not None: Complete_CRISPR_Cas_inContig_list.append(";".join(complete_CRISPR_Cas_systems))
                                else: Complete_CRISPR_Cas_inContig_list.append(np.nan)
                                if row["Protein ID"] in dic_pfam and pd.isna(row["Potential Aca"])==True:
                                    pfam_list.append(dic_pfam[row["Protein ID"]])
                                else: pfam_list.append(np.nan)
                                contig_length_list.append(len(fna_dic[row["Contig ID"]].seq))
                            df_outTable["Operon Number"]=operon_number_list
                            df_outTable["Complete CRISPR-Cas and STSS"] = Complete_CRISPR_Cas_inContig_list
                            df_outTable["Operon in Prophage"] = prophage_withOperon_location_list
                            df_outTable["Protein Sequence"]=faa_list
                            df_outTable["Pfam"]=pfam_list
                            df_outTable["Contig Length (nt)"]=contig_length_list
                            df_outTable = df_outTable[["Operon Number","Protein ID","Contig ID","Contig Length (nt)","Strand","Protein Length (nt)","Start","End","Acr Homolog","Potential Aca","AcaHMM HIT","Pfam","Complete CRISPR-Cas and STSS", "Operon in Prophage","Protein Sequence"]]
                            df_allResult=pd.concat([df_allResult,df_outTable],ignore_index=True)
                        df_allResult_withIRana=IR_find(df_allResult,self.fna,sub_outputfolder_path)
                        df_allResult_withIRana.to_csv(os.path.join(sub_outputfolder_path,"All_Aca_operons.csv"),index=False)

                else:
                    print("%s have %s number of identified Acr homologs, but not in the identified short gene operons, but no Aca candidate"% (self.faa, len(protein_NP_list)))

            elif len(Acr_homolog_candidate_list) == 0:
                print("%s have %s number of identified Acr homologs, but not in the identified short gene operons, "
                      "maybe try and adjust the 'Short Gene Operon Identification' parameters" % (self.faa, len(protein_NP_list)))
        else:
            print("No Acr homologs in annotated protein file '%s' under the current set condition :)" %self.faa)
