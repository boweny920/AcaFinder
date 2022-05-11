
# The gff to be used must be GFF3, format explaination: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

import argparse
from Annotation import annotation_prodigal
import os
import sys

def is_non_zero_file(fpath):
    #Check if file is empty or does not exsit
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def make_output_folder(GCF_Name):
    os.makedirs(GCF_Name)
    return os.path.dirname(GCF_Name)

parser = argparse.ArgumentParser(description='This script will find Acr homologs with diamond, Find potential Acr-Aca loci and find perform a hmmscan with the 89 published HTH domain againest those loci. '
                                             'Output of the script: Acr homologs fasta file; Potential Acr-Aca loci file; Candidate Aca fasta file; Diamond output file hmmscan output file'
                                             'look for file "!!!!! Potential Acr-Aca loci found and saved in xxx"; "$$$$$ Aca candidate found and saved in file xxx" for final output results')
parser.add_argument('-n','--FNA_file',nargs='?',default=None,help="Genome FNA file")
parser.add_argument('-p','--FAA_file', nargs='?',default=None ,help="Genome annotated FAA file")
parser.add_argument('-g','--GFF_file', nargs='?',default=None, help="Genome annotated GFF file")
parser.add_argument('-m','--mode_prodiagal', nargs='?', default="meta",choices=["single","meta"],help="mode prodigal will be run")
parser.add_argument('-o','--outputFolder',type=str,default="AcaFinder_Output",help="folder containing all output results of all ran GCFs")
parser.add_argument('-a','--Acr_alignment_evalue', nargs='?', default="1e-3",help="evalue cut-off for Acr homolog search")
parser.add_argument('-c','--Acr_alignment_coverage', nargs='?', default="0.6",help="coverage cut-off for Acr homolog search")
parser.add_argument('-t','--HTH_alignment_evalue', nargs='?', default="1e-3",help="evalue cut-off for HTH domian hummer search")
parser.add_argument('-v','--HTH_alignment_coverage', nargs='?', default="0.5",help="coverage cut-off for HTH domian hummer search")
parser.add_argument('-l','--all_protein_length_in_AcrAca_operon', nargs='?', default=600,type=int,help="max proten lenght in Acr-Aca operon")
parser.add_argument('-i','--intergenic_dist_in_AcrAca_operon',nargs='?',default=250,type=int,help="Maximum Intergenic distance in Acr-Aca operon")
parser.add_argument('-r','--Acr_protein_database',nargs='?',default="AcrDatabase.faa",help="The Acr proteins that will be used search for Acas, default are the published Acrs + AcrHub predicted Acrs + 2500 high confident Acr prediction of AcrCatalog")
parser.add_argument('-e','--HTH_hmm',nargs='?',default="HMM/Pfam-A.clan_HTH_le150bp.hmm",help="The hmmdb build from the HTH hmms used to search Acas, default are the Pfam hmms from clan_HTH with model length < 150bp")
parser.add_argument('-b','--Acr_Aca_inBetweenGenes', nargs='?', default=4,help="Maximum number of genes allowed between Aca and Acr proteins + 1 (e.g if the input is 4, then maximum 3 genes are allowed between the potental Aca genes to its closest Acr homolog)")
parser.add_argument('-w', '--Virus', action='store_true',help="Provide option -w/--Virus if input data is of viral origin")
parser.add_argument('-d','--threads',nargs='?', default="1", help="Number of cpu cores to run the program with")
parser.add_argument('-z','--phamDir', nargs='?', default="all_pFam_hmm", help="Directory of all pfam hmm files with .dat files and other binaries")
parser.add_argument('-y','--published_acaHMM', nargs='?', default="HMM/published_13_Aca.hmm",help="HMM for all 13 publsihed Aca proteins, recommended to use the default hmm provided from us")
parser.add_argument('-x','--acaHMM_evalue', nargs='?', default="1e-10",help="evalue cut-off for acaHMM hummer search, recommended to use default")
parser.add_argument('-u','--acaHMM_cov', nargs='?', default=0.9,type=float,help="Coverage cut-off for acaHMM hummer search, recommended to use default")
args=parser.parse_args()

if os.path.isdir(args.outputFolder) is not True:
    os.makedirs(args.outputFolder)

if (args.FAA_file == None or args.GFF_file == None) and args.FNA_file != None:
    print("No annotations provided, using prodigal to annotate genome :)")
    prodigal_outDir=os.path.join(args.outputFolder,os.path.basename(args.FNA_file)+".prodigalOUT")
    gff,faa = annotation_prodigal(args.FNA_file, args.mode_prodiagal,prodigal_outDir).run_prodigal()
    fna=args.FNA_file
    from AcaFind_process_verProdigal import Aca_Find_process

elif args.FAA_file != None and args.GFF_file != None and args.FNA_file != None:
    gff=args.GFF_file
    faa=args.FAA_file
    fna=args.FNA_file
    from AcaFind_process import Aca_Find_process

elif (args.FAA_file is None or args.GFF_file is None) and args.FNA_file == None:
    sys.exit("##No FNA, GFF and FAA detected##\nPlease provide FNA file or an annotated protein FAA file with an associated GFF file of you sequence of interest.")

Aca_Find_process(fna,gff,faa,
                 args.outputFolder, args.Acr_alignment_evalue, args.Acr_alignment_coverage,
                 args.HTH_alignment_evalue, args.HTH_alignment_coverage, args.all_protein_length_in_AcrAca_operon,
                 args.intergenic_dist_in_AcrAca_operon,args.Acr_protein_database,args.HTH_hmm,args.Acr_Aca_inBetweenGenes,
                 args.Virus, args.threads,args.phamDir,args.published_acaHMM,args.acaHMM_evalue,args.acaHMM_cov).run_process()
