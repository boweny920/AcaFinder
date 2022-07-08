<img src="http://bcb.unl.edu/AcaFinder/img/AcaFinder-logos.jpeg" width="400" height="400">

<center>(c) <a href='http://bcb.unl.edu'>Yin Lab</a>@<a href='https://www.unl.edu'>UNL</a>2022</center>


## Contents:

<a href='#installation'>I. Installation / Dependencies</a>

<a href='#about'>II. About</a>

<a href='#using_acafinder'>III. Using AcaFinder</a>


<a href='#examples'>V. Examples</a>

<a href='#workflow'>VI. Workflow</a>

<a href='#faq'>VII. FAQ</a>

****

<div id='installation' />

## I. Installation / Dependencies

### Dependences 

Program expects these versions and using other versions can result in unexpected behavior.

`Python` - Python version should be >=3.8. It is recommended that you create a conda enviroment and install all the below dependences within said enviroment:
```sh
## python version 3.9 is used here as an example
conda create -n TestAcaFinder python=3.9
```

`VIBRANT` - Used to search for potential prophage regions from input genomic sequences

Version used v1.2.0. We recommend installing VIBRANT using Anaconda, but you may also install VIBRANT with other methods from https://github.com/AnantharamanLab/VIBRANT

To install using Anaconda, 
Install dependencies. See Requirements section https://github.com/AnantharamanLab/VIBRANT.
Install directly to $PATH using bioconda. 
```sh
conda install -c bioconda vibrant
```

Download and setup databases. This will take some time due to file sizes, but it only needs to be run once. This step requires ~20GB of temporary storage space and ~11GB of final storage space. To do this, run download-db.sh which should be in your system's $PATH. download-db.sh
```sh
download-db.sh
```

`Cctyper` - Used for complete CRISPR-Cas system search

Can be installed either through conda or pip.
It is advised to use conda, since this installs CRISPRCasTyper and all dependencies, and downloads the database in one go.

To use conda:
```sh
conda install -c conda-forge -c bioconda -c russel88 cctyper
```

To use pip:
python -m pip install cctyper

When installing with pip, you need to download the database manually:

```sh
# Download and unpack
svn checkout https://github.com/Russel88/CRISPRCasTyper/trunk/data
tar -xvzf data/Profiles.tar.gz
mv Profiles/ data/
rm data/Profiles.tar.gz

# Tell CRISPRCasTyper where the data is:
# either by setting an environment variable (has to be done for each terminal session, or added to .bashrc):
export CCTYPER_DB="/path/to/data/"
# or by using the --db argument each time you run CRISPRCasTyper:
cctyper input.fa output --db /path/to/data/
```

Detailed information can be found at : https://github.com/Russel88/CRISPRCasTyper#install

`PfamScan/Pfam Database` - For protein annotation 

Install PfamScan with Anaconda:
```sh
conda install -c bioconda pfam_scan
#or
conda install -c bioconda/label/cf201901 pfam_scan
```

Downloading Pfam database
```sh
# a. Go to the all_pFam_hmm/ folder

# b. Download database: 
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

#c. Unpack 
gunzip Pfam-A.hmm.gz

#d. Prepare HMM files 
hmmpress Pfam-A.hmm
```

****


<div id='about' />

## II. About 

### AcaFinder allows for automated genome mining for reliable Acas. 

To more confidently identify Acas given a genome or metagenome assembled genome, we implemented two approaches. The first approach is based on guilt-by-association (GBA), meaning that we identify homologs of Acrs first and then search for HTH-containing proteins in the acr gene neighborhood. The second approach is to build an HMM (hidden markov model) database using training data of the 12 known Aca families, and then search for Aca homologs with this Aca-HMMdb instead of Pfam HTH HMMs. In addition to the two implemented approaches, AcaFinder also integrates a CRISPR-Cas search tool (CRISPRCasTyper), a prophage search tool (VIBRANT), and in-house a Self-targeting spacer (STSS) searching tool, providing users with detailed information vital to the assessment of Aca predictions

****

<div id='using_acafinder' />

## **III. <span style='color:RebeccaPurple'>Using AcaFinder</span>**

#### Input 

AcrFinder needs **.fna**, **.gff** and **.faa** as input. Only **.fna** file as input is also acceptable; in that case, the **.gff** and **.faa** file will be generated by running <a href='https://github.com/hyattpd/Prodigal'>Prodigal</a>.

#### List of Options

| Option | Alternative | Purpose                     |
| -----  | ----------- | --------------------------- |
| -h     | --help      | Shows all available options |
| -n     | --FNA_file  | <span style="color:red">Required</span> fna file |
| -g     | --GFF_file     | <span style="color:red">Required</span> Path to gff file to use/parse |
| -p     | --FAA_file     | <span style="color:red">Required</span> Path to faa file to use/parse |
| -m     | --mode_prodiagal  | Mode prodigal will be run choices=["single","meta"], default=meta |
| -o     | --outputFolder      | Folder containing all output results, default=AcaFinder_Output |
| -a     | --Acr_alignment_evalue      | Evalue cut-off for Acr homolog search, default=1e-3 |
| -c     | --Acr_alignment_coverage      | Coverage cut-off for Acr homolog search, default=0.6 |
| -t     | --HTH_alignment_evalue    | Evalue cut-off for HTH domian hummer search, default=1e-3 |
| -v     | --HTH_alignment_coverage    | Coverage cut-off for HTH domian hummer search, default="0.4 |
| -l     | --all_protein_length_in_AcrAca_operon    | Max proten length in Acr-Aca operon when length of Acr homolog < 200aa, default=600 |
| -i     | --intergenic_dist_in_AcrAca_operon      | Maximum Intergenic distance in Acr-Aca operon, default=250 |
| -r     | --Acr_protein_database     | The Acr proteins that will be used search for Acas, default are the published Acrs + AcrHub predicted Acrs + 2500 high confident Acr prediction of AcrCatalog, default=AcrDatabase.faa |
| -e     | --HTH_hmm_strict     | Provide option -e/--HTH_hmm_strict to use the more strict HTH HMM database (HTH-HMM_strict) for potential Aca protein search in Aca-Acr operons |
| -b     | --Acr_Aca_inBetweenGenes     | Maximum number of genes allowed between Aca and Acr proteins + 1 (e.g if the input is 4, then maximum 3 genes are allowed between the potental Aca genes to its closest Acr homolog), default=4 |
| -w     | --Virus     | Provide option -w/--Virus if input data is of viral origin |
| -d     | --threads     | Number of cpu cores to run the program with, default=1 |
| -z     | --phamDir     | Directory of all pfam hmm files with .dat files and other binaries, default=AcaFinder/all_pFam_hmm | 
| -y     | --published_acaHMM     | HMM for all 13 publsihed Aca proteins, recommended to use the default hmm provided from us, default=AcaFinder/HMM/AcaHMMs |
| -x     | --acaHMM_evalue     | Evalue cut-off for acaHMM hummer search, recommended to use default, default=1e-10 |
| -u     | --acaHMM_cov     | Coverage cut-off for acaHMM hummer search, recommended to use default, default=0.6 |


#### Output files

| Name                 | Meaning               |
| -------------------- | --------------------- |
|*<output_dir>*/Acr_homolog_positive_Short_Gene_Operons   | Folder containing intermediate Acr-Aca operon result files |
|*<output_dir>*/CRISPR_Cas_Found   | CCtyper direct output folder |
|*<output_dir>*/VIBRANT_*<input_ID>*_genomic   | VIBRANT direct output folder |
|*<output_dir>*/Aca-like_protein.csv    | Final results from the AcaHMM search approach |
|*<output_dir>*/Aca-like_protein.faa    | Final results from the AcaHMM search approach, the Aca-like protein fasta sequences | 
|*<output_dir>*/Aca_HMM_hits.hmmOut    | AcaHMM hmmersearch output | 
|*<output_dir>*/log_acaHMM.hmm    | AcaHMM hmmersearch log file | 
|*<output_dir>*/Acr_homologs.faa    | Protein seuqnces of Acr homologs found from input genomic sequences | 
|*<output_dir>*/All_Aca_operons.csv    | Final results from the GBA search approach |
|*<output_dir>*/CRISPR-Cas_found.csv   | Summary of Complete CRISPR-Cas systems discovered |
|*<output_dir>*/diamond_blastp_result.txt   | Diamond blastp output from Acr homolog search |
|*<output_dir>*/diamond_blastp_result.coverageParsed.txt   | Diamond blastp output from Acr homolog search, parsed for coverage | 
|*<output_dir>*/log_HTH.hmm    | HTH domain hmmerscan log file | 
|*<output_dir>*/prophage_locations.csv    | Summary of prophage regions discovered |
|*<output_dir>*/Acr_homolog_positive_Short_Gene_Operons/Aca_candidates_within_Acr_Homolog_poisitve_SGO_OperonNumber-*<operon_ID>*.faa   | Protein fast file of final Aca-Acr operon | 
|*<output_dir>*/Acr_homolog_positive_Short_Gene_Operons/GBA_identified_AcrAca_loci_OperonNumber-*<operon_ID>*.check_Result  | Summary file of final Aca-Acr operon | 
|*<output_dir>*/Acr_homolog_positive_Short_Gene_Operons/Acr_Homolog_poisitve_SGO_OperonNumber-*<operon_ID>*.faa | Protein fast file of short-Gene-Operons with Acr homologs |
|*<output_dir>*/Acr_homolog_positive_Short_Gene_Operons/Acr_Homolog_poisitve_SGO_OperonNumber-*<operon_ID>*.faa.hmmout | hmmscan output of HTH search within short-Gene-Operons with Acr homologs |
|*<output_dir>*/Acr_homolog_positive_Short_Gene_Operons/Acr_Homolog_poisitve_SGO_OperonNumber-*<operon_ID>*.faa.hmmout.Coverage_parsed | hmmscan output of HTH search within short-Gene-Operons with Acr homologs |
|*<output_dir>*/Acr_homolog_positive_Short_Gene_Operons/Acr_Homolog_poisitve_SGO_OperonNumber-*<operon_ID>*.faa.hmmout.Coverage_parsed.new_found_ACA.faa | HTH positive proteins fasta sequences within short-Gene-Operons with Acr homologs |
*<output_dir>*/Acr_homolog_positive_Short_Gene_Operons/Acr_Homolog_poisitve_SGO_OperonNumber-*<operon_ID>*.faa.pfamScanOut | Pfam annoatations of proteins of short-Gene-Operons with Acr homologs |

****

<div id='examples' />

## **IV. <span style='color:RebeccaPurple'>Examples</span>**

```sh
python3 AcaFind_runner.py --FNA_file sample_organism/GCF_000381965.1_ASM38196v1_genomic.fna --GFF_file sample_organism/GCF_000381965.1_ASM38196v1_genomic.gff --FAA_file sample_organism/GCF_000381965.1_ASM38196v1_protein.faa -o [output_dir] 
```
or you can only use **.fna** file as input.

```sh
python3 AcaFind_runner.py --FNA_file sample_organism/GCF_000381965.1_ASM38196v1_genomic.fna -o [output_dir] 
```

You will see the output result in output_dir/. If you dont specifiy an output_dir, result will be in AcaFinder_Output/


****

<div id="workflow" />

## **V. <span style='color:RebeccaPurple'>Workflow of AcaFinder</span>**

<img src="http://bcb.unl.edu/AcaFinder/img/Figure2_Pipeline.png">


With provided input, AcaFinder proceed to 2 Aca screening routes:

i) HMM approach Aca-like protein find (purple box)

ii) GBA (guilt-by-association) approach Aca protein/operon find (orange box), generating predicted Aca operons/Aca proteins. 

Complete CRISPR-Cas along with STSS (blue box) and prophage regions (green box) will also be searched from input genomic sequences. 

**All generated information will be associated together, and provided to the users as tables.** 

****
