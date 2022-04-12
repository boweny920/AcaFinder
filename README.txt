1. VIBRANT
Anaconda (currently running v1.2.0)
	1. Install dependencies. See Requirements section above.
	2. Install directly to $PATH using bioconda.
conda install -c bioconda vibrant==1.2.0
	3. Download and setup databases. This will take some time due to file sizes, but it only needs to be run once. This step requires ~20GB of temporary storage space and ~11GB of final storage space. To do this, run download-db.sh which should be in your system's $PATH.
download-db.sh

2. Cctyper

version xgboost>=1.4

python -m pip install cctyper
Or 
Use conda to install cctyper: conda install -c russel88 cctyper

3. Pfam database
a. Go to the all_pFam_hmm/ folder
b. wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
c. gunzip Pfam-A.hmm.gz
d. hmmpress Pfam-A.hmm
 

