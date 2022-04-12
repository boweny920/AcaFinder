import subprocess
import os

class annotation_prodigal:
    def __init__(self,fna_file,run_mode,prodigal_outDir):
        self.fna=fna_file
        self.mode=run_mode
        self.out_dir=prodigal_outDir
    def run_prodigal(self):
        os.makedirs(self.out_dir)
        subprocess.Popen(["prodigal","-a",os.path.join(self.out_dir,os.path.basename(str(self.fna))+"_prodigal_out.faa"),"-p",self.mode,"-f","gff",
                          "-i",self.fna,"-d",os.path.join(self.out_dir,os.path.basename(str(self.fna))+"_prodigal_out.fna"),
                          "-o",os.path.join(self.out_dir,os.path.basename(str(self.fna))+"_prodigal_out.gff")]).wait()
        print("Genome annotated by Prodigal, files stored in %s"%self.out_dir)
        #returned gff, faa, fna
        return os.path.join(self.out_dir,os.path.basename(str(self.fna))+"_prodigal_out.gff"), os.path.join(self.out_dir,os.path.basename(str(self.fna))+"_prodigal_out.faa")


