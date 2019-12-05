#!/usr/bin/python  
#-*- coding:utf-8 -*-  
############################  

from subprocess import *
from tempfile import mkdtemp
import sys
import os 
import time
start = time.clock()

infilename = sys.argv[1]
outdir = sys.argv[2]
vcf_version = sys.argv[3]
job_id = infilename.split("/")[-1].split(".")[0]
python_path = os.path.abspath('Tools/') + "/anaconda2/bin/python " 
R_path = os.path.abspath('Tools/') + "/anaconda2/bin/Rscript "

cpoutdir=True
check_call(['mkdir','-p',outdir])
precent_dir=os.getcwd()

if infilename.endswith('vcf'):
    try:
        tempdir = mkdtemp()
        
        check_call(["sed 's/chr//g' " + infilename + " | sed 's/^/chr/g' | sort -k1,1 -k2n > " + tempdir + '/infile.tmp.vcf'],shell=True)
        check_call([R_path + "01_check_format.R " + tempdir+"/infile.tmp.vcf "+ tempdir + ' ' + vcf_version],shell=True)

        print "Successfully copied input vcf file to working directory " + tempdir
#===========================GRCh38 sequence orientated annotation =============== 
        ''' 
        Get Fasta file
        '''
        check_call([R_path + "02_coor2fasta.R " + tempdir + "/infile.hg38.vcf 300"],shell=True)
    
        check_call([python_path + "03_fasta2wtandmu.py " + tempdir + "/infile.hg38.vcf.wt200.fasta 300"],shell=True)

        check_call(["paste -d '\t' " + tempdir + "/infile.hg38.vcf " + tempdir + "/infile.hg38.vcf.wt200.fasta.ref.fasta " + tempdir + "/infile.hg38.vcf.mut200.fasta.ref.fasta > " + tempdir + "/infile.hg38.vcf.merged.tmp"],shell=True)
    
        check_call(["cat merge_file_title.txt " + tempdir + "/infile.hg38.vcf.merged.tmp > " + tempdir + "/infile.hg38.vcf.merged.txt"],shell=True)
        '''
        Get annotations with sequence orientated features
        '''
        os.chdir("./DeepOcean_web/")
        check_call(["bash run_predict.sh " + tempdir + "/infile.hg38.vcf.merged.txt " + job_id + " " + python_path],shell=True)
        check_call(["mv ./output/" + job_id + "/infile.log.foldchange " + tempdir + "/infile.log.foldchange"],shell=True)
        call(['rm',"./output/" + job_id,'-r'])
        
#=============================GRCh37 physicalchemical and evolutionary features===============
        os.chdir(precent_dir)
        check_call([R_path + "04_coor2fasta.R " + tempdir + "/infile.hg19.vcf"],shell=True)
        check_call([python_path + "05_fasta2input.py " + tempdir + "/infile.hg19.vcf.wt200.fasta"],shell=True)
        check_call(["paste -d '\t' " + tempdir + "/infile.hg19.vcf.wt200.fasta.ref.vcf " + tempdir + "/infile.hg19.vcf.wt200.fasta.ref.fasta " + tempdir + "/infile.hg19.vcf.mut200.fasta.ref.fasta | cut -f 1,6,7 > " + tempdir + "/infile.hg19.wt200.txt"],shell=True)

        print "Successfully got fasta to get physicalchemical and evolutionary features"


        '''
        physicalchemical and OH properties
        '''

        check_call([python_path + "06_12features_annotation.py " + tempdir + "/infile.hg19.wt200.txt"], shell=True)
        check_call([python_path + "07_OHfeature_annotation.py " + tempdir + "/infile.hg19.wt200.txt"], shell=True)

        
        '''
        conservation score
        '''
        os.chdir("./Conservation/")
        check_call(["bash Conservation_parallele.sh " + tempdir + "/infile.hg19.bed " + tempdir + "/"], shell=True)
        
        os.chdir(precent_dir)

        for coniterm in ["phastCons100way","phastCons46wayPlacental","phastCons46wayPrimates","phastCons46way","phyloP100wayAll","phyloP46wayAll","phyloP46wayPlacental","phyloP46wayPrimates"]:
            check_call([python_path+"08_conservation_score_merge.py " + tempdir + "/" + coniterm + ".txt " + tempdir + "/infile.hg19.vcf"], shell=True)

        os.chdir(tempdir)
        check_call(["paste -d '\t' phastCons100way.format.txt phastCons46wayPlacental.format.txt phastCons46wayPrimates.format.txt phastCons46way.format.txt phyloP100wayAll.format.txt phyloP46wayAll.format.txt phyloP46wayPlacental.format.txt phyloP46wayPrimates.format.txt > infile.conservation.tmp"],shell=True)
        check_call(["""awk '{for (i=2;i<NF;i++) { if (i%2 == 0) printf "%s\\t",$i }; if (NF%2 == 0) {print $NF} else printf "\\n"}' infile.conservation.tmp > infile.conservation.txt"""],shell=True)
        check_call(["cat " + precent_dir + "/conservation_title_2424.txt " + tempdir + "/infile.conservation.txt > " + tempdir + "/infile.conservation.title.txt"],shell=True)

#==================================Merged annotations and make prediction======================
        check_call(["less infile.hg19.vcf | cut -f 2 | sed '1i hg19pos' > hg19_pos.txt"],shell=True)

        check_call(["paste -d '\t' hg19_pos.txt infile.log.foldchange infile.hg19.wt200_PC_features.txt infile.hg19.wt200_OH_features.txt infile.conservation.title.txt | sed 's/# //g'  > infile.annotation.txt"],shell=True)
        os.chdir(precent_dir)
        check_call([python_path + "09_merge_all_features.py " + tempdir + "/infile.annotation.txt " + tempdir + "/infile.annotation.689.hdf5 " + tempdir + "/infile.annotation.1190.hdf5"], shell=True)
        check_call([python_path + "10_pinpoint_emVar.py " + tempdir + "/infile.annotation.689.hdf5 " + tempdir + "/infile.annotation.1190.hdf5"],shell=True)
        check_call(["cp", tempdir + "/infile.annotation.txt" ,outdir ])
        check_call(["cp", tempdir + "/infile.annotation.CARMEN.predict.csv" ,outdir ])
    except:
        raise Exception("VCF format error.")

    finally:
        print "Finished creating output file. Now clean up..."
        call(['rm',tempdir,'-r'])
        print "Everything done!" 
