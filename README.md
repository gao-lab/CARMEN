# Installation Guide

## Stage 1: Download and uncompress the core package

[user@host ~]$ gzip -dc CARMEN.tar.gz | tar x

[user@host ~]$ cd CARMEN

## Stage 2: Download code and data for annotation module (a LARGE file, ~6GB)

[user@host ~]$ wget http://carmen.gao-lab.org/download/DeepOcean_web.tar.gz

[user@host ~]$ tar -zxvf DeepOcean_web.tar.gz

## Stage 3: Setup 3rd-party toolchain

### 3.1 Anaconda

We recommend Anaconda2-5.0-1.

[user@host ~]$ cd ./Tools/

[user@host ~]$ wget https://repo.continuum.io/archive/Anaconda2-5.0.1-Linux-x86_64.sh 

[user@host ~]$ bash Anaconda2-5.0.1-Linux-x86_64.sh

"""
Do you accept the license terms?
[1]yes

Anaconda2 will now be installed into this location:
[2]${abs_path}+/anaconda2

Do you wish the installer to prepend the Anaconda2 install location?
[3]no
"""

[user@host ~]$ ./anaconda2/bin/conda install mkl-rt

[user@host ~]$ ./anaconda2/bin/pip install tensorflow==1.4.0

[user@host ~]$ ./anaconda2/bin/pip install scikit-learn==0.19.1

[user@host ~]$ ./anaconda2/bin/pip install biopython

[user@host ~]$ ./anaconda2/bin/conda install -c r r=3.2.2 #INSTALL R-3.2.2

[user@host ~]$ ./anaconda2/bin/R

"Install within R"

source("http://bioconductor.org/biocLite.R")

biocLite("BSgenome.Hsapiens.UCSC.hg19")

biocLite("BSgenome.Hsapiens.UCSC.hg38")

"If you did not use the python in anaconda that we recommended, please remember to change the path of python in run-carmen.py"
    
### 3.2 Conservation via UCSC genome browse.

[user@host ~]$ wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgWiggle

[user@host ~]$ chmod +x hgWiggle

[user@host ~]$ cd ../Conservation/

[user@host ~]$ bash Conservation_download.sh


## Stage 4: Run Demo

[user@host ~]$ cd ..

[user@host ~]$ bash run-carmen_demo.sh

Congratuate if you can see output like

```
Finsied prediction! The prediction cost 1 mins.
```

Or, please just feel free to let us know by email carmen@mail.cbi.pku.edu.cn :)

# Notes and Tips

While the core package needs ~50MB, more than 50GB space will be required eventually.

With all pre-required infrastructure set up correctly, software installation could be done in minutes.

We also recommend download conservation data and installation packages at the same time to save the installation time.

# Contact:

Website: http://carmen.gao-lab.org

Email: carmen@mail.cbi.pku.edu.cn
