#!/bin/bash  
#-*- coding:utf-8 -*-  
############################  
#File Name: Conservation.sh
#Author: ShiFY
#Mail: shify@mail.cbi.pku.edu.cn  
#Created Time: 2018-05-30 13:29:28
############################  

bed_file=$1

data_type="phastCons100way"
output=${2}${data_type}".txt"
hgWiggle -db=hg19 -bedFile=${bed_file} phastCons100way  >> ${output}

data_type="phyloP46wayAll"
output=$2${data_type}".txt"
hgWiggle -db=hg19 -bedFile=${bed_file} phyloP46wayAll  >> ${output}


data_type="phyloP46wayPrimates"
output=$2${data_type}".txt"
hgWiggle -db=hg19 -bedFile=${bed_file} phyloP46wayPrimates  >> ${output}



data_type="phastCons46wayPlacental"
output=$2${data_type}".txt"
hgWiggle -db=hg19 -bedFile=${bed_file} phastCons46wayPlacental  >> ${output}

data_type="phastCons46way"
output=$2${data_type}".txt"
hgWiggle -db=hg19 -bedFile=${bed_file} phastCons46way  >> ${output}

data_type="phastCons46wayPrimates"
output=$2${data_type}".txt"
hgWiggle -db=hg19 -bedFile=${bed_file} phastCons46wayPrimates  >> ${output}

data_type="phyloP46wayPlacental"
output=$2${data_type}".txt"
hgWiggle -db=hg19 -bedFile=${bed_file} phyloP46wayPlacental  >> ${output}

data_type="phyloP100way"
output=$2${data_type}".txt"
hgWiggle -db=hg19 -bedFile=${bed_file} phyloP100wayAll  >> ${output}

