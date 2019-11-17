#!/bin/bash  
#-*- coding:utf-8 -*-  
############################  
#File Name: Conservation.sh
#Author: ShiFY
#Mail: shify@mail.cbi.pku.edu.cn  
#Created Time: 2018-05-30 13:29:28
############################  

bed_file=$1

for iterm in {phastCons100way,phyloP46wayAll,phyloP46wayPrimates,phastCons46wayPlacental,phastCons46way,phastCons46wayPrimates,phyloP46wayPlacental,phyloP100wayAll}
do
{
data_type=${iterm}
output=${2}${data_type}".txt"
../Tools/hgWiggle -bedFile=${bed_file} ${iterm}  >> ${output}
}&
done
wait
