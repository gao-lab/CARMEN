#!/bin/bash  
#-*- coding:utf-8 -*-  
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
