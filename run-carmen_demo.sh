#!/bin/bash  
#-*- coding:utf-8 -*-  
############################  
#Reference version "GRCh37" or "GRCh38".
#Annotation version "all" or "basic".
#Option "all" : make annotations with 2424 features.
#Option "basic" : make annotations with 1393 features. 

start=$(date +%s)

./Tools/anaconda2/bin/ipython run-carmen.py ./example/example.vcf example-result GRCh37 basic >> ./LOG/example-$(date "+%y%m%d")_log.txt 2>&1

end=$(date +%s)
time_cost=$(($end - $start))
echo "Finsied prediction! The prediction cost "$[$time_cost/60]" mins."
