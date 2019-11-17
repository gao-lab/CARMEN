#!/bin/bash  
#-*- coding:utf-8 -*-  
############################  

start=$(date +%s)

./Tools/anaconda2/bin/ipython run-carmen.py ./example/example.vcf example-result GRCh37 >> ./LOG/example-$(date "+%y%m%d")_log.txt 2>&1

end=$(date +%s)
time_cost=$(($end - $start))
echo "Finsied prediction! The prediction cost "$[$time_cost/60]" mins."
