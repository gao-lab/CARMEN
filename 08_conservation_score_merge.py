#!/usr/bin/python  
#-*- coding:utf-8 -*-  
############################  
#File Name: 08_conservation_score_merge.py
#Author: ShiFY
#Mail: shify@mail.cbi.pku.edu.cn  
#Created Time: 2017-05-22 14:07:20
############################  
import sys
import re
import numpy as np
import pandas as pd

input_file = sys.argv[1]


# if has # continue
# if has v i = i + 1

file = open(input_file, 'r')
i = 0
final_list = []
for line in file:
    line = line.strip()
    anno=re.match(r'^\#',line)
    chr_num=re.match(r'^variableStep',line)
    if anno:
        continue
    if chr_num:
        chro = re.findall(r'chr\d{1,}|chr[X]|chr[Y]',line)[0]
        continue
    pos = line.split('\t')[0]
    score = line.split('\t')[1]
    #bed file has 0 base and vcf file has 1 base ,so the vcf-1 can get right conservation score.
    iterm = [chro,str(int(pos)+1),score]
    final_list.append(iterm)
    #print 'chr'+str(i)+'\t'+str(pos)+'\t'+str(score)
file.close()

final_form = np.array(final_list)

input_id = map(lambda chro,pos : '_'.join((chro,pos)),final_form[:,0],final_form[:,1])
input_df_1 = pd.DataFrame(data=final_form[:,2],columns=['conservation_score'])
input_df_2 = pd.concat((pd.DataFrame(input_id),input_df_1),axis=1)
input_df = input_df_2.drop_duplicates(subset = 0)
input_df = input_df.set_index(input_df[0])

direct_file = sys.argv[2]
direct_data = np.loadtxt(direct_file,dtype=str,delimiter='\t')
if direct_data.shape == (6,):
    direct_data = direct_data.reshape(1,6)
direct_id = map(lambda chro,pos : '_'.join((chro,pos)),direct_data[:,0],direct_data[:,1])

out_file = input_file.replace('txt','format.txt')
direct_out = input_df.loc[direct_id]
direct_out.to_csv(out_file,sep='\t',columns=None,header=False,index=False) 
