#!/usr/bin/python  
#-*- coding:utf-8 -*-  
############################  
#File Name: 02_OHfeature_annotation.py
#Author: ShiFY
#Mail: shify@mail.cbi.pku.edu.cn  
#Created Time: 2017-05-16 16:51:37
############################  
"""
This script is used to caluclate the OH feature, fasta file is also from 04_coor2fasta.R 
"""

import numpy as np
import pandas as pd
import multiprocessing
import sys

Feature_OH_score = pd.read_table("OH_feature.txt",index_col=0,header=None) #Read the annotation matrix score from supplementary file OH.

def get_output(sequence):
    
    if "N" in sequence[1][93:106] or "N" in sequence[2][93:106]:
        dalta_oh = np.zeros(1)
    else:
        reverse_dict={'A': 'T','C': 'G','G': 'C','T': 'A'} #Construct a dict to make the complement nucletide
        seq_ori = sequence[1][93:106] # original sequence
        seq_var = sequence[2][93:106] # variant sequence
        if ((len(seq_ori) == 13) & (len(seq_var) == 13)):
            seq_ori_rev = "".join(map(lambda x : reverse_dict[x], seq_ori))[::-1] # make the reverse of complement nucletide
            seq_var_rev = "".join(map(lambda x : reverse_dict[x], seq_var))[::-1]
        
            def get_OH_vetor_all(seq_o,seq_r):
                '''This function is used to calculate the matrix of OH score, is 10x13 '''
                def get_OH_vector(seq):
                    seqlen=len(seq_o)
                    valuearray=np.zeros([seqlen-3,seqlen])
                    for j in range(0,seqlen-3):#plus
                            tetramer=seq[j:j+4]
                            valuearray[j][j:j+4]=Feature_OH_score.ix[tetramer][1:5]
                    return valuearray
                
                ohvaluearray = get_OH_vector(seq_o) # get the original sequence matrix
                r_ohvaluearray = get_OH_vector(seq_r) # get the reverse sequence matrix
                
                [m,n]=np.shape(ohvaluearray)
                fohdata =[]
                for i in range(0,n):
                    fohdata.append(sum(ohvaluearray[:,i])/len(np.nonzero(ohvaluearray[:,i])[0])) # average each column in sequence matrix
                ohdata=fohdata
                rohdata =[]
                for i in range(0,n):
                    rohdata.append(sum(r_ohvaluearray[:,i])/len(np.nonzero(r_ohvaluearray[:,i])[0])) # average each column
                for i in range(3,n):
                    ohdata[i]=(fohdata[i]+rohdata[n+2-i])/2 # only the upstream 3bp and downstream 6bp is effected by the variants
                return ohdata
            ohdata1 = np.array(get_OH_vetor_all(seq_ori,seq_ori_rev))
            ohdata2 = np.array(get_OH_vetor_all(seq_var,seq_var_rev))
            dalta_oh = np.sqrt(np.dot((ohdata1-ohdata2),(ohdata1-ohdata2).reshape(ohdata2.shape[0],1)))
        else:
            dalta_oh = np.zeros(1)
    return dalta_oh
if __name__=='__main__':
    in_file = sys.argv[1]
    inputs = np.loadtxt(in_file,dtype=str)
    if inputs.shape == (3,):
        inputs = inputs.reshape(1,3)
    pool = multiprocessing.Pool(processes=15)
    pool_outputs = pool.map(get_output,inputs)
    pool.close()
    pool_outputs_matrix = np.array(pool_outputs)
    np.savetxt(in_file.replace('.txt','_OH_features.txt'),pool_outputs_matrix,delimiter='\n',header='OH')

