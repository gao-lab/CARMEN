#!/usr/bin/python  
#-*- coding:utf-8 -*-  
############################  
#File Name: 01_12features_annotation.py
#Author: ShiFY
#Mail: shify@mail.cbi.pku.edu.cn  
#Created Time: 2017-05-16 14:09:12
############################  
"""
This script is used to calculate the 12 PC features from sequence, fasta file is from 04_coor2fasta.R 
"""

import numpy as np
import pandas as pd
import math
import multiprocessing
import sys

Feature_12_score = pd.read_table("12features_PC.txt") #Read the annotation matrix score from supplementary file.

def get_output(sequence):
    seq_ori = sequence[1][97:102]
    seq_var = sequence[2][97:102]
    
    if "N" in seq_ori or "N" in seq_var:
        delta_T_12 = np.zeros(12)

    else:
    	def generate_seed(seq):
        	seq_list = map(lambda x : seq[x:x+2],range(len(seq)-1))
        	seq_matrix = Feature_12_score[seq_list]
        	T_score = map(lambda x : sum(seq_matrix.iloc[x]),range(12))
        	return T_score

        T_ori = generate_seed(seq_ori)
        T_var = generate_seed(seq_var)
    	
        def calculate_ori_var_score(TA,TB):
    	    if TA*TB > 0 or TA*TB == 0:
    		    delta_T = math.log(abs(TA)+0.00001) - math.log(abs(TB)+0.00001)
    	    else:
    		    delta_T = math.log(abs(TA)+0.00001) + math.log(abs(TB)+0.00001)
    	    return delta_T
    	
        delta_T_12 = map(lambda x,y : calculate_ori_var_score(x,y),T_ori,T_var)
    return delta_T_12


if __name__=='__main__':
    in_file = sys.argv[1]
    inputs = np.loadtxt(in_file,dtype=str)
    if inputs.shape == (3,):
	inputs = inputs.reshape(1,3)
    pool = multiprocessing.Pool(processes=15)
    pool_outputs = pool.map(get_output,inputs)
    pool.close()
    pool_outputs_matrix = np.array(pool_outputs) 
    np.savetxt(in_file.replace('.txt','_PC_features.txt'),pool_outputs_matrix,delimiter = '\t',fmt='%.6e',header='Enthalpy\tEntropy\tInclination\tMinor_Groove_Size\tProbability_contacting_nucleosome_core\tRise\tSlide_rise\tSlide_slide\tTilt_stiffness\tTilt_shift\tTilt_tilt\tTwist_stiffness')

