import numpy as np
import pandas as pd
import h5py
import os
import sys

#=================================Training file with label===========================

datafile = sys.argv[1]

f = open(datafile,'r')
final_list = []

for line in f:
    eline = line.strip()
    line_list = eline.split('\t')
    if len(line_list) == 2431:
        final_list.append(line_list)
    else:
        continue
f.close()

final_form = pd.DataFrame(final_list).iloc[1:,:]


information = np.array(final_form.iloc[:,0:7],dtype=str)
label = np.array(final_form.iloc[:,6],dtype=str)
training = np.array(final_form.iloc[:,7:2431],dtype='float32')

feature_list = np.loadtxt('feature_list_689.txt')

feature_list = feature_list.astype(int)
feature_list_1 = feature_list -1
training_new = training[:,feature_list_1]

save_file_path = sys.argv[2]
training_data_out = h5py.File(save_file_path,'w')
training_data_out.create_dataset('data_tensor',data=training_new)
training_data_out.create_dataset('information_tensor',data=information)
training_data_out.create_dataset('label_tensor',data=label)
training_data_out.close()

feature_list = np.loadtxt('feature_list_1190.txt')

feature_list = feature_list.astype(int)
feature_list_1 = feature_list -1
training_new = training[:,feature_list_1]

save_file_path = sys.argv[3]
training_data_out = h5py.File(save_file_path,'w')
training_data_out.create_dataset('data_tensor',data=training_new)
training_data_out.create_dataset('information_tensor',data=information)
training_data_out.create_dataset('label_tensor',data=label)
training_data_out.close()
