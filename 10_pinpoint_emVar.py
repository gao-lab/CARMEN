import os
import sys
import h5py
import numpy as np
import numpy
import pandas as pd
import time
import itertools
import gzip
from itertools import cycle
from sklearn.externals import joblib
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier


def makdir(path):
        isExists=os.path.exists(path)
        if not isExists:
                os.makedirs(path)
                return(True)
        else:
                return(False)

# load data #
def load_data(data_path):
        data = h5py.File(data_path,'r')
        training_data = data['data_tensor'].value
        information = data['information_tensor'].value
        label = data['label_tensor'].value
        return([information,training_data,label])



#==================================Model prediction=====================
data_path_689 = sys.argv[1]
data_689 = h5py.File(data_path_689)
info_689 = pd.DataFrame(data_689['information_tensor'].value)
anno_in_689 = data_689['data_tensor'].value

clf = joblib.load('./emVar_model_parameters/AdaBoost_800_estimators2.pkl')
predict_Y_689 = pd.DataFrame(clf.predict_proba(anno_in_689)[:,1])

data_path_1190 = sys.argv[2]
data_1190 = h5py.File(data_path_1190)
info_1190 = pd.DataFrame(data_1190['information_tensor'].value)
anno_in_1190 = data_1190['data_tensor'].value

with gzip.GzipFile('./emVar_model_parameters/RF_1200_1KG-HGMDDM-sampleweight_estimators4.pkl.gz', 'rb') as f:
        rlf = joblib.load(f)

predict_Y_1190 = pd.DataFrame(rlf.predict_proba(anno_in_1190)[:,1])

#===================================Bayesian correct proba==============

CARMEN_score = predict_Y_689*predict_Y_1190/(predict_Y_689*predict_Y_1190+(1-predict_Y_689)*(1-predict_Y_1190))

#===================================Output file=========================

out = pd.concat([info_689,predict_Y_689,predict_Y_1190,CARMEN_score],axis=1)

out.columns = ['GRCh37','Chr','GRCh38','rs_id','Ref','Alt','Comment','CARMEN-E','CARMEN-F','CARMEN']
out = out[['Chr','GRCh37','GRCh38','rs_id','Ref','Alt','Comment','CARMEN-E','CARMEN-F','CARMEN']]

out.to_csv(data_path_689.replace(".689.hdf5",".CARMEN.predict.csv"),sep=',')
