# -*- coding: utf-8 -*-
"""
Created on Sun Nov 1 12:21:00 2020

@author: Caterina Brighi
"""

#%% Import functions 

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import datetime
import os
import glob
import shutil
import xlsxwriter
import time
from scipy.stats.stats import pearsonr

def createList(r1, r2): 
    return [item for item in range(r1, r2+1)] 

#%% Set Working directory
        
# Path of the main folder 
DataSupraDir = 'Path to users supradirectory'+'/Analysis results/Similarity analysis/Tumor/'
# DataSupraDir = 'Path to users supradirectory'+'/Analysis results/Similarity analysis/Brain/'

#Read three comparisons spreadsheets into dataframes
EC_AK_results_dic = pd.read_excel(DataSupraDir+'AthenaKam-EdwardChen_Tumor_Similarity_Results.xlsx', sheet_name=None, index_col='Subject_ID')
EC_SD_results_dic = pd.read_excel(DataSupraDir+'EdwardChen-SameerDave_Tumor_Similarity_Results.xlsx', sheet_name=None, index_col='Subject_ID')
SD_AK_results_dic = pd.read_excel(DataSupraDir+'AthenaKam-SameerDave_Tumor_Similarity_Results.xlsx', sheet_name=None, index_col='Subject_ID')
# EC_AK_results_dic = pd.read_excel(DataSupraDir+'AthenaKam-EdwardChen_Brain_Similarity_Results.xlsx', sheet_name=None, index_col='Subject_ID')
# EC_SD_results_dic = pd.read_excel(DataSupraDir+'SameerDave-EdwardChen_Brain_Similarity_Results.xlsx', sheet_name=None, index_col='Subject_ID')
# SD_AK_results_dic = pd.read_excel(DataSupraDir+'AthenaKam-SameerDave_Brain_Similarity_Results.xlsx', sheet_name=None, index_col='Subject_ID')

Imgs_list = list(EC_AK_results_dic.keys())

Subj_list = EC_AK_results_dic.get(Imgs_list[0]).index.tolist()

metric_list = ['MSE', 'RMSE', 'NRMSE', 'PSNR', 'SSIM']

storeResults = {metric: pd.DataFrame(columns=['Comparison type', 'Subject_ID', 'rBV_noLC', 'rBV_LC', 'rBF_LC', 'MTT_LC', 'n_rBV_noLC', 'n_rBV_LC', 'n_rBF_LC', 'n_MTT_LC']) for metric in metric_list}

comp_1, comp_2, comp_3 ='Expert-Non expert','Expert-Trainee','Trainee-Non expert'

for metric in metric_list:
    
    for subj in Subj_list:

        storeResults[metric] = storeResults[metric].append({'Comparison type': comp_1,
                                                            'Subject_ID': subj,
                                                            'rBV_noLC':EC_AK_results_dic.get('rBV_noLC').loc[subj, metric], 
                                                            'n_rBV_noLC':EC_AK_results_dic.get('n_rBV_noLC').loc[subj, metric], 
                                                            'rBV_LC':EC_AK_results_dic.get('rBV_LC').loc[subj, metric], 
                                                            'n_rBV_LC':EC_AK_results_dic.get('n_rBV_LC').loc[subj, metric], 
                                                            'rBF_LC':EC_AK_results_dic.get('rBF_LC').loc[subj, metric], 
                                                            'n_rBF_LC':EC_AK_results_dic.get('n_rBF_LC').loc[subj, metric], 
                                                            'MTT_LC':EC_AK_results_dic.get('MTT_LC').loc[subj, metric], 
                                                            'n_MTT_LC':EC_AK_results_dic.get('n_MTT_LC').loc[subj, metric]
                                                         }, ignore_index=True)
    
    for subj in Subj_list:

        storeResults[metric] = storeResults[metric].append({'Comparison type': comp_2,
                                                            'Subject_ID': subj,
                                                            'rBV_noLC':EC_SD_results_dic.get('rBV_noLC').loc[subj, metric], 
                                                            'n_rBV_noLC':EC_SD_results_dic.get('n_rBV_noLC').loc[subj, metric], 
                                                            'rBV_LC':EC_SD_results_dic.get('rBV_LC').loc[subj, metric], 
                                                            'n_rBV_LC':EC_SD_results_dic.get('n_rBV_LC').loc[subj, metric], 
                                                            'rBF_LC':EC_SD_results_dic.get('rBF_LC').loc[subj, metric], 
                                                            'n_rBF_LC':EC_SD_results_dic.get('n_rBF_LC').loc[subj, metric], 
                                                            'MTT_LC':EC_SD_results_dic.get('MTT_LC').loc[subj, metric], 
                                                            'n_MTT_LC':EC_SD_results_dic.get('n_MTT_LC').loc[subj, metric]
                                                         }, ignore_index=True)
        
    for subj in Subj_list:

        storeResults[metric] = storeResults[metric].append({'Comparison type': comp_3,
                                                            'Subject_ID': subj,
                                                            'rBV_noLC':SD_AK_results_dic.get('rBV_noLC').loc[subj, metric], 
                                                            'n_rBV_noLC':SD_AK_results_dic.get('n_rBV_noLC').loc[subj, metric], 
                                                            'rBV_LC':SD_AK_results_dic.get('rBV_LC').loc[subj, metric], 
                                                            'n_rBV_LC':SD_AK_results_dic.get('n_rBV_LC').loc[subj, metric], 
                                                            'rBF_LC':SD_AK_results_dic.get('rBF_LC').loc[subj, metric], 
                                                            'n_rBF_LC':SD_AK_results_dic.get('n_rBF_LC').loc[subj, metric], 
                                                            'MTT_LC':SD_AK_results_dic.get('MTT_LC').loc[subj, metric], 
                                                            'n_MTT_LC':SD_AK_results_dic.get('n_MTT_LC').loc[subj, metric]
                                                         }, ignore_index=True)
        

writer = pd.ExcelWriter(DataSupraDir +'Simil_anal_reordered.xlsx', engine='xlsxwriter')

for name, df in storeResults.items():

    df.to_excel(writer, sheet_name=name, index=False)
writer.save()

            
            
            
            
