# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 11:58:00 2022

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
import gzip
import shutil
import xlsxwriter
import numpy as np
from statistics import mean
from scipy.integrate import simps
from scipy.interpolate import UnivariateSpline


def FWHM_curve(x,y):
    '''This function returns the full width half maximum of a function y of x having a peak.'''
    spline = UnivariateSpline(x, y-np.max(y)/2, s=0)
    a = spline.roots() # find the roots
    r1, r2 = a[0], a[1]
    return abs(r1-r2)

def AIF_shape_metrics(Time, AIF):
    '''This function calculates the following curve shape parameters of an AIF:
        maximum peak: MP
        full width at half maximum: FWHM
        time to peak: TTP
        area under the curve: AUC
        M=MP/(TTP x FWHM): M
        -----------------------
        Time, AIF : ndarray
        '''   
    MP = max(AIF) 
    FWHM = FWHM_curve(Time,AIF)
    TTP = Time[AIF.argmax()] 
    AUC = simps(AIF, Time)
    M = MP/(TTP*FWHM)
    
    return {'MP':round(MP,2), 'FWHM':round(FWHM,2), 'TTP':round(TTP,2), 'AUC':round(AUC,2), 'M':round(M,2)}

#%% Set Working directory
        
# Path of the main folder 
DataSupraDir = 'Path to users supradirectory'

reters_path = [ f.path for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of the paths to the raters directories
raters_name = [ f.name for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of raters names
raters_name.remove('Analysis results')

Results_AIF = {raters_name[0]: pd.DataFrame(columns=['Subject_ID', 'MP', 'FWHM', 'TTP', 'AUC', 'M']), raters_name[1]: pd.DataFrame(columns=['Subject_ID','MP', 'FWHM', 'TTP', 'AUC', 'M']), raters_name[2]: pd.DataFrame(columns=['Subject_ID','MP', 'FWHM', 'TTP', 'AUC', 'M'])}
IntraWriter = pd.ExcelWriter(DataSupraDir +'Analysis results/AIF analysis/AIF_curve_parameters.xlsx', engine='xlsxwriter')

for rater in raters_name:

    n_raters = len(raters_name) #Total number of raters

    # Extract number of patient folders as a list 
    subjs_name = [ f.name for f in os.scandir(DataSupraDir+rater) if f.is_dir() ] #Create a list of raters names
    n_subjs = len(subjs_name) #Total number of subjects
    
    #%%Create a for loop to perform image analysis on each rater sequentially
    
    for subject in subjs_name:
        
        print('Analysing patient '+subject)
        
        rater_dir = DataSupraDir+rater          
                
        # Set path to patient subfolders
        patientFolderPath = os.path.join(rater_dir,subject)
            
        # Set path to AIF.xlsx file
        AIFpath = os.path.join(patientFolderPath,'AIF.xlsx')
            
        # Read AIF data in a dataframe
        AIF_df = pd.read_excel(AIFpath, sheet_name='Sheet1')
        Time = np.array(AIF_df.loc[:,'Time (s)'])
        AIF = np.array(AIF_df.loc[:, ' Average AIF'])
                
        # Calculate curve shape metrics for AIF    
        AIF_param_dic = AIF_shape_metrics(Time, AIF)
        Results_AIF[rater] = Results_AIF[rater].append({'Subject_ID': subject,'MP':AIF_param_dic.get('MP'), 'FWHM':AIF_param_dic.get('FWHM'), 'TTP':AIF_param_dic.get('TTP'), 'AUC':AIF_param_dic.get('AUC'), 'M':AIF_param_dic.get('M')}, ignore_index=True)
                    
#%%Save all dataframes to excel files here
print('Save all results to excel files')
    
for name, df in Results_AIF.items():
    df.to_excel(IntraWriter, sheet_name=name, index=False)
IntraWriter.save()
