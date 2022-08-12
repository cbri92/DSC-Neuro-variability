# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 14:37:53 2022

@author: Caterina Brighi

This script evaluates inter-rater repeatability in a ROI between images of the same subject generated from multiple raters.
Repeatability is assessed via ICC(2,K) based on a mean-rating (k= 3) from multiple raters, absolute-agreement, 2-way random effects model. 
"""
#%% Import functions 

import matplotlib.pyplot as plt
import matplotlib as mpl
import SimpleITK as sitk
import numpy as np
import pandas as pd
import datetime
import os
import glob
import gzip
import shutil
import xlsxwriter
from scipy.stats.stats import pearsonr
from multiprocessing.pool import ThreadPool
from functools import partial
from radiomics import featureextractor
import six, numpy as np
from statistics import pvariance
from statistics import mean
from scipy.stats import normaltest
import seaborn as sns
import pingouin as pg
from alive_progress import alive_bar
from time import sleep
    
    
#%% Set Working directory
        
# Path of the main folder 
DataSupraDir = 'Path to users supradirectory'


Results_ICC = {'rBV_noLC': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'rBV_LC': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'rBF_LC': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'MTT_LC': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'n_rBV_noLC': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'n_rBV_LC': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'n_rBF_LC': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'n_MTT_LC': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%'])}
IntraWriter = pd.ExcelWriter(DataSupraDir +'Brain_ICCResults.xlsx', engine='xlsxwriter')


reters_path = [ f.path for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of the paths to the raters directories
raters_name = [ f.name for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of raters names

n_raters = len(raters_name) #Total number of raters


# Extract number of patient folders as a list 
subjs_name = [ f.name for f in os.scandir(DataSupraDir+raters_name[0]) if f.is_dir() ] #Create a list of raters names
n_subjs = len(subjs_name) #Total number of subjects

#%%Create a for loop to perform image analysis on each rater sequentially

for subject in subjs_name:
    
    print('Analysisng patient '+subject)
    
    #Create empty dataframes to populate as going through the loop
    
    rBV_noLC_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'rBV_noLC'], dtype=float)
    rBV_LC_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'rBV_LC'], dtype=float)
    rBF_LC_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'rBF_LC'], dtype=float)
    MTT_LC_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'MTT_LC'], dtype=float)    
    
    n_rBV_noLC_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'n_rBV_noLC'], dtype=float)
    n_rBV_LC_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'n_rBV_LC'], dtype=float)
    n_rBF_LC_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'n_rBF_LC'], dtype=float)
    n_MTT_LC_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'n_MTT_LC'], dtype=float)   
    
    # Loop through each rater 
    for idx, rater in enumerate(raters_name):
        print(idx, rater)
        rater_dir = DataSupraDir+rater

    
        print('Calculating ICC for rater:',rater)
        with alive_bar(6) as bar:
            
            # Set path to patient subfolders
            patientFolderPath = os.path.join(rater_dir,subject)
            BETFolderPath = os.path.join(patientFolderPath,'BETimages')
            normalised_path = os.path.join(patientFolderPath,'BET_normalised')
            sleep(0.03) #update progress bar
            bar()
            
            #Read brain and tumor masks
            brainMask = sitk.ReadImage(os.path.join(BETFolderPath,'MaskBrain.nii'), sitk.sitkUInt16)
            # tumorMask = sitk.ReadImage(os.path.join(BETFolderPath,'MaskTumor.nii'), sitk.sitkUInt16)
            
            #Read BET images
            rBV_noLC = sitk.ReadImage(os.path.join(BETFolderPath,'BET_rBV_noLC.nii'))
            rBV_LC = sitk.ReadImage(os.path.join(BETFolderPath,'BET_rBV_LC.nii'))
            rBF_LC = sitk.ReadImage(os.path.join(BETFolderPath,'BET_rBF_LC.nii'))
            MTT_LC = sitk.ReadImage(os.path.join(BETFolderPath,'BET_MTT_LC.nii'))
            
            #Read Normalised images
            n_rBV_noLC = sitk.ReadImage(os.path.join(normalised_path,'normalised_rBV_noLC.nii'))
            n_rBV_LC = sitk.ReadImage(os.path.join(normalised_path,'normalised_rBV_LC.nii'))
            n_rBF_LC = sitk.ReadImage(os.path.join(normalised_path,'normalised_rBF_LC.nii'))
            n_MTT_LC = sitk.ReadImage(os.path.join(normalised_path,'normalised_MTT_LC.nii'))
                       
            #Extract voxel values
            rBV_noLC_vxls = np.vstack(allVoxInt(rBV_noLC, brainMask))
            rBV_LC_vxls = np.vstack(allVoxInt(rBV_LC, brainMask))
            rBF_LC_vxls = np.vstack(allVoxInt(rBF_LC, brainMask))
            MTT_LC_vxls = np.vstack(allVoxInt(MTT_LC, brainMask))
            
            n_rBV_noLC_vxls = np.vstack(allVoxInt(n_rBV_noLC, brainMask))
            n_rBV_LC_vxls = np.vstack(allVoxInt(n_rBV_LC, brainMask))
            n_rBF_LC_vxls = np.vstack(allVoxInt(n_rBF_LC, brainMask))
            n_MTT_LC_vxls = np.vstack(allVoxInt(n_MTT_LC, brainMask))
            
            # Create a patient voxels intensity input dataframe and append to general

            rBV_noLC_d = {'Voxel_ID':[i for i in range(len(rBV_noLC_vxls))], 'Measurement':[idx for i in range(len(rBV_noLC_vxls))], 'rBV_noLC': [i[0] for i in rBV_noLC_vxls.tolist()]}
            rBV_noLC_df = pd.DataFrame(data=rBV_noLC_d)
            rBV_noLC_input_df=rBV_noLC_input_df.append(rBV_noLC_df)
            
            rBV_LC_d = {'Voxel_ID':[i for i in range(len(rBV_LC_vxls))], 'Measurement':[idx for i in range(len(rBV_LC_vxls))], 'rBV_LC': [i[0] for i in rBV_LC_vxls.tolist()]}
            rBV_LC_df = pd.DataFrame(data=rBV_LC_d)
            rBV_LC_input_df=rBV_LC_input_df.append(rBV_LC_df)
            
            rBF_LC_d = {'Voxel_ID':[i for i in range(len(rBF_LC_vxls))], 'Measurement':[idx for i in range(len(rBF_LC_vxls))], 'rBF_LC': [i[0] for i in rBF_LC_vxls.tolist()]}
            rBF_LC_df = pd.DataFrame(data=rBF_LC_d)
            rBF_LC_input_df=rBF_LC_input_df.append(rBF_LC_df)
            
            MTT_LC_d = {'Voxel_ID':[i for i in range(len(MTT_LC_vxls))], 'Measurement':[idx for i in range(len(MTT_LC_vxls))], 'MTT_LC': [i[0] for i in MTT_LC_vxls.tolist()]}
            MTT_LC_df = pd.DataFrame(data=MTT_LC_d)
            MTT_LC_input_df=MTT_LC_input_df.append(MTT_LC_df)
            
            
            
            n_rBV_noLC_d = {'Voxel_ID':[i for i in range(len(n_rBV_noLC_vxls))], 'Measurement':[idx for i in range(len(n_rBV_noLC_vxls))], 'n_rBV_noLC': [i[0] for i in n_rBV_noLC_vxls.tolist()]}
            n_rBV_noLC_df = pd.DataFrame(data=n_rBV_noLC_d)
            n_rBV_noLC_input_df=n_rBV_noLC_input_df.append(n_rBV_noLC_df)
            
            n_rBV_LC_d = {'Voxel_ID':[i for i in range(len(n_rBV_LC_vxls))], 'Measurement':[idx for i in range(len(n_rBV_LC_vxls))], 'n_rBV_LC': [i[0] for i in n_rBV_LC_vxls.tolist()]}
            n_rBV_LC_df = pd.DataFrame(data=n_rBV_LC_d)
            n_rBV_LC_input_df=n_rBV_LC_input_df.append(n_rBV_LC_df)
            
            n_rBF_LC_d = {'Voxel_ID':[i for i in range(len(n_rBF_LC_vxls))], 'Measurement':[idx for i in range(len(n_rBF_LC_vxls))], 'n_rBF_LC': [i[0] for i in n_rBF_LC_vxls.tolist()]}
            n_rBF_LC_df = pd.DataFrame(data=n_rBF_LC_d)
            n_rBF_LC_input_df=n_rBF_LC_input_df.append(n_rBF_LC_df)
            
            n_MTT_LC_d = {'Voxel_ID':[i for i in range(len(n_MTT_LC_vxls))], 'Measurement':[idx for i in range(len(n_MTT_LC_vxls))], 'n_MTT_LC': [i[0] for i in n_MTT_LC_vxls.tolist()]}
            n_MTT_LC_df = pd.DataFrame(data=n_MTT_LC_d)
            n_MTT_LC_input_df=n_MTT_LC_input_df.append(n_MTT_LC_df)
            

    #Calculate voxel-wise ICC in tumor
         
    print('Calculating rBV_noLC ICC...')
            
    ICC_rBV_noLC = pg.intraclass_corr(data=rBV_noLC_input_df, targets='Voxel_ID', raters='Measurement', ratings='rBV_noLC').round(3)  
    Results_ICC['rBV_noLC'] = Results_ICC['rBV_noLC'].append({'Subject_ID': subject,
                                                                  'Type':ICC_rBV_noLC.loc[4, 'Type'],
                                                                  'Description':ICC_rBV_noLC.loc[4, 'Description'], 
                                                                  'ICC':ICC_rBV_noLC.loc[4, 'ICC'], 
                                                                  'F':ICC_rBV_noLC.loc[4, 'F'],
                                                                  'df1':ICC_rBV_noLC.loc[4, 'df1'],
                                                                  'df2':ICC_rBV_noLC.loc[4, 'df2'],
                                                                  'pval':ICC_rBV_noLC.loc[4, 'pval'],
                                                                  'CI95%':ICC_rBV_noLC.loc[4, 'CI95%']
                                                                  }, ignore_index=True)
            
    print('Calculating rBV_LC ICC...')
            
    ICC_rBV_LC = pg.intraclass_corr(data=rBV_LC_input_df, targets='Voxel_ID', raters='Measurement', ratings='rBV_LC').round(3)  
    Results_ICC['rBV_LC'] = Results_ICC['rBV_LC'].append({'Subject_ID': subject,
                                                                  'Type':ICC_rBV_LC.loc[4, 'Type'],
                                                                  'Description':ICC_rBV_LC.loc[4, 'Description'], 
                                                                  'ICC':ICC_rBV_LC.loc[4, 'ICC'], 
                                                                  'F':ICC_rBV_LC.loc[4, 'F'],
                                                                  'df1':ICC_rBV_LC.loc[4, 'df1'],
                                                                  'df2':ICC_rBV_LC.loc[4, 'df2'],
                                                                  'pval':ICC_rBV_LC.loc[4, 'pval'],
                                                                  'CI95%':ICC_rBV_LC.loc[4, 'CI95%']
                                                                  }, ignore_index=True)
            
                    
    print('Calculating rBF_LC ICC...')
            
    ICC_rBF_LC = pg.intraclass_corr(data=rBF_LC_input_df, targets='Voxel_ID', raters='Measurement', ratings='rBF_LC').round(3)  
    Results_ICC['rBF_LC'] = Results_ICC['rBF_LC'].append({'Subject_ID': subject,
                                                                  'Type':ICC_rBF_LC.loc[4, 'Type'],
                                                                  'Description':ICC_rBF_LC.loc[4, 'Description'], 
                                                                  'ICC':ICC_rBF_LC.loc[4, 'ICC'], 
                                                                  'F':ICC_rBF_LC.loc[4, 'F'],
                                                                  'df1':ICC_rBF_LC.loc[4, 'df1'],
                                                                  'df2':ICC_rBF_LC.loc[4, 'df2'],
                                                                  'pval':ICC_rBF_LC.loc[4, 'pval'],
                                                                  'CI95%':ICC_rBF_LC.loc[4, 'CI95%']
                                                                  }, ignore_index=True)
                    
    print('Calculating MTT_LC ICC...')
            
    ICC_MTT_LC = pg.intraclass_corr(data=MTT_LC_input_df, targets='Voxel_ID', raters='Measurement', ratings='MTT_LC').round(3)  
    Results_ICC['MTT_LC'] = Results_ICC['MTT_LC'].append({'Subject_ID': subject,
                                                                  'Type':ICC_MTT_LC.loc[4, 'Type'],
                                                                  'Description':ICC_MTT_LC.loc[4, 'Description'], 
                                                                  'ICC':ICC_MTT_LC.loc[4, 'ICC'], 
                                                                  'F':ICC_MTT_LC.loc[4, 'F'],
                                                                  'df1':ICC_MTT_LC.loc[4, 'df1'],
                                                                  'df2':ICC_MTT_LC.loc[4, 'df2'],
                                                                  'pval':ICC_MTT_LC.loc[4, 'pval'],
                                                                  'CI95%':ICC_MTT_LC.loc[4, 'CI95%']
                                                                  }, ignore_index=True)
            
            
    print('Calculating n_rBV_noLC ICC...')
            
    ICC_n_rBV_noLC = pg.intraclass_corr(data=n_rBV_noLC_input_df, targets='Voxel_ID', raters='Measurement', ratings='n_rBV_noLC').round(3)  
    Results_ICC['n_rBV_noLC'] = Results_ICC['n_rBV_noLC'].append({'Subject_ID': subject,
                                                                  'Type':ICC_n_rBV_noLC.loc[4, 'Type'],
                                                                  'Description':ICC_n_rBV_noLC.loc[4, 'Description'], 
                                                                  'ICC':ICC_n_rBV_noLC.loc[4, 'ICC'], 
                                                                  'F':ICC_n_rBV_noLC.loc[4, 'F'],
                                                                  'df1':ICC_n_rBV_noLC.loc[4, 'df1'],
                                                                  'df2':ICC_n_rBV_noLC.loc[4, 'df2'],
                                                                  'pval':ICC_n_rBV_noLC.loc[4, 'pval'],
                                                                  'CI95%':ICC_n_rBV_noLC.loc[4, 'CI95%']
                                                                  }, ignore_index=True)
            
    print('Calculating n_rBV_LC ICC...')
            
    ICC_n_rBV_LC = pg.intraclass_corr(data=n_rBV_LC_input_df, targets='Voxel_ID', raters='Measurement', ratings='n_rBV_LC').round(3)  
    Results_ICC['n_rBV_LC'] = Results_ICC['n_rBV_LC'].append({'Subject_ID': subject,
                                                                  'Type':ICC_n_rBV_LC.loc[4, 'Type'],
                                                                  'Description':ICC_n_rBV_LC.loc[4, 'Description'], 
                                                                  'ICC':ICC_n_rBV_LC.loc[4, 'ICC'], 
                                                                  'F':ICC_n_rBV_LC.loc[4, 'F'],
                                                                  'df1':ICC_n_rBV_LC.loc[4, 'df1'],
                                                                  'df2':ICC_n_rBV_LC.loc[4, 'df2'],
                                                                  'pval':ICC_n_rBV_LC.loc[4, 'pval'],
                                                                  'CI95%':ICC_n_rBV_LC.loc[4, 'CI95%']
                                                                  }, ignore_index=True)
            
                    
    print('Calculating n_rBF_LC ICC...')
            
    ICC_n_rBF_LC = pg.intraclass_corr(data=n_rBF_LC_input_df, targets='Voxel_ID', raters='Measurement', ratings='n_rBF_LC').round(3)  
    Results_ICC['n_rBF_LC'] = Results_ICC['n_rBF_LC'].append({'Subject_ID': subject,
                                                                  'Type':ICC_n_rBF_LC.loc[4, 'Type'],
                                                                  'Description':ICC_n_rBF_LC.loc[4, 'Description'], 
                                                                  'ICC':ICC_n_rBF_LC.loc[4, 'ICC'], 
                                                                  'F':ICC_n_rBF_LC.loc[4, 'F'],
                                                                  'df1':ICC_n_rBF_LC.loc[4, 'df1'],
                                                                  'df2':ICC_n_rBF_LC.loc[4, 'df2'],
                                                                  'pval':ICC_n_rBF_LC.loc[4, 'pval'],
                                                                  'CI95%':ICC_n_rBF_LC.loc[4, 'CI95%']
                                                                  }, ignore_index=True)
                    
    print('Calculating n_MTT_LC ICC...')
            
    ICC_n_MTT_LC = pg.intraclass_corr(data=n_MTT_LC_input_df, targets='Voxel_ID', raters='Measurement', ratings='n_MTT_LC').round(3)  
    Results_ICC['n_MTT_LC'] = Results_ICC['n_MTT_LC'].append({'Subject_ID': subject,
                                                                  'Type':ICC_n_MTT_LC.loc[4, 'Type'],
                                                                  'Description':ICC_n_MTT_LC.loc[4, 'Description'], 
                                                                  'ICC':ICC_n_MTT_LC.loc[4, 'ICC'], 
                                                                  'F':ICC_n_MTT_LC.loc[4, 'F'],
                                                                  'df1':ICC_n_MTT_LC.loc[4, 'df1'],
                                                                  'df2':ICC_n_MTT_LC.loc[4, 'df2'],
                                                                  'pval':ICC_n_MTT_LC.loc[4, 'pval'],
                                                                  'CI95%':ICC_n_MTT_LC.loc[4, 'CI95%']
                                                                  }, ignore_index=True)
            
#%%Save all dataframes to excel files here
print('Save all results to excel files')

for name, df in Results_ICC.items():
    df.to_excel(IntraWriter, sheet_name=name, index=False)
IntraWriter.save()
