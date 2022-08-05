# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 14:37:53 2022

@author: Caterina Brighi
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
from ImageSimilarityMetricsFunctions import *

    
def generate_mask(image, roi):
    
    '''Returns the masked image of the roi applied to the image.
    Remember to save the masked image file after applying this function.'''
    
    masking = sitk.MaskImageFilter() 
    mask = masking.Execute(image, roi)
    return mask

    
#%% Set Working directory
        
# Path of the main folder 
DataSupraDir = 'Path to users supradirectory'


reters_path = [ f.path for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of the paths to the raters directories
raters_name = [ f.name for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of raters names

for i in raters_name:
    raters_name = [ f.name for f in os.scandir(DataSupraDir) if f.is_dir() ]
    raters_name.remove(i)

    n_raters = len(raters_name) #Total number of raters
    
    Results_simil = {'rBV_noLC': pd.DataFrame(columns=['Subject_ID','MSE', 'RMSE','NRMSE', 'PSNR', 'SSIM']), 'rBV_LC': pd.DataFrame(columns=['Subject_ID','MSE', 'RMSE','NRMSE', 'PSNR', 'SSIM']), 'rBF_LC': pd.DataFrame(columns=['Subject_ID','MSE', 'RMSE','NRMSE', 'PSNR', 'SSIM']), 'MTT_LC': pd.DataFrame(columns=['Subject_ID','MSE', 'RMSE','NRMSE', 'PSNR', 'SSIM']), 'n_rBV_noLC': pd.DataFrame(columns=['Subject_ID','MSE', 'RMSE','NRMSE', 'PSNR', 'SSIM']), 'n_rBV_LC': pd.DataFrame(columns=['Subject_ID','MSE', 'RMSE','NRMSE', 'PSNR', 'SSIM']), 'n_rBF_LC': pd.DataFrame(columns=['Subject_ID','MSE', 'RMSE','NRMSE', 'PSNR', 'SSIM']), 'n_MTT_LC': pd.DataFrame(columns=['Subject_ID','MSE', 'RMSE','NRMSE', 'PSNR', 'SSIM'])}
    IntraWriter = pd.ExcelWriter(DataSupraDir +raters_name[0]+'-'+raters_name[1]+'_Brain_Similarity_Results.xlsx', engine='xlsxwriter')
    
    # Extract number of patient folders as a list 
    subjs_name = [ f.name for f in os.scandir(DataSupraDir+raters_name[0]) if f.is_dir() ] #Create a list of raters names
    n_subjs = len(subjs_name) #Total number of subjects
    
    #%%Create a for loop to perform image analysis on each rater sequentially
    
    for subject in subjs_name:
        
        print('Analysing patient '+subject)
        
        rBV_noLC_dic = {'Subject_ID':subject}
        rBV_LC_dic = {'Subject_ID':subject}
        rBF_LC_dic = {'Subject_ID':subject}
        MTT_LC_dic = {'Subject_ID':subject}
        
        n_rBV_noLC_dic = {'Subject_ID':subject}
        n_rBV_LC_dic = {'Subject_ID':subject}
        n_rBF_LC_dic = {'Subject_ID':subject}
        n_MTT_LC_dic = {'Subject_ID':subject}
        
        brainMask_dic = {'Subject_ID':subject}
    
    
        # Loop through each rater 
        for idx, rater in enumerate(raters_name):
            print(idx, rater)
            rater_dir = DataSupraDir+rater
            
                
            # Set path to patient subfolders
            patientFolderPath = os.path.join(rater_dir,subject)
            BETFolderPath = os.path.join(patientFolderPath,'BETimages')
            normalised_path = os.path.join(patientFolderPath,'BET_normalised')
                
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
                
            #Append image this rater to dictionary
            rBV_noLC_dic['Image for rater'+str(idx+1)]=rBV_noLC
            rBV_LC_dic['Image for rater'+str(idx+1)]=rBV_LC
            rBF_LC_dic['Image for rater'+str(idx+1)]=rBF_LC
            MTT_LC_dic['Image for rater'+str(idx+1)]=MTT_LC
                
            n_rBV_noLC_dic['Image for rater'+str(idx+1)]=n_rBV_noLC
            n_rBV_LC_dic['Image for rater'+str(idx+1)]=n_rBV_LC
            n_rBF_LC_dic['Image for rater'+str(idx+1)]=n_rBF_LC
            n_MTT_LC_dic['Image for rater'+str(idx+1)]=n_MTT_LC
                
            brainMask_dic['Image for rater'+str(idx+1)]=brainMask
                
 
#%% Calculate similarity metrics within ROI
        
        print('Calculate similarity metrics within ROI for subject '+subject)
            
        rBV_noLC_mtrcs = similarity_metrics_inROI(rBV_noLC_dic.get('Image for rater1'), rBV_noLC_dic.get('Image for rater2'),  brainMask_dic.get('Image for rater1'))
        rBV_LC_mtrcs = similarity_metrics_inROI(rBV_LC_dic.get('Image for rater1'), rBV_LC_dic.get('Image for rater2'),  brainMask_dic.get('Image for rater1'))
        rBF_LC_mtrcs = similarity_metrics_inROI(rBF_LC_dic.get('Image for rater1'), rBF_LC_dic.get('Image for rater2'),  brainMask_dic.get('Image for rater1'))
        MTT_LC_mtrcs = similarity_metrics_inROI(MTT_LC_dic.get('Image for rater1'), MTT_LC_dic.get('Image for rater2'),  brainMask_dic.get('Image for rater1'))
            
        n_rBV_noLC_mtrcs = similarity_metrics_inROI(n_rBV_noLC_dic.get('Image for rater1'), n_rBV_noLC_dic.get('Image for rater2'),  brainMask_dic.get('Image for rater1'))
        n_rBV_LC_mtrcs = similarity_metrics_inROI(n_rBV_LC_dic.get('Image for rater1'), n_rBV_LC_dic.get('Image for rater2'),  brainMask_dic.get('Image for rater1'))
        n_rBF_LC_mtrcs = similarity_metrics_inROI(n_rBF_LC_dic.get('Image for rater1'), n_rBF_LC_dic.get('Image for rater2'),  brainMask_dic.get('Image for rater1'))
        n_MTT_LC_mtrcs = similarity_metrics_inROI(n_MTT_LC_dic.get('Image for rater1'), n_MTT_LC_dic.get('Image for rater2'),  brainMask_dic.get('Image for rater1'))
    
    #%% Export subjects stats into excel spreadsheet
    
        print('Export subjects '+subject+' metrics into excel spreadsheet')
            

        Results_simil['rBV_noLC'] = Results_simil['rBV_noLC'].append({'Subject_ID': subject,'MSE':rBV_noLC_mtrcs.get('MSE'), 'RMSE':rBV_noLC_mtrcs.get('RMSE'),'NRMSE':rBV_noLC_mtrcs.get('NRMSE'), 'PSNR':rBV_noLC_mtrcs.get('PSNR'), 'SSIM':rBV_noLC_mtrcs.get('SSIM')}, ignore_index=True)
        Results_simil['rBV_LC'] = Results_simil['rBV_LC'].append({'Subject_ID': subject,'MSE':rBV_LC_mtrcs.get('MSE'), 'RMSE':rBV_LC_mtrcs.get('RMSE'),'NRMSE':rBV_LC_mtrcs.get('NRMSE'), 'PSNR':rBV_LC_mtrcs.get('PSNR'), 'SSIM':rBV_LC_mtrcs.get('SSIM')}, ignore_index=True)
        Results_simil['rBF_LC'] = Results_simil['rBF_LC'].append({'Subject_ID': subject,'MSE':rBF_LC_mtrcs.get('MSE'), 'RMSE':rBF_LC_mtrcs.get('RMSE'),'NRMSE':rBF_LC_mtrcs.get('NRMSE'), 'PSNR':rBF_LC_mtrcs.get('PSNR'), 'SSIM':rBF_LC_mtrcs.get('SSIM')}, ignore_index=True)
        Results_simil['MTT_LC'] = Results_simil['MTT_LC'].append({'Subject_ID': subject,'MSE':MTT_LC_mtrcs.get('MSE'), 'RMSE':MTT_LC_mtrcs.get('RMSE'),'NRMSE':MTT_LC_mtrcs.get('NRMSE'), 'PSNR':MTT_LC_mtrcs.get('PSNR'), 'SSIM':MTT_LC_mtrcs.get('SSIM')}, ignore_index=True)
                
        Results_simil['n_rBV_noLC'] = Results_simil['n_rBV_noLC'].append({'Subject_ID': subject,'MSE':n_rBV_noLC_mtrcs.get('MSE'), 'RMSE':n_rBV_noLC_mtrcs.get('RMSE'),'NRMSE':n_rBV_noLC_mtrcs.get('NRMSE'), 'PSNR':n_rBV_noLC_mtrcs.get('PSNR'), 'SSIM':n_rBV_noLC_mtrcs.get('SSIM')}, ignore_index=True)
        Results_simil['n_rBV_LC'] = Results_simil['n_rBV_LC'].append({'Subject_ID': subject,'MSE':n_rBV_LC_mtrcs.get('MSE'), 'RMSE':n_rBV_LC_mtrcs.get('RMSE'),'NRMSE':n_rBV_LC_mtrcs.get('NRMSE'), 'PSNR':n_rBV_LC_mtrcs.get('PSNR'), 'SSIM':n_rBV_LC_mtrcs.get('SSIM')}, ignore_index=True)
        Results_simil['n_rBF_LC'] = Results_simil['n_rBF_LC'].append({'Subject_ID': subject,'MSE':n_rBF_LC_mtrcs.get('MSE'), 'RMSE':n_rBF_LC_mtrcs.get('RMSE'),'NRMSE':n_rBF_LC_mtrcs.get('NRMSE'), 'PSNR':n_rBF_LC_mtrcs.get('PSNR'), 'SSIM':n_rBF_LC_mtrcs.get('SSIM')}, ignore_index=True)
        Results_simil['n_MTT_LC'] = Results_simil['n_MTT_LC'].append({'Subject_ID': subject,'MSE':n_MTT_LC_mtrcs.get('MSE'), 'RMSE':n_MTT_LC_mtrcs.get('RMSE'),'NRMSE':n_MTT_LC_mtrcs.get('NRMSE'), 'PSNR':n_MTT_LC_mtrcs.get('PSNR'), 'SSIM':n_MTT_LC_mtrcs.get('SSIM')}, ignore_index=True)
                
    #%%Save all dataframes to excel files here
    print('Save all results to excel files')
    
    for name, df in Results_simil.items():
        df.to_excel(IntraWriter, sheet_name=name, index=False)
    IntraWriter.save()
