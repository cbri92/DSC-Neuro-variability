# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 14:52:43 2021

@author: Caterina Brighi

This script calculates repeatability metrics in a ROI between images of the same subject genrated by multiple raters.
Voxel-wise metrics include: 'Mean intensity in ROI rater1', 'Mean intensity in ROI rater2', 
                            'Mean intensity in ROI rater3', 'Between voxels variance', 'Within voxels variance', 
                            'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'
ROI-wise metrics include: 'BMS: Between subjects mean'
                          'WMS: Within subject mean'
                          'tSD: Total standard deviation'
                          'bSD: Between subjects standard deviation'
                          'wSD: Within subject standard deviation'
                          'RC: Repeatibility coefficient'
                          'RC Lower: Lower 95% CI'
                          'RC Upper: Upper 95% CI'
                          'wCV: Within subject coefficient of variation'
                          'ICC: intra-correlation coefficient'
"""


#%% Import functions 

import matplotlib.pyplot as plt
import math
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
import scipy
from scipy.stats.stats import pearsonr
from multiprocessing.pool import ThreadPool
from functools import partial
from radiomics import featureextractor
import six, numpy as np
from statistics import pvariance
from statistics import mean
from scipy.stats import normaltest
import seaborn as sns
from alive_progress import alive_bar
from time import sleep

def allVoxInt(image, roi):
    
    '''Returns a flatten array (in 2D) of the intensity values from all pixels in roi applied to an image.'''
    
    image = sitk.GetArrayFromImage(image)
    roi = sitk.GetArrayFromImage(roi)
    mask = np.array(roi,dtype='bool')
    mask = np.invert(mask)
    masked = image[~np.array(mask)]
    return masked

def Repeatibility_metrics(K, repeat1, repeat2, repeat3):
    '''This function returns a dictionary with the following repeatibility metrics:
        'BMS: Between subjects mean'
        'WMS: Within subject mean'
        'tSD: Total standard deviation'
        'bSD: Between subjects standard deviation'
        'wSD: Within subject standard deviation'
        'RC: Repeatibility coefficient'
        'RC Lower: Lower 95% CI'
        'RC Upper: Upper 95% CI'
        'wCV: Within subject coefficient of variation'
        'ICC: intra-correlation coefficient'
        '''
    mean_of_repeats = (repeat1+repeat2+repeat3)/K
    n=len(repeat1)
    overall_mean = (sum(repeat1)+sum(repeat2)+sum(repeat3))/(n*K)
    x=0
    for i in mean_of_repeats:
        x=x+(((i-overall_mean)**2)/n)
    BMS=K*x
    
    y=0
    for i,r1,r2,r3 in zip(mean_of_repeats, repeat1, repeat2, repeat3):
        y=y+(((r1-i)**2)/(n*(K-1))+((r2-i)**2)/(n*(K-1))+((r3-i)**2)/(n*(K-1)))
    WMS=y
    
    tSD = math.sqrt(((BMS+(K-1)*WMS)/K))
    bSD = math.sqrt((abs((BMS-WMS)/K)))
    wSD = math.sqrt(WMS)
    
    bsVar = bSD**2
    wsVar = wSD**2
    
    df=n*(K-1)
    #find Chi-Square critical value
    chi2_975=scipy.stats.chi2.ppf(0.975, df=df)
    chi2_025=scipy.stats.chi2.ppf(0.025, df=df)
    RC = 2.77*wSD
    RC_Lower = 2.77*math.sqrt((n*(K-1)*WMS)/chi2_975)
    RC_Upper = 2.77*math.sqrt((n*(K-1)*WMS)/chi2_025)
    wCV=wSD/overall_mean
    
    ICC = (tSD**2-wSD**2)/tSD**2
    
    return {'BMS':BMS, 'WMS':WMS, 'tSD':tSD, 'bSD':bSD, 'wSD':wSD, 'Within subject variance':wsVar, 'Between subject variance':bsVar, 'RC':RC, 'RC Lower':RC_Lower, 'RC Upper':RC_Upper, 'wCV':wCV, 'ICC':ICC}

    
    
    
#%% Set Working directory
        
# Path of the main folder 
DataSupraDir = 'Path to users supradirectory'

#Create empty dataframes to populate as going through the loop

rBV_noLC_df = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI rater1', 'Mean intensity in ROI rater2', 'Mean intensity in ROI rater3', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
rBV_LC_df = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI rater1', 'Mean intensity in ROI rater2', 'Mean intensity in ROI rater3', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
rBF_LC_df = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI rater1', 'Mean intensity in ROI rater2', 'Mean intensity in ROI rater3', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
MTT_LC_df = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI rater1', 'Mean intensity in ROI rater2', 'Mean intensity in ROI rater3', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])

n_rBV_noLC_df = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI rater1', 'Mean intensity in ROI rater2', 'Mean intensity in ROI rater3', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
n_rBV_LC_df = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI rater1', 'Mean intensity in ROI rater2', 'Mean intensity in ROI rater3', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
n_rBF_LC_df = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI rater1', 'Mean intensity in ROI rater2', 'Mean intensity in ROI rater3', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
n_MTT_LC_df = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI rater1', 'Mean intensity in ROI rater2', 'Mean intensity in ROI rater3', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])


Repeatibility_statsdf = pd.DataFrame(columns=['Image type', 'WMS', 'BMS','wSD', 'bSD', 'tSD', 'RC', 'RC Upper', 'RC Lower', 'wCV', 'ICC'])


#Create an excel spreadsheet with the stats results
Group_stats_results = pd.ExcelWriter(DataSupraDir +'Brain_repeatability_results.xlsx')

reters_path = [ f.path for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of the paths to the raters directories
raters_name = [ f.name for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of raters names

n_raters = len(raters_name) #Total number of raters


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


    # Loop through each rater 
    for idx, rater in enumerate(raters_name):
        print(idx, rater)
        rater_dir = DataSupraDir+rater
        
    
        print('Calculating stats for rater:',rater)
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
           
            #Append voxels intensity in ROI for this rater to dictionary
            rBV_noLC_dic['Voxel intensity in ROI rater'+str(idx+1)]=rBV_noLC_vxls
            rBV_LC_dic['Voxel intensity in ROI rater'+str(idx+1)]=rBV_LC_vxls
            rBF_LC_dic['Voxel intensity in ROI rater'+str(idx+1)]=rBF_LC_vxls
            MTT_LC_dic['Voxel intensity in ROI rater'+str(idx+1)]=MTT_LC_vxls
            
            n_rBV_noLC_dic['Voxel intensity in ROI rater'+str(idx+1)]=n_rBV_noLC_vxls
            n_rBV_LC_dic['Voxel intensity in ROI rater'+str(idx+1)]=n_rBV_LC_vxls
            n_rBF_LC_dic['Voxel intensity in ROI rater'+str(idx+1)]=n_rBF_LC_vxls
            n_MTT_LC_dic['Voxel intensity in ROI rater'+str(idx+1)]=n_MTT_LC_vxls

            #Calculate the mean intensity in ROI
            rBV_noLC_mean=rBV_noLC_vxls.mean()
            rBV_LC_mean=rBV_LC_vxls.mean()
            rBF_LC_mean=rBF_LC_vxls.mean()
            MTT_LC_mean=MTT_LC_vxls.mean()
            
            n_rBV_noLC_mean=rBV_noLC_vxls.mean()
            n_rBV_LC_mean=rBV_LC_vxls.mean()
            n_rBF_LC_mean=rBF_LC_vxls.mean()
            n_MTT_LC_mean=MTT_LC_vxls.mean()
            
            #Append Mean intensity in ROI for this rater to dictionary
            rBV_noLC_dic['Mean intensity in ROI rater'+str(idx+1)]=rBV_noLC_mean
            rBV_LC_dic['Mean intensity in ROI rater'+str(idx+1)]=rBV_LC_mean
            rBF_LC_dic['Mean intensity in ROI rater'+str(idx+1)]=rBF_LC_mean
            MTT_LC_dic['Mean intensity in ROI rater'+str(idx+1)]=MTT_LC_mean
            
            n_rBV_noLC_dic['Mean intensity in ROI rater'+str(idx+1)]=n_rBV_noLC_mean
            n_rBV_LC_dic['Mean intensity in ROI rater'+str(idx+1)]=n_rBV_LC_mean
            n_rBF_LC_dic['Mean intensity in ROI rater'+str(idx+1)]=n_rBF_LC_mean
            n_MTT_LC_dic['Mean intensity in ROI rater'+str(idx+1)]=n_MTT_LC_mean



   
#%% Calculate ROI-based group repeatibility statistics

    #Calculate repeatability metrics within ROI
        
    print('Calculate repeatability metrics within ROI for subject '+subject)
        
    rBV_noLC_mtrcs = Repeatibility_metrics(3, rBV_noLC_dic.get('Voxel intensity in ROI rater1'), rBV_noLC_dic.get('Voxel intensity in ROI rater2'), rBV_noLC_dic.get('Voxel intensity in ROI rater3'))
    rBV_LC_mtrcs = Repeatibility_metrics(3, rBV_LC_dic.get('Voxel intensity in ROI rater1'), rBV_LC_dic.get('Voxel intensity in ROI rater2'), rBV_LC_dic.get('Voxel intensity in ROI rater3'))
    rBF_LC_mtrcs = Repeatibility_metrics(3, rBF_LC_dic.get('Voxel intensity in ROI rater1'), rBF_LC_dic.get('Voxel intensity in ROI rater2'), rBF_LC_dic.get('Voxel intensity in ROI rater3'))
    MTT_LC_mtrcs = Repeatibility_metrics(3, MTT_LC_dic.get('Voxel intensity in ROI rater1'), MTT_LC_dic.get('Voxel intensity in ROI rater2'), MTT_LC_dic.get('Voxel intensity in ROI rater3'))
        
    n_rBV_noLC_mtrcs = Repeatibility_metrics(3, n_rBV_noLC_dic.get('Voxel intensity in ROI rater1'), n_rBV_noLC_dic.get('Voxel intensity in ROI rater2'), n_rBV_noLC_dic.get('Voxel intensity in ROI rater3'))
    n_rBV_LC_mtrcs = Repeatibility_metrics(3, n_rBV_LC_dic.get('Voxel intensity in ROI rater1'), n_rBV_LC_dic.get('Voxel intensity in ROI rater2'), n_rBV_LC_dic.get('Voxel intensity in ROI rater3'))
    n_rBF_LC_mtrcs = Repeatibility_metrics(3, n_rBF_LC_dic.get('Voxel intensity in ROI rater1'), n_rBF_LC_dic.get('Voxel intensity in ROI rater2'), n_rBF_LC_dic.get('Voxel intensity in ROI rater3'))
    n_MTT_LC_mtrcs = Repeatibility_metrics(3, n_MTT_LC_dic.get('Voxel intensity in ROI rater1'), n_MTT_LC_dic.get('Voxel intensity in ROI rater2'), n_MTT_LC_dic.get('Voxel intensity in ROI rater3'))

#%% Export subjects stats into excel spreadsheet

    print('Export subjects '+subject+' stats into excel spreadsheet')
        
    rBV_noLC_stats = {'Subject_ID':subject, 'Mean intensity in ROI rater1':rBV_noLC_dic.get('Mean intensity in ROI rater1'), 'Mean intensity in ROI rater2':rBV_noLC_dic.get('Mean intensity in ROI rater2'), 'Mean intensity in ROI rater3':rBV_noLC_dic.get('Mean intensity in ROI rater3'),'Between voxels variance':rBV_noLC_mtrcs.get('Between subject variance'), 'Within voxels variance':rBV_noLC_mtrcs.get('Within subject variance'), 'ROI ICC':rBV_noLC_mtrcs.get('ICC'), 'Within voxel CoV':rBV_noLC_mtrcs.get('wCV'), 'Repeatability coefficient':rBV_noLC_mtrcs.get('RC'), 'RC Upper':rBV_noLC_mtrcs.get('RC Upper'), 'RC Lower':rBV_noLC_mtrcs.get('RC Lower')}
    rBV_LC_stats = {'Subject_ID':subject, 'Mean intensity in ROI rater1':rBV_LC_dic.get('Mean intensity in ROI rater1'), 'Mean intensity in ROI rater2':rBV_LC_dic.get('Mean intensity in ROI rater2'), 'Mean intensity in ROI rater3':rBV_LC_dic.get('Mean intensity in ROI rater3'),'Between voxels variance':rBV_LC_mtrcs.get('Between subject variance'), 'Within voxels variance':rBV_LC_mtrcs.get('Within subject variance'), 'ROI ICC':rBV_LC_mtrcs.get('ICC'), 'Within voxel CoV':rBV_LC_mtrcs.get('wCV'), 'Repeatability coefficient':rBV_LC_mtrcs.get('RC'), 'RC Upper':rBV_LC_mtrcs.get('RC Upper'), 'RC Lower':rBV_LC_mtrcs.get('RC Lower')}
    rBF_LC_stats = {'Subject_ID':subject, 'Mean intensity in ROI rater1':rBF_LC_dic.get('Mean intensity in ROI rater1'), 'Mean intensity in ROI rater2':rBF_LC_dic.get('Mean intensity in ROI rater2'), 'Mean intensity in ROI rater3':rBF_LC_dic.get('Mean intensity in ROI rater3'),'Between voxels variance':rBF_LC_mtrcs.get('Between subject variance'), 'Within voxels variance':rBF_LC_mtrcs.get('Within subject variance'), 'ROI ICC':rBF_LC_mtrcs.get('ICC'), 'Within voxel CoV':rBF_LC_mtrcs.get('wCV'), 'Repeatability coefficient':rBF_LC_mtrcs.get('RC'), 'RC Upper':rBF_LC_mtrcs.get('RC Upper'), 'RC Lower':rBF_LC_mtrcs.get('RC Lower')}
    MTT_LC_stats = {'Subject_ID':subject, 'Mean intensity in ROI rater1':MTT_LC_dic.get('Mean intensity in ROI rater1'), 'Mean intensity in ROI rater2':MTT_LC_dic.get('Mean intensity in ROI rater2'), 'Mean intensity in ROI rater3':MTT_LC_dic.get('Mean intensity in ROI rater3'),'Between voxels variance':MTT_LC_mtrcs.get('Between subject variance'), 'Within voxels variance':MTT_LC_mtrcs.get('Within subject variance'), 'ROI ICC':MTT_LC_mtrcs.get('ICC'), 'Within voxel CoV':MTT_LC_mtrcs.get('wCV'), 'Repeatability coefficient':MTT_LC_mtrcs.get('RC'), 'RC Upper':MTT_LC_mtrcs.get('RC Upper'), 'RC Lower':MTT_LC_mtrcs.get('RC Lower')}
        
    n_rBV_noLC_stats = {'Subject_ID':subject, 'Mean intensity in ROI rater1':n_rBV_noLC_dic.get('Mean intensity in ROI rater1'), 'Mean intensity in ROI rater2':n_rBV_noLC_dic.get('Mean intensity in ROI rater2'), 'Mean intensity in ROI rater3':n_rBV_noLC_dic.get('Mean intensity in ROI rater3'),'Between voxels variance':n_rBV_noLC_mtrcs.get('Between subject variance'), 'Within voxels variance':n_rBV_noLC_mtrcs.get('Within subject variance'), 'ROI ICC':n_rBV_noLC_mtrcs.get('ICC'), 'Within voxel CoV':n_rBV_noLC_mtrcs.get('wCV'), 'Repeatability coefficient':n_rBV_noLC_mtrcs.get('RC'), 'RC Upper':n_rBV_noLC_mtrcs.get('RC Upper'), 'RC Lower':n_rBV_noLC_mtrcs.get('RC Lower')}
    n_rBV_LC_stats = {'Subject_ID':subject, 'Mean intensity in ROI rater1':n_rBV_LC_dic.get('Mean intensity in ROI rater1'), 'Mean intensity in ROI rater2':n_rBV_LC_dic.get('Mean intensity in ROI rater2'), 'Mean intensity in ROI rater3':n_rBV_LC_dic.get('Mean intensity in ROI rater3'),'Between voxels variance':n_rBV_LC_mtrcs.get('Between subject variance'), 'Within voxels variance':n_rBV_LC_mtrcs.get('Within subject variance'), 'ROI ICC':n_rBV_LC_mtrcs.get('ICC'), 'Within voxel CoV':n_rBV_LC_mtrcs.get('wCV'), 'Repeatability coefficient':n_rBV_LC_mtrcs.get('RC'), 'RC Upper':n_rBV_LC_mtrcs.get('RC Upper'), 'RC Lower':n_rBV_LC_mtrcs.get('RC Lower')}
    n_rBF_LC_stats = {'Subject_ID':subject, 'Mean intensity in ROI rater1':n_rBF_LC_dic.get('Mean intensity in ROI rater1'), 'Mean intensity in ROI rater2':n_rBF_LC_dic.get('Mean intensity in ROI rater2'), 'Mean intensity in ROI rater3':n_rBF_LC_dic.get('Mean intensity in ROI rater3'),'Between voxels variance':n_rBF_LC_mtrcs.get('Between subject variance'), 'Within voxels variance':n_rBF_LC_mtrcs.get('Within subject variance'), 'ROI ICC':n_rBF_LC_mtrcs.get('ICC'), 'Within voxel CoV':n_rBF_LC_mtrcs.get('wCV'), 'Repeatability coefficient':n_rBF_LC_mtrcs.get('RC'), 'RC Upper':n_rBF_LC_mtrcs.get('RC Upper'), 'RC Lower':n_rBF_LC_mtrcs.get('RC Lower')}
    n_MTT_LC_stats = {'Subject_ID':subject, 'Mean intensity in ROI rater1':n_MTT_LC_dic.get('Mean intensity in ROI rater1'), 'Mean intensity in ROI rater2':n_MTT_LC_dic.get('Mean intensity in ROI rater2'), 'Mean intensity in ROI rater3':n_MTT_LC_dic.get('Mean intensity in ROI rater3'),'Between voxels variance':n_MTT_LC_mtrcs.get('Between subject variance'), 'Within voxels variance':n_MTT_LC_mtrcs.get('Within subject variance'), 'ROI ICC':n_MTT_LC_mtrcs.get('ICC'), 'Within voxel CoV':n_MTT_LC_mtrcs.get('wCV'), 'Repeatability coefficient':n_MTT_LC_mtrcs.get('RC'), 'RC Upper':n_MTT_LC_mtrcs.get('RC Upper'), 'RC Lower':n_MTT_LC_mtrcs.get('RC Lower')}
        
        
    rBV_noLC_df = rBV_noLC_df.append(rBV_noLC_stats, ignore_index=True)
    rBV_LC_df = rBV_LC_df.append(rBV_LC_stats, ignore_index=True)
    rBF_LC_df = rBF_LC_df.append(rBF_LC_stats, ignore_index=True)
    MTT_LC_df = MTT_LC_df.append(MTT_LC_stats, ignore_index=True)
        
    n_rBV_noLC_df = n_rBV_noLC_df.append(n_rBV_noLC_stats, ignore_index=True)
    n_rBV_LC_df = n_rBV_LC_df.append(n_rBV_LC_stats, ignore_index=True)
    n_rBF_LC_df = n_rBF_LC_df.append(n_rBF_LC_stats, ignore_index=True)
    n_MTT_LC_df = n_MTT_LC_df.append(n_MTT_LC_stats, ignore_index=True)
    
    
    
#%%Calculate group arrays

#Calculate rebeatability metrics

print('Calculate repeatability metrics on entire group')

rBV_noLC_ROIrepeatibility=Repeatibility_metrics(3, rBV_noLC_df['Mean intensity in ROI rater1'], rBV_noLC_df['Mean intensity in ROI rater2'], rBV_noLC_df['Mean intensity in ROI rater3'])
rBV_LC_ROIrepeatibility=Repeatibility_metrics(3, rBV_LC_df['Mean intensity in ROI rater1'], rBV_LC_df['Mean intensity in ROI rater2'], rBV_LC_df['Mean intensity in ROI rater3'])
rBF_LC_ROIrepeatibility=Repeatibility_metrics(3, rBF_LC_df['Mean intensity in ROI rater1'], rBF_LC_df['Mean intensity in ROI rater2'], rBF_LC_df['Mean intensity in ROI rater3'])
MTT_LC_ROIrepeatibility=Repeatibility_metrics(3, MTT_LC_df['Mean intensity in ROI rater1'], MTT_LC_df['Mean intensity in ROI rater2'], MTT_LC_df['Mean intensity in ROI rater3'])

n_rBV_noLC_ROIrepeatibility=Repeatibility_metrics(3, n_rBV_noLC_df['Mean intensity in ROI rater1'], n_rBV_noLC_df['Mean intensity in ROI rater2'], n_rBV_noLC_df['Mean intensity in ROI rater3'])
n_rBV_LC_ROIrepeatibility=Repeatibility_metrics(3, n_rBV_LC_df['Mean intensity in ROI rater1'], n_rBV_LC_df['Mean intensity in ROI rater2'], n_rBV_LC_df['Mean intensity in ROI rater3'])
n_rBF_LC_ROIrepeatibility=Repeatibility_metrics(3, n_rBF_LC_df['Mean intensity in ROI rater1'], n_rBF_LC_df['Mean intensity in ROI rater2'], n_rBF_LC_df['Mean intensity in ROI rater3'])
n_MTT_LC_ROIrepeatibility=Repeatibility_metrics(3, n_MTT_LC_df['Mean intensity in ROI rater1'], n_MTT_LC_df['Mean intensity in ROI rater2'], n_MTT_LC_df['Mean intensity in ROI rater3'])

#Add image modality to dictionary
rBV_noLC_ROIrepeatibility['Image type']= 'rBV_noLC'
rBV_LC_ROIrepeatibility['Image type']= 'rBV_LC'
rBF_LC_ROIrepeatibility['Image type']= 'rBF_LC'
MTT_LC_ROIrepeatibility['Image type']= 'MTT_LC'

n_rBV_noLC_ROIrepeatibility['Image type']= 'n_rBV_noLC'
n_rBV_LC_ROIrepeatibility['Image type']= 'n_rBV_LC'
n_rBF_LC_ROIrepeatibility['Image type']= 'n_rBF_LC'
n_MTT_LC_ROIrepeatibility['Image type']= 'n_MTT_LC'



#Append metrics to repeatibility sheet
Repeatibility_statsdf = Repeatibility_statsdf.append(rBV_noLC_ROIrepeatibility, ignore_index=True)
Repeatibility_statsdf = Repeatibility_statsdf.append(rBV_LC_ROIrepeatibility, ignore_index=True)
Repeatibility_statsdf = Repeatibility_statsdf.append(rBF_LC_ROIrepeatibility, ignore_index=True)
Repeatibility_statsdf = Repeatibility_statsdf.append(MTT_LC_ROIrepeatibility, ignore_index=True)

Repeatibility_statsdf = Repeatibility_statsdf.append(n_rBV_noLC_ROIrepeatibility, ignore_index=True)
Repeatibility_statsdf = Repeatibility_statsdf.append(n_rBV_LC_ROIrepeatibility, ignore_index=True)
Repeatibility_statsdf = Repeatibility_statsdf.append(n_rBF_LC_ROIrepeatibility, ignore_index=True)
Repeatibility_statsdf = Repeatibility_statsdf.append(n_MTT_LC_ROIrepeatibility, ignore_index=True)



#%%Save all dataframes to excel files here
print('Save all results to excel files')

rBV_noLC_df.to_excel(Group_stats_results, sheet_name='rBV_noLC', index=False)
rBV_LC_df.to_excel(Group_stats_results, sheet_name='rBV_LC', index=False)
rBF_LC_df.to_excel(Group_stats_results, sheet_name='rBF_LC', index=False)
MTT_LC_df.to_excel(Group_stats_results, sheet_name='MTT_LC', index=False)

n_rBV_noLC_df.to_excel(Group_stats_results, sheet_name='n_rBV_noLC', index=False)
n_rBV_LC_df.to_excel(Group_stats_results, sheet_name='n_rBV_LC', index=False)
n_rBF_LC_df.to_excel(Group_stats_results, sheet_name='n_rBF_LC', index=False)
n_MTT_LC_df.to_excel(Group_stats_results, sheet_name='n_MTT_LC', index=False)


Repeatibility_statsdf.to_excel(Group_stats_results, sheet_name='Group Repeatibility metrics', index=False)

Group_stats_results.save()   

#%%Plot results of ICC

print('Plot results of group ICC')

rBV_noLC_ICC = rBV_noLC_df['ROI ICC'].to_numpy()
rBV_LC_ICC = rBV_LC_df['ROI ICC'].to_numpy()
rBF_LC_ICC = rBF_LC_df['ROI ICC'].to_numpy()
MTT_LC_ICC = MTT_LC_df['ROI ICC'].to_numpy()

n_rBV_noLC_ICC = n_rBV_noLC_df['ROI ICC'].to_numpy()
n_rBV_LC_ICC = n_rBV_LC_df['ROI ICC'].to_numpy()
n_rBF_LC_ICC = n_rBF_LC_df['ROI ICC'].to_numpy()
n_MTT_LC_ICC = n_MTT_LC_df['ROI ICC'].to_numpy()


ICC_data = [rBV_noLC_ICC, n_rBV_noLC_ICC, rBV_LC_ICC, n_rBV_LC_ICC, rBF_LC_ICC, n_rBF_LC_ICC, MTT_LC_ICC, n_MTT_LC_ICC]
ICC_labels = ['rBV_noLC', 'n_rBV_noLC', 'rBV_LC', 'n_rBV_LC', 'rBF_LC', 'n_rBF_LC', 'MTT_LC', 'n_MTT_LC']

fig, ax = plt.subplots(figsize=(10, 6))
ax.boxplot(ICC_data, labels=ICC_labels)
ax.set_ylim(0, 1)
plt.xlabel('Image modality', fontsize=20)
plt.ylabel('ICC', fontsize=20)
plt.setp(ax.get_xticklabels(), rotation = 45, fontsize=12)
plt.setp(ax.get_yticklabels(), fontsize=16)
fig.tight_layout(pad=2.0)
plt.savefig(DataSupraDir+'Brain_ICC.png')     
plt.close()



print('Plot results of group wCV')

rBV_noLC_wCV = rBV_noLC_df['Within voxel CoV'].to_numpy()
rBV_LC_wCV = rBV_LC_df['Within voxel CoV'].to_numpy()
rBF_LC_wCV = rBF_LC_df['Within voxel CoV'].to_numpy()
MTT_LC_wCV = MTT_LC_df['Within voxel CoV'].to_numpy()

n_rBV_noLC_wCV = n_rBV_noLC_df['Within voxel CoV'].to_numpy()
n_rBV_LC_wCV = n_rBV_LC_df['Within voxel CoV'].to_numpy()
n_rBF_LC_wCV = n_rBF_LC_df['Within voxel CoV'].to_numpy()
n_MTT_LC_wCV = n_MTT_LC_df['Within voxel CoV'].to_numpy()

wCV_data = [rBV_noLC_wCV, n_rBV_noLC_wCV, rBV_LC_wCV, n_rBV_LC_wCV, rBF_LC_wCV, n_rBF_LC_wCV, MTT_LC_wCV, n_MTT_LC_wCV]
wCV_labels = ['rBV_noLC', 'n_rBV_noLC', 'rBV_LC', 'n_rBV_LC', 'rBF_LC', 'n_rBF_LC', 'MTT_LC', 'n_MTT_LC']

fig, ax = plt.subplots(figsize=(10, 6))
ax.boxplot(wCV_data, labels=wCV_labels)
ax.set_ylim(0, 1)
plt.xlabel('Image modality', fontsize=20)
plt.ylabel('Within voxel CoV', fontsize=20)
plt.setp(ax.get_xticklabels(), rotation = 45, fontsize=12)
plt.setp(ax.get_yticklabels(), fontsize=16)
fig.tight_layout(pad=2.0)
plt.savefig(DataSupraDir+'Brain_wCV.png')     
plt.close()
