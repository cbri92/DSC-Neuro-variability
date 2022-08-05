# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 14:21:54 2022

@author: Caterina Brighi and Athena Kam
"""

# This script will loop through the RDS to create normalised maps for everything in the BETimages folder
import SimpleITK as sitk
import os 
import numpy as np

# Progress bar 
from alive_progress import alive_bar
from time import sleep

# Define masking function 


def normalise(image_path, roi):
    '''This function takes in the perfusion map and the mask as parameters. It returns the mean value of the ROI.'''
    #Read the image 
    image = sitk.ReadImage(image_path)
    stats = sitk.LabelIntensityStatisticsImageFilter()    
    stats.Execute(roi, image)    
    meanVal = stats.GetMean(1) 
    normalised = image/meanVal
    return normalised


# Path of the main folder 
DataSupraDir = 'Path to users supradirectory'

reters_path = [ f.path for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of the paths to the raters directories
raters_name = [ f.name for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of raters names
raters_name.remove("AthenaKam")

n_raters = len(raters_name) #Total number of raters

#%%Create a for loop to perform image analysis on each rater sequentially

for rater in raters_name:
    rater_dir = DataSupraDir+rater

    # Extract number of patient folders as a list 
    patientFolderList = [f for f in os.listdir(rater_dir) if os.path.isdir(os.path.join(rater_dir,f))]
    
    # # Remove non-patient folders 
    # folderToremove = []
    # for n in range(len(patientFolderList)):
    #     if patientFolderList[n][0:3]!= "QIN":
    #         folderToremove.append(patientFolderList[n])
    # for n in range(len(folderToremove)):
    #     patientFolderList.remove(folderToremove[n])
    
    
    # Loop through each patient 
    print(patientFolderList)
    for n in range(len(patientFolderList)): 
        print('Processing:',patientFolderList[n])
        with alive_bar(6) as bar:
    
            # Access the folder and Brain Mask file 
            patientFolderPath = os.path.join(rater_dir,patientFolderList[n])
            BETFolderPath = os.path.join(patientFolderPath,'BETimages')
            sleep(0.03) #update progress bar
            bar()
    
            # List of all files in nii folder 
            filesInBET = [f for f in os.listdir(BETFolderPath) if os.path.isfile(os.path.join(BETFolderPath,f))]
    
            # Create a new folder in the patient folder 
            normalised_path = os.path.join(patientFolderPath,'BET_normalised')
            if not os.path.isdir(normalised_path):
                os.makedirs(normalised_path)
    
            # Extract maps
            rBV_map_path = os.path.join(BETFolderPath,'BET_rBV_LC.nii')
            rBV_noLC_map_path = os.path.join(BETFolderPath,'BET_rBV_noLC.nii')
            rBF_map_path = os.path.join(BETFolderPath,'BET_rBF_LC.nii')
            MTT_map_path = os.path.join(BETFolderPath,'BET_MTT_LC.nii')
            sleep(0.03)
            bar()
            
            # Extract GM and WM masks and read as an image 
            gm_wm_path = os.path.join(BETFolderPath,'MaskGMWM.nii')
            gm_wm_mask = sitk.ReadImage(gm_wm_path)     
            rBV_map = sitk.ReadImage(rBV_map_path)
    
            # Transform the mask 
            gm_wm_mask.SetSpacing(rBV_map.GetSpacing())
            gm_wm_mask.SetOrigin(rBV_map.GetOrigin())
            gm_wm_mask.SetDirection(rBV_map.GetDirection()) 
            sleep(0.03)
            bar()
    
            # calculate the normalised maps for each and export the mape to the new folder 
            rBV_norm = normalise(rBV_map_path,gm_wm_mask)
            sitk.WriteImage(rBV_norm,os.path.join(normalised_path,'normalised_rBV_LC.nii'))
            sleep(0.03)
            bar()
    
            rBV_noLC_norm = normalise(rBV_noLC_map_path,gm_wm_mask)
            sitk.WriteImage(rBV_norm,os.path.join(normalised_path,'normalised_rBV_noLC.nii'))
            sleep(0.03)
            bar()
    
            rBF_norm = normalise(rBF_map_path,gm_wm_mask)
            sitk.WriteImage(rBV_norm,os.path.join(normalised_path,'normalised_rBF_LC.nii'))
            sleep(0.03)
            bar()

            MTT_norm = normalise(MTT_map_path,gm_wm_mask)
            sitk.WriteImage(MTT_norm,os.path.join(normalised_path,'normalised_MTT_LC.nii'))
            sleep(0.03)
            bar()
    
