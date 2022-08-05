# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 12:54:41 2022

@author: Caterina Brighi and Athena Kam
"""

# This script will loop through the RDS to create masked maps for everything in the patient folders 
import SimpleITK as sitk
import os 
import numpy as np

# Progress bar 
from alive_progress import alive_bar
from time import sleep

# Define masking function 


def generate_mask(image, roi):
    
    '''Returns the masked image of the roi applied to the image.
    Remember to save the masked image file after applying this function.'''
    
    masking = sitk.MaskImageFilter() 
    mask = masking.Execute(image, roi)
    return mask


# Path of the main folder 
DataSupraDir = 'Path to users supradirectory'

reters_path = [ f.path for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of the paths to the raters directories
raters_name = [ f.name for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of raters names
# raters_name.remove("AthenaKam")
# raters_name.remove("EdwardChen")

n_raters = len(raters_name) #Total number of raters

#%%Create a for loop to perform image analysis on each rater sequentially

for rater in raters_name:
    rater_dir = DataSupraDir+rater

    # Extract number of patient folders as a list 
    patientFolderList = [f for f in os.listdir(rater_dir) if os.path.isdir(os.path.join(rater_dir,f))]
    
    # Define maps
    mapList = ['DelayMap.nii','LeakageMap.nii','MTT_LC.nii','rBF_LC.nii','rBV_LC.nii','rBV_noLC.nii']
    
    # Loop through each patient 
    print(patientFolderList)
    for n in range(len(patientFolderList)): 
        print('Processing:',patientFolderList[n])
        with alive_bar(6) as bar:
    
            # Set path to the masks
            patientFolderPath = os.path.join(rater_dir,patientFolderList[n])
            
            brainMaskPath = os.path.join(patientFolderPath,'MaskBrain.nii')
            brainMask = sitk.ReadImage(brainMaskPath, sitk.sitkUInt16)
    
            # Create a new folder in the patient folder 
            masked_path = os.path.join(patientFolderPath,'BETimages')
            if not os.path.isdir(masked_path):
                os.makedirs(masked_path)
             
            # Loop through the images that need to be masked in the patient folder 
            for m in range(len(mapList)):
                mapFile = os.path.join(patientFolderPath,mapList[m])
                image = sitk.ReadImage(mapFile)
                image.SetSpacing(brainMask.GetSpacing())
                image.SetOrigin(brainMask.GetOrigin())
                image.SetDirection(brainMask.GetDirection())
                map_masked = generate_mask(image, brainMask)
    
                # Save as a new map 
                filename = 'BET_'+ mapList[m]
                outPath = os.path.join(masked_path,filename)
                sitk.WriteImage(map_masked,outPath)
            
            #BET also the GMWM mask
            GMWMMaskPath = os.path.join(patientFolderPath,'MaskGMWM.nii')
            GMWMMask = sitk.ReadImage(GMWMMaskPath)
            GMWMMask.SetSpacing(brainMask.GetSpacing())
            GMWMMask.SetOrigin(brainMask.GetOrigin())
            GMWMMask.SetDirection(brainMask.GetDirection())
            GMWM_masked = generate_mask(GMWMMask, brainMask)
            sitk.WriteImage(GMWM_masked,os.path.join(masked_path,'MaskGMWM.nii'))
            
            #Save also other images in BET folder
            tumornMaskPath = os.path.join(patientFolderPath,'MaskTumor.nii')
            tumorMask = sitk.ReadImage(tumornMaskPath)
            tumorMask.SetSpacing(brainMask.GetSpacing())
            tumorMask.SetOrigin(brainMask.GetOrigin())
            tumorMask.SetDirection(brainMask.GetDirection())
            sitk.WriteImage(tumorMask,os.path.join(masked_path,'MaskTumor.nii'))
            
            sitk.WriteImage(brainMask,os.path.join(masked_path,'MaskBrain.nii'))
            
            AnatomicPath = os.path.join(patientFolderPath,'Anatomic.nii')
            Anatomic = sitk.ReadImage(AnatomicPath)
            Anatomic.SetSpacing(brainMask.GetSpacing())
            Anatomic.SetOrigin(brainMask.GetOrigin())
            Anatomic.SetDirection(brainMask.GetDirection())
            sitk.WriteImage(Anatomic,os.path.join(masked_path,'Anatomic.nii'))
            
          
            sleep(0.03)
            bar() 
    
