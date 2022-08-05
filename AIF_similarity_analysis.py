# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 11:58:00 2022

@author: Caterina Brighi
"""

#%% Import functions 


import numpy as np
import pandas as pd
import os
from skimage.metrics import mean_squared_error as mse
from skimage.metrics import normalized_root_mse as nrmse
from matplotlib import pyplot as plt



def rmse(vector1, vector2):
    '''This function returns the root mean squared error between two vectors.
    --------------------------
    vector1, vector2 : ndarray'''
    return np.sqrt(mse(vector1, vector2))    


#%% Set Working directory
        
# Path of the main folder 
DataSupraDir = 'Path to users supradirectory'


reters_path = [ f.path for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of the paths to the raters directories
raters_name = [ f.name for f in os.scandir(DataSupraDir) if f.is_dir() ] #Create a list of raters names
raters_name.remove('Analysis results')

subjs_name = [ f.name for f in os.scandir(DataSupraDir+raters_name[0]) if f.is_dir() ] 

for i in subjs_name:
    os.mkdir(DataSupraDir+'Analysis results/AIF analysis/'+i)



IntraWriter = pd.ExcelWriter(DataSupraDir +'Analysis results/AIF analysis/AIF_Similarity_Results.xlsx', engine='xlsxwriter')


for i in raters_name:
    raters_name = [ f.name for f in os.scandir(DataSupraDir) if f.is_dir() ]
    raters_name.remove('Analysis results')
    raters_name.remove(i)

    n_raters = len(raters_name) #Total number of raters
   
    # Extract number of patient folders as a list 
    subjs_name = [ f.name for f in os.scandir(DataSupraDir+raters_name[0]) if f.is_dir() ] #Create a list of raters names
    n_subjs = len(subjs_name) #Total number of subjects
    
    df = pd.DataFrame(columns=['Subject_ID', 'RMSE', 'NRMSE'])
    
    #%%Create a for loop to perform image analysis on each rater sequentially
    
    for subject in subjs_name:
        
        print('Analysing patient '+subject)
        
        AIF_dic = {'Subject_ID':subject}
          
    
        # Loop through each rater 
        for idx, rater in enumerate(raters_name):
            print(idx, rater)
            rater_dir = DataSupraDir+rater          
                
            # Set path to patient subfolders
            patientFolderPath = os.path.join(rater_dir,subject)
            
            # Set path to AIF.xlsx file
            AIFpath = os.path.join(patientFolderPath,'AIF.xlsx')
            
            # Read AIF data in a dataframe
            AIF_df = pd.read_excel(AIFpath, sheet_name='Sheet1')
            AIF = np.array(AIF_df.loc[:, ' Average AIF'])
            Time = np.array(AIF_df.loc[:, 'Time (s)'])
                
            #Append AIF data for this rater to dictionary
            AIF_dic['AIF for rater'+str(idx+1)]=AIF
            AIF_dic['Time for rater'+str(idx+1)]=Time
            
            if rater == 'EdwardChen':
                AIF_dic['Rater'+str(idx+1)+'name']='Expert'
            elif rater == 'SameerDave':
                AIF_dic['Rater'+str(idx+1)+'name']='Trainee'
            elif rater == 'AthenaKam':
                AIF_dic['Rater'+str(idx+1)+'name']='Non-expert'
            

#%% Calculate similarity metrics within ROI
        
        print('Calculate similarity metrics between AIF for subject '+subject)
            
        RMSE = rmse(AIF_dic.get('AIF for rater1'), AIF_dic.get('AIF for rater2'))
        NRMSE = nrmse(AIF_dic.get('AIF for rater1'), AIF_dic.get('AIF for rater2'))
        
        #Append these metrics to AIF dictionary
        AIF_dic['RMSE']=RMSE
        AIF_dic['NRMSE']=NRMSE
        
        #Plot AIF for two raters
        fig, ax = plt.subplots(figsize=(10, 6))
        plt.plot(AIF_dic.get('Time for rater1'), AIF_dic.get('AIF for rater1'), label='AIF '+ AIF_dic.get('Rater1name'))
        plt.plot(AIF_dic.get('Time for rater1'), AIF_dic.get('AIF for rater2'), label='AIF '+ AIF_dic.get('Rater2name'))
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('AIF')
        fig.tight_layout(pad=2.0)
        plt.savefig(DataSupraDir+'Analysis results/AIF analysis/'+subject+'/AIF '+AIF_dic.get('Rater1name')+'-'+AIF_dic.get('Rater2name')+'.png')     
        plt.close()
        
        del AIF_dic['AIF for rater1']
        del AIF_dic['AIF for rater2']
        del AIF_dic['Time for rater1']
        del AIF_dic['Time for rater2']
        del AIF_dic['Rater1name']
        del AIF_dic['Rater2name']
        
    #%% Export subjects stats into excel spreadsheet
    
        print('Export subjects '+subject+' metrics into excel spreadsheet')
            
        df = df.append(AIF_dic, ignore_index=True)
        
    df.to_excel(IntraWriter, sheet_name=raters_name[0]+'-'+raters_name[1], index=False)
        
IntraWriter.save()