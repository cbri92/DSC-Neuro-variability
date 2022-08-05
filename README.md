# DSC Neuro variability

**Author:** *Caterina Brighi*

The code available on this repository was created for a project aiming at evaluating the inter-rater variability in selecting AIF and generating parametric maps of perfusion from DSC MRI data in brain cancer patients.

## Setup/Build/Install

To use this code you will need the following python packages installed: matplotlib, SimpleITK, numpy, pandas, datetime, os, glob, gzip, shutil, xlsxwriter, scipy, functools, statistics, seaborn, pingouin, alive-progress.

## Usage

Before running the scripts, make sure you edit the correct path as indicated at the top of each script. The scripts are to be run in the following order:

1. *BrainExtraction.py* - This script will perform brain extraction on all images for all the raters in the main data folder.

2. *ImageNormalisation.py* - This script will normalise all the brain extracted images wrt the mean intensity in a ROI placed in the contralateral healthy brain.

3. *ICCAnalysis_brain.py* and *ICCAnalysis_tumor.py* - These scripts will calculate the per-voxel ICC mean-rating (k= 3) from multiple raters, absolute-agreement, 2-way random effects model in the brain and in the tumour masks, respectively.

4. *Brain_Repeatability_Analysis.py* and *Tumor_Repeatability_Analysis.py* - These scripts will calculate the per-voxel repeatability metrics in the brain and in the tumour masks, respectively.

5. *Brain_Similarity_Analysis.py* and *Tumor_Similarity_Analysis.py* - These scripts will calculate similarity metrics (MSE, RMSE, NRMSE, PSNR, SSIM) between pairs of raters in the brain and in the tumour masks, respectively.

6. *Similarity_metrics_reformatting.py* - This script will reformat the results of the similarity metrics (MSE, RMSE, NRMSE, PSNR, SSIM) between pairs of raters in a format that facilitates statistical analysis.

7. *AIF_curve_parameters.py* - This script calculates AIF curve parameters (MP, FWHM, TTP, AUC, M) for each subject for each rater.

8. *AIF_similarity_analysis.py* - This script calculates similarity metrics (RMSE and NRMSE) for AIF between pairs of raters.

NOTE: The *ImageSimilarityMetricsFunctions.py* has to be in the same folder of the other scripts while running them, as it contains functions to calculate similarity metrics.


## Directory Structure

The folder structure for the data should be the following:

    Data Supradirectory
          Rater1
              Subject001
                   AIF.xlsx
                   MaskBrain.nii
                   MaskGMWM.nii
                   rBV_noLC.nii
                   rBV_LC.nii
                   rBF_LC.nii
                   MTT_LC.nii
              SubjectXXX
              ...
          RaterX
          ...
              
    
