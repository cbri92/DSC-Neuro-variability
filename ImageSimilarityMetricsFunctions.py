# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 16:18:10 2022

@author: Caterina Brighi
"""


import numpy as np
import SimpleITK as sitk
from skimage.metrics import mean_squared_error as mse
from skimage.metrics import normalized_root_mse as nrmse
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import peak_signal_noise_ratio as psnr

#%%Similarity metrics between two images within an roi


def mse_inROI(image0, image1, ROI):
    
    """
    Compute the mean-squared error between two images within a ROI.
    Parameters
    ----------
    image0, image1, ROI : ndarray
        Images.  Any dimensionality, must have same shape.
    Returns
    -------
    mse : float
        The mean-squared error (MSE) metric calculated within the specified ROI.
    """
  
    image0 = sitk.GetArrayFromImage(image0)
    roi = sitk.GetArrayFromImage(ROI)
    mask = np.array(roi,dtype='bool')
    mask = np.invert(mask)
    image0 = image0[~np.array(mask)]
    
    image1 = sitk.GetArrayFromImage(image1)
    image1 = image1[~np.array(mask)]
        
    return mse(image0, image1)


def rmse_inROI(image0, image1, ROI):
    """
    Compute the root mean-squared error between two images within a ROI.
    Parameters
    ----------
    image0, image1, ROI : ndarray
        Images.  Any dimensionality, must have same shape.
    Returns
    -------
    mse : float
        The root mean-squared error (RMSE) metric calculated within the specified ROI.
    """
    return np.sqrt(mse_inROI(image0, image1, ROI))



def nrmse_inROI(image_true, image_test, ROI):
    """
    Compute the normalized root mean-squared error (NRMSE) between two
    images in a ROI.
    Parameters
    ----------
    image_true : ndarray
        Ground-truth image, same shape as im_test.
    image_test : ndarray
        Test image.
    ROI : ndarray
        ROI
    normalization : {'euclidean', 'min-max', 'mean'}, optional
        Controls the normalization method to use in the denominator of the
        NRMSE.  There is no standard method of normalization across the
        literature [1]_.  The methods available here are as follows:
        - 'euclidean' : normalize by the averaged Euclidean norm of
          ``im_true``::
              NRMSE = RMSE * sqrt(N) / || im_true ||
          where || . || denotes the Frobenius norm and ``N = im_true.size``.
          This result is equivalent to::
              NRMSE = || im_true - im_test || / || im_true ||.
        - 'min-max'   : normalize by the intensity range of ``im_true``.
        - 'mean'      : normalize by the mean of ``im_true``
    Returns
    -------
    nrmse : float
        The NRMSE metric.
    """

    image_true = sitk.GetArrayFromImage(image_true)
    roi = sitk.GetArrayFromImage(ROI)
    mask = np.array(roi,dtype='bool')
    mask = np.invert(mask)
    image_true = image_true[~np.array(mask)]
    
    image_test = sitk.GetArrayFromImage(image_test)
    image_test = image_test[~np.array(mask)]
    
    return nrmse(image_true, image_test)
    

def psnr_inROI(image_true, image_test, ROI):
    """
    Compute the peak signal to noise ratio (PSNR) for an image in a ROI.
    Parameters
    ----------
    image_true : ndarray
        Ground-truth image, same shape as im_test.
    image_test : ndarray
        Test image.
    ROI : ndarray
        ROI
    data_range : int, optional
        The data range of the input image (distance between minimum and
        maximum possible values).  By default, this is estimated from the image
        data-type.
    Returns
    -------
    psnr : float
        The PSNR metric.
    """
    image_true = sitk.GetArrayFromImage(image_true)
    roi = sitk.GetArrayFromImage(ROI)
    mask = np.array(roi,dtype='bool')
    mask = np.invert(mask)
    image_true = image_true[~np.array(mask)]
    
    image_test = sitk.GetArrayFromImage(image_test)
    image_test = image_test[~np.array(mask)]
    true_min, true_max = np.min(image_true), np.max(image_true)
    drange = true_max-true_min
    return psnr(image_true, image_test, data_range=drange)


def ssim_inROI(image0, image1, ROI):
    """
    Compute the mean structural similarity index between two images in a ROI.
    Parameters
    ----------
    im1, im2, ROI : ndarray
        Images. Any dimensionality with same shape.
    win_size : int or None, optional
        The side-length of the sliding window used in comparison. Must be an
        odd value. If `gaussian_weights` is True, this is ignored and the
        window size will depend on `sigma`.
    gradient : bool, optional
        If True, also return the gradient with respect to im2.
    data_range : float, optional
        The data range of the input image (distance between minimum and
        maximum possible values). By default, this is estimated from the image
        data-type.
    channel_axis : int or None, optional
        If None, the image is assumed to be a grayscale (single channel) image.
        Otherwise, this parameter indicates which axis of the array corresponds
        to channels.
        .. versionadded:: 0.19
           ``channel_axis`` was added in 0.19.
    multichannel : bool, optional
        If True, treat the last dimension of the array as channels. Similarity
        calculations are done independently for each channel then averaged.
        This argument is deprecated: specify `channel_axis` instead.
    gaussian_weights : bool, optional
        If True, each patch has its mean and variance spatially weighted by a
        normalized Gaussian kernel of width sigma=1.5.
    full : bool, optional
        If True, also return the full structural similarity image.
    Other Parameters
    ----------------
    use_sample_covariance : bool
        If True, normalize covariances by N-1 rather than, N where N is the
        number of pixels within the sliding window.
    K1 : float
        Algorithm parameter, K1 (small constant, see [1]_).
    K2 : float
        Algorithm parameter, K2 (small constant, see [1]_).
    sigma : float
        Standard deviation for the Gaussian when `gaussian_weights` is True.
    Returns
    -------
    mssim : float
        The mean structural similarity index over the image.
    grad : ndarray
        The gradient of the structural similarity between im1 and im2 [2]_.
        This is only returned if `gradient` is set to True.
    S : ndarray
        The full SSIM image.  This is only returned if `full` is set to True.
    Notes
    -----
    To match the implementation of Wang et. al. [1]_, set `gaussian_weights`
    to True, `sigma` to 1.5, and `use_sample_covariance` to False.
    .. versionchanged:: 0.16
        This function was renamed from ``skimage.measure.compare_ssim`` to
        ``skimage.metrics.structural_similarity``.
    References
    ----------
    .. [1] Wang, Z., Bovik, A. C., Sheikh, H. R., & Simoncelli, E. P.
       (2004). Image quality assessment: From error visibility to
       structural similarity. IEEE Transactions on Image Processing,
       13, 600-612.
       https://ece.uwaterloo.ca/~z70wang/publications/ssim.pdf,
       :DOI:`10.1109/TIP.2003.819861`
    .. [2] Avanaki, A. N. (2009). Exact global histogram specification
       optimized for structural similarity. Optical Review, 16, 613-621.
       :arxiv:`0901.0065`
       :DOI:`10.1007/s10043-009-0119-z`
    """
    image0 = sitk.GetArrayFromImage(image0)
    roi = sitk.GetArrayFromImage(ROI)
    mask = np.array(roi,dtype='bool')
    mask = np.invert(mask)
    image0 = image0[~np.array(mask)]
    
    image1 = sitk.GetArrayFromImage(image1)
    image1 = image1[~np.array(mask)]
    
    return ssim(image0, image1)


def similarity_metrics_inROI(image1, image2, ROI):
    '''This function calculates the following similarity metrics between two images in a defined ROI:
        mean square error: mse
        root mean square error: rmse
        normalised root mean square error: nrmse
        peak signal to noise ratio: psnr
        structural similarity index: ssim
        -----------------------
        im1, im2, ROI : ndarray
        Images. Any dimensionality with same shape.'''
    
    mse = mse_inROI(image1, image2, ROI)
    rmse = rmse_inROI(image1, image2, ROI)
    nrmse = nrmse_inROI(image1, image2, ROI)
    psnr = psnr_inROI(image1, image2, ROI)
    ssim = ssim_inROI(image1, image2, ROI)
    
    
    return {'MSE':mse, 'RMSE':rmse, 'NRMSE':nrmse, 'PSNR':psnr, 'SSIM':ssim}
