#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 14:14:35 2020

@author: smith
"""

import os
from skimage import io
from skimage import morphology
from skimage import img_as_uint
import numpy as np
import time
import argparse

#If you are planning to run this straight from the command line, delete or comment out the following 
#4 lines (everything before the function definition). See bottom of script for CLI options, or just input
#'python3 TrailMap_Skeletonize3D --help' into the terminal
directory='/d2/studies/ClearMap/ROC_ChR2/ROC31/seg-ChR2/'
if not os.path.exists(directory):
    raise FileNotFoundError("Directory does not exist, check path")
crop=None #{'width':1900,'height':2300,'x':135,'y':260} #ROC31

def skel(directory, crop=None, flip='y', convertTo16Bit=False, debug=False):
    """Skeletonize TrailMap results as described in Friedmann et al., PNAS 2020 (PMID: 32358193)
    
    Copied from methods section:
        "Probabilistic volumes from the output of TrailMap were binarized at 8 separate thresholds, from 0.2 to 
        0.9. Each of these eight volumes was skeletonized in three dimensions (skimage-morphology-skeletonize_3d).
        The logical binary skeleton outputs were each weighted by the initial probability threshold used to generate
        them and subsequently were summed together"
    
    The above is done in lines 61-67 of this script (starting with the line 'for i in range(2,10,1):')
    
    Parameters
    ----------
    directory : string
        Path to directory with segmented data.
    crop : dict (optional, default None)
        Dictionary with ImageJ-format cropping coordinates: {width:int, height:int, x:int, y:int}
    flip : string (optional, default 'y')
        Option to flip (reverse) axis, can be any combination of 'xyz'.
    convertTo16Bit : bool
        Whether to convert the image to 16 bit before saving. Have not double checked that this gives the same results as converting in ImageJ yet though.
    debug : bool
        Run in debug mode, process only a small substack of data, 50 planes in the center of the stack.
    """
    sample = os.path.realpath(directory).split('/')[-2] #the number at the end of this line may change based on the file path of directory, for the current path -2 is correct (sample variable is ROC_9)"
    print("Started " + str(sample) + " at " + time.ctime())
    ims = io.ImageCollection(load_pattern=os.path.join(directory, '*.tif'), load_func=io.imread)
    data = ims.concatenate()
    if crop:
        rawshape=data.shape
        data = data[:,crop['y']:crop['y']+crop['height'],crop['x']:crop['x']+crop['width']]
        print("Cropped data from " + str(rawshape) + " to " + str(data.shape) + " at " + time.ctime())
    if debug:
        center = data.shape[0]//2
        data = data[center-25:center+25,:,:]
    cat = np.zeros(shape=(data.shape), dtype='float32')
    for i in range(2,10,1):
        print(str(i) + " started at " + time.ctime())
        i=i/10
        im = (data>i).astype('float32')
        skel = morphology.skeletonize_3d(im).astype('float32')*i
        print(str(i) + " completed at " + time.ctime())
        cat = cat+skel
    if flip:
        if 'y' in flip:
            cat = np.flip(cat, axis=1)
        if 'x' in flip:
            cat = np.flip(cat, axis=2)
        if 'z' in flip:
            cat = np.flip(cat, axis=0)
    if not convertTo16Bit and not debug:
        io.imsave(os.path.join(directory, sample + '_ThresholdedSkeleton3D.tif'), cat, check_contrast=False)
    elif not convertTo16Bit and debug:
        io.imsave(os.path.join(directory, sample + '_ThresholdedSkeleton3D_DEBUG.tif'), cat, check_contrast=False)
    elif convertTo16Bit and not debug:
        cat = cat/cat.max()
        cat = img_as_uint(cat)
        io.imsave(os.path.join(directory, sample+ '_ThresholdedSkeleton3D_16bit.tif'), cat, check_contrast=False)
    elif convertTo16Bit and debug:
        cat = cat/cat.max()
        cat = img_as_uint(cat)
        io.imsave(os.path.join(directory, sample + '_ThresholdedSkeleton3D_16bit_DEBUG.tif'), cat, check_contrast=False)
    print("Finished " + sample + ' ' + time.ctime())
    return cat

cat = skel(directory, crop=crop, flip='y', convertTo16Bit=True, debug=False)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--directory',type=str,default=os.getcwd(),help='Directory with segmented images.')
    p.add_argument('--crop',type=dict,default=None,help='Dict with ImageJ-format cropping coordinates.')
    p.add_argument('--flip',type=str,default=None,help='xyz axes to flip')
    p.add_argument('--convertTo16Bit',type=bool,default=False,help='Whether to convert to 16 bit image before saving')
    p.add_argument('--debug',type=bool,default=False,help='If True, only process 50 planes from center of the stack for faster debugging')
    args = p.parse_args()
    skel(args.directory, crop=args.crop, flip=args.flip, convertTo16Bit=args.convertTo16Bit, debug=args.debug)
    
    
