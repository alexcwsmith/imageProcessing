#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 14:14:35 2020

@author: smith
"""

import os
from skimage import io
from skimage import morphology
import numpy as np
import time

directory='/d2/studies/ClearMap/FosTRAP_ChR2_March2021/'
sample='F3_RT'
crop = {'x':110,'y':140,'width':1900,'height':2400}

def skel(directory, sample, crop=None, flip=None):
    """Threshold and skeletonize 3D axons as described in Friedmann et al., PNAS 2020 (PMID: 32358193).
    
    Parameters
    ----------
    directory : str
        Path to base directory containing all sample subdirectories.
    sample : str
        Name of sample (must be name of subdirectory in directory).
    crop : dict (optional, default None)
        Dictionary of {x:int, y:int, width:int, height:int} to give cropping coordinates (ImageJ selection format).
    flip : str (optional, default None)
        Axis to flip/invert after completing. Valid arguments are 'x', 'y', 'z'.
        
    Returns
    -------
    None.
    """
    print("Started " + sample + ' ' + time.ctime())
    path = os.path.join(directory, sample + '/seg-ChR2/')
    ims = io.ImageCollection(path+'*.ome.tif', load_func=io.imread)
    data = ims.concatenate()
    if crop:
        rawshape=data.shape
        x=crop['x']
        y=crop['y']
        width=crop['width']
        height=crop['height']
        data = data[:,y:y+height,x:x+width]
        print("Cropped data from " + str(rawshape) + " to " + str(data.shape) + " at " + time.ctime())
    cat = np.zeros(shape=(data.shape), dtype='float32')
    for i in range(2,10,1):
        i=i/10
        im = (data>i).astype('float32')
        skel = morphology.skeletonize_3d(im).astype('float32')*i
        print(str(i) + " completed at " + time.ctime())
        cat = cat+skel
    if flip=='y':
        cat = np.flip(cat, axis=1)
    elif flip=='x':
        cat = np.flip(cat, axis=2)
    elif flip=='z':
        cat = np.flip(cat, axis=0)
    io.imsave(os.path.join(directory, sample + '/' + sample + '_ThresholdedSkeleton3D.tif'), cat, check_contrast=False)
    print("Finished " + sample + ' ' + time.ctime())

skel(directory, sample, crop, flip='y')


