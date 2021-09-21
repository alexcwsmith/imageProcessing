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


def Skeletonize3D(directory, start=2, save='planes', crop=False, flip='y', dtype=False, debug=False):
    """Skeletonize TrailMap results.
    
    Parameters
    ----------
    directory : string
        Path to directory with segmented data.
    start : int
        Integer for starting threshold for scaling/concatenation.
    save : str (optional, default 'stack')
        Specify 'stack' or 'planes' to save as single image stack or individual files with each Z-plane.
    crop : dict (optional, default False)
        Dictionary with ImageJ-format cropping coordinates ({'width':, 'height':, 'X':, 'Y':,})
    flip : string (optional, default 'y')
        Option to flip axis, can be any combination of 'xyz'.
    dtype : str (optional, default 'uint16')
        Specify 'uint16' to automatically convert to 16-bit. None or False remains 32 bit.
    debug : bool (optional, default False)
        If set to True, only process 50 planes in center of stack for quicker debugging.
    """
    #Load Data:
    if not directory.endswith('/'):
        directory = directory+'/'
    sample = directory.split('/')[-3]
    print("Started " + time.ctime())
    ims = io.ImageCollection(os.path.join(directory, '*.tif'), load_func=io.imread)
    data = ims.concatenate()
    #If in DEBUG mode, select only 50 plans in center of stack for quick processing.
    if debug:
        center = data.shape[0]//2
        data = data[center-25:center+25,:,:]
    #Optionally crop:
    if crop:
        rawshape=data.shape
        data = data[:,crop['y']:crop['y']+crop['height'],crop['x']:crop['x']+crop['width']]
        print("Cropped data from " + str(rawshape) + " to " + str(data.shape) + " at " + time.ctime())
    elif not crop:
        print("Using full data of shape " + str(data.shape))
    cat = np.zeros(shape=(data.shape), dtype='float32') #Create output array
    #Loop through thresholds 0.2 -> 0.9, extract signal, scale, and combine
    for i in range(start,10,1):
        print(str(i) + " started at " + time.ctime())
        i=i/10
        im = (data>i).astype('float32')
        print("Data thresholded for " + str(i) + " at " + time.ctime())
        skel = morphology.skeletonize_3d(im).astype('float32')*i
        print(str(i) + " completed at " + time.ctime())
        cat = cat+skel
    #Optionally flip along the x, y, or z axis:
    if flip:
        if 'y' in flip:
            cat = np.flip(cat, axis=1)
        if 'x' in flip:
            cat = np.flip(cat, axis=2)
        if 'z' in flip:
            cat = np.flip(cat, axis=0)
    if dtype:
        cat = cat.astype(dtype) #have not tested that this results in same pixel values as changing image type in ImageJ.
    #Save the result image stack:
    try:
        if save=='stack' and not debug:
            io.imsave(os.path.join(directory, '../' + sample + '_ThresholdedSkeleton3D_Start'+str(start)+'.tif'), cat, check_contrast=False, bigtiff=True)
        elif save=='stack' and debug:
            io.imsave(os.path.join(directory, '../' + sample + '_ThresholdedSkeleton3D_Start'+str(start)+'_DEBUG.tif'), cat, check_contrast=False, bigtiff=True)
        elif save=='planes':
            planeDir = os.path.join(directory, '../'+sample+'_ThresholdedSkeleton_Planes')
            if not os.path.exists(planeDir):
                os.mkdir(planeDir)
            for plane in list(range(cat.shape[0])):
                planeNum = str(plane).zfill(4)
                if dtype=='uint16':
                    if not debug:
                        io.imsave(os.path.join(planeDir, sample+'_ThresholdedSkeleton3D_start'+str(int(start*10))+'_Z'+str(planeNum)+'_16bit.tif'), img_as_uint(cat[plane]/cat[plane].max()), check_contrast=False)
                    elif debug:
                        io.imsave(os.path.join(planeDir, sample+'_ThresholdedSkeleton3D_start'+str(int(start*10))+'_Z'+str(planeNum)+'_16bit_DEBUG.tif'), img_as_uint(cat[plane]/cat[plane].max()), check_contrast=False)
                else:
                    if not debug:
                        io.imsave(os.path.join(planeDir, sample+'_ThresholdedSkeleton3D_start'+str(int(start*10))+'Z'+str(planeNum)+'.tif'), cat[plane], check_contrast=False)                        
                    elif debug:
                        io.imsave(os.path.join(planeDir, sample+'_ThresholdedSkeleton3D_start'+str(int(start*10))+'Z'+str(planeNum)+'_DEBUG.tif'), cat[plane], check_contrast=False)                        
    except PermissionError:
        homedir = os.path.expanduser('~/')
        print("You do not have write permissions for " + str(directory) + '\n' + "Saving stack to " + str(homedir) + " instead.")
        io.imsave(os.path.join(homedir, sample + '_ThresholdedSkeleton3D_Start'+str(start)+'.tif'), cat, check_contrast=False, bigtiff=True)
    print("Finished " + sample + ' ' + time.ctime())
    return cat

def clipBelowThreshold(directory, threshold=0.6, save='planes', crop=False, flip=False, dtype='uint16', debug=False):
    """Clip all pixel values below segmentation confidence threshold.
    
    Parameters
    ----------
    directory : string
        Path to Trailmap segmented image files.
    threshold : float
        Float (0 to 1) of values to clip pixels below.
    save : string (optional, default 'planes')
        'planes' or 'stack', sets mode for saving.
    crop : dict (optional, default False)
        Dictionary with ImageJ-format cropping coordinates ({'width':, 'height':, 'X':, 'Y':,})
    flip : str (optional, default False)
        Any combination of 'xyz', to flip image about specified axis.
    dtype : str (optional, default 'uint16')
        Specify 'uint16' to automatically convert to 16-bit. None or False remains 32 bit.
    debug : bool (optional, default False)
        If set to True, only process 50 planes in center of stack for quicker debugging.
    """
    if not directory.endswith('/'):
        directory = directory+'/'
    sample = directory.split('/')[-3]
    print("Started " + sample + " at " + time.ctime())
    ims = io.ImageCollection(os.path.join(directory, '*.tif'), load_func=io.imread)
    data = ims.concatenate()
    #If in DEBUG mode, select only 50 plans in center of stack for quick processing.
    if debug:
        print("Running in DEBUG mode. Processing 50 planes from center of stack.")
        center = data.shape[0]//2
        data = data[center-25:center+25,:,:]
    #Optionally crop:
    if crop:
        rawshape=data.shape
        data = data[:,crop['y']:crop['y']+crop['height'],crop['x']:crop['x']+crop['width']]
        print("Cropped data from " + str(rawshape) + " to " + str(data.shape) + " at " + time.ctime())
    elif not crop:
        print("Using full data of shape " + str(data.shape))
    print("Clipping values below " + str(threshold) + " started at " + time.ctime())
    im = (data>threshold).astype('float32')
    skel = morphology.skeletonize_3d(im).astype('float32')
    print("Values below " + str(threshold) + " clipped at " + time.ctime())
    try:
        if save=='stack':
            if debug:
                io.imsave(os.path.join(directory, '../' + sample + '_ThresholdedSkeleton3D_Thresh'+str(int(threshold*10))+'_DEBUG.tif'), skel, check_contrast=False, bigtiff=True)
            elif not debug:
                io.imsave(os.path.join(directory, '../' + sample + '_ThresholdedSkeleton3D_Thresh'+str(int(threshold*10))+'.tif'), skel, check_contrast=False, bigtiff=True)
        elif save=='planes':
            planeDir = os.path.join(directory, '../'+sample+'_ThresholdedSkeleton'+str(int(threshold*10))+'_Planes', check_contrast=False)
            if not os.path.exists(planeDir):
                os.mkdir(planeDir)
            for plane in list(range(skel.shape[0])):
                planeNum = str(plane).zfill(4)
                if dtype=='uint16':
                    if debug:
                        io.imsave(os.path.join(planeDir, sample+'_ThresholdedSkeleton3D_Z'+str(planeNum)+'_DEBUG.tif'), img_as_uint(skel[plane]/(skel[plane].max())), check_contrast=False)
                    elif not debug:
                        io.imsave(os.path.join(planeDir, sample+'_ThresholdedSkeleton3D_Z'+str(planeNum)+'.tif'), img_as_uint(skel[plane]/(skel[plane].max())), check_contrast=False)
                else:
                    if debug:
                        io.imsave(os.path.join(planeDir, sample+'_ThresholdedSkeleton3D_Z'+str(planeNum)+'_DEBUG.tif'), skel[plane], check_contrast=False)                        
                    elif not debug:
                        io.imsave(os.path.join(planeDir, sample+'_ThresholdedSkeleton3D_Z'+str(planeNum)+'.tif'), skel[plane], check_contrast=False)                        
    except PermissionError:
        homedir = os.path.expanduser('~/')
        print("You do not have write permissions for " + str(directory) + '\n' + "Saving to " + str(homedir) + " instead.")
        io.imsave(os.path.join(homedir, sample + '_ThresholdedSkeleton3D_Thresh'+str(int(threshold*10))+'.tif'), skel, check_contrast=False, bigtiff=True)
    print("Finished " + sample + ' ' + time.ctime())
    return skel


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--skel',type=bool,default=False,help='Whether to run skeletonization function')
    p.add_argument('--clip',type=bool,default=False,help='Whether to clip pixels values below a threshold')
    p.add_argument('--directory',type=str,default=os.getcwd(),help='Directory with segmented images.')
    p.add_argument('--save',type=str,default='planes',help="Whether to save a 'stack' or 'planes'.")
    p.add_argument('--start',type=int,default=2,help='Starting threshold for skeletonization. Integer between 2 and 10')
    p.add_argument('--threshold',type=float,default=0.7,help='Minimum pixel value for clipping function, pixels under this threshold will be clipped.')
    p.add_argument('--dtype',type=str,default='uint16',help="'uint16' converts data automatically. False or None leaves data as float32")
    p.add_argument('--flip',type=str,default=None,help='xyz axes to flip')
    p.add_argument('--crop',type=bool,default=False,help='Whether to crop the data prior to processing.')
    p.add_argument('--width',type=int,default=None,help="Width to crop to. Only used if crop=True")
    p.add_argument('--height',type=int,default=None,help="Height to crop to. Only used if crop=True")
    p.add_argument('--X',type=int,default=None,help='Starting X coordinate for cropping. Only used if crop=True')
    p.add_argument('--Y',type=int,default=None,help='Starting Y coordinate for cropping. Only used if crop=True')
    p.add_argument('--debug',type=bool,default=False,help="If true, only process 50 planes from center of stack for quick debugging")
    args = p.parse_args()
    if not args.skel and not args.clip:
        raise KeyError("Either skel or clip argument must be set to True.")
    if args.skel and args.clip:
        raise KeyError("Only one of skel or clip must be set to True.")#TODO allow performing both functions in-line.
    if args.crop:
        cropCoords = {}
        cropCoords = {'width':args.width,'height':args.height,'x':args.X,'y':args.Y}
        args.crop=cropCoords
    elif not args.crop:
        args.crop = None
    if args.skel:
        Skeletonize3D(args.directory, start=args.start, save=args.save, crop=args.crop, flip=args.flip, dtype=args.dtype, debug=args.debug)
    elif args.clip:
        clipBelowThreshold(args.directory, threshold=args.threshold, save=args.save, crop=args.crop, flip=args.flip, dtype=args.dtype)
    
