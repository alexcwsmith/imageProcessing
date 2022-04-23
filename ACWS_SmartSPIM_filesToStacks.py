#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 12:58:26 2022

@author: smith
"""

from skimage import io
from aicsimageio.writers import OmeTiffWriter
from aicsimageio.types import PhysicalPixelSizes as ps
import numpy as np
import os
import argparse
import time

def get_dir(directory, ext=None):
    if ext:
        return sorted(os.path.join(directory, x) for x in os.listdir(directory) if x.endswith(ext))
    elif not ext:
        return sorted(os.path.join(directory, x) for x in os.listdir(directory))

def createStack(directory, ext='.tif', crop=None, overwrite=False,
                pixelSizes=(4,4,4), channelName=None, flip=None):
    """Create ome.tif stack from stitched image series produced by SmartSPIM.
    
    Parameters
    ----------
    directory : str
        Path containing image series.
    pixelSizes : tuple or int/float
        Physical pixel sizes (ints or floats). If single number is given, it is assumed that is pixel size in all 3 dimensions.
    ext : str
        File extension of image series
    channelName : str
        Name of channel in folder. If None laser wavelength is used.
    """
    #file/directory organization
    print("Started " + str(time.ctime()))
    if directory.endswith('/'):
        dirname = directory.split('/')[-2]
        sampleName = directory.split('/')[-3].split('_')[0]
    else:
        dirname = directory.split('/')[-1]
        sampleName = directory.split('/')[-2].split('_')[0]
    baseDir = directory.split('Ex_')[0]
    laser = dirname.split('_')[1]
    if not channelName and int(laser)==488:
        channelName='auto'
    elif not channelName and int(laser)!=488:
        channelName=str(laser)
    saveDir=os.path.join(baseDir, 'stacks')
    if not os.path.exists(saveDir):
        os.mkdir(saveDir)

    if not overwrite:
        if os.path.exists(os.path.join(saveDir, sampleName+'_Stitched_'+channelName+'.ome.tif')):
            overwrite = input("File already exists, do you want to overwrite (y/n)?")
            if overwrite.lower().startswith('n'):
                return

    if not isinstance(pixelSizes, tuple):
        pixelSizes=(pixelSizes, pixelSizes, pixelSizes)
    pix = ps(Z=pixelSizes[0], Y=pixelSizes[1], X=pixelSizes[2])
    if crop and not isinstance(crop, dict):
        raise ValueError("Crop argument must be dictionary of width, height, x, y. (ImageJ format cropping coordinates)")
    
    #read images
    files = get_dir(directory, ext)
    print("Concatenating " + str(len(files)) + " images")
    imc = io.imread_collection(files)
    im = io.concatenate_images(imc)
    print("Data loaded at " + str(time.ctime()) + ", saving image stack")
    if crop:
        im = im[:,crop['y']:crop['y']+crop['height'], crop['x']:crop['x']+crop['width']]
    if flip:
        if 'y' in flip:
            im = np.flip(im, 1)
        if 'x' in flip:
            im = np.flip(im, 2)
        if 'z' in flip:
            im = np.flip(im, 0)
    
    #write image stack:
    OmeTiffWriter.save(data=im, 
                    uri=os.path.join(saveDir, sampleName+'_Stitched_'+channelName+'.ome.tif'), 
                    physical_pixel_sizes=pix, 
                    dim_order="ZYX")
    print("Finished at " + str(time.ctime()))


if __name__=='__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-d', '--directory', type=str, default=os.getcwd(), help='Path to directory containing stitched image files')
    p.add_argument('-p', '--pixelSizes', type=tuple, default=(4,4,4), help='ZYX pixel sizes')
    p.add_argument('-e', '--extension', type=str, default='.tif', help='Extension of image files')
    p.add_argument('-n', '--name', type=str, default=None, help='Name to use for saving channel stack, e.g. fos')
    p.add_argument('-f', '--flip', type=str, default=None, help='Combination of xyz axes to flip')
    p.add_argument('-o', '--overwrite', type=str, default=False, help='Force overwrite file if it already exists')
    args = p.parse_args()
    createStack(args.directory, pixelSizes=args.pixelSizes, ext=args.extension, channelName=args.name, flip=args.flip, overwrite=args.overwrite)