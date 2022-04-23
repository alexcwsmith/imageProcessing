#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 09:03:58 2022

@author: smith
"""


from skimage import io
import os
import numpy as np
import glob
import argparse

def cropZeissTiles(directory, level=3, axis='X', cropDim=256, side='left', plugin='skimage', pixelSizes=None):
    """Crop tiles from the Zeiss Z1 light-sheet.
    
    Parameters
    ----------
    directory : str
        Path to directory containing TIF stacks for each tile.
    level : int (optional, default 3)
        Index of row or column to crop.
    axis : str (optional, default 'X')
        Whether level refers to X or Y axis.
    cropDim : int (optional, default 256)
        How many pixels to crop off.
    side : str (optional, default 'left')
        Which side to crop pixels off of, can be 'left', 'right', top', bottom'
    plugin : str (optional, default 'skimage')
        image IO plugin to use, either 'skimage' or 'aicsimageio'. aicsimageio supports writing .ome.tif with pixel size metadata
    pixelSizes : tuple (floats)
        Tuple of ZYX pixel sizes, if using aicsimageio plugin
    """
    saveDir = os.path.join(directory, 'cropped/')
    if not os.path.exists(saveDir):
        os.mkdir(saveDir)
    if axis.lower()=='x':
        files = glob.glob(os.path.join(directory, '*'+'0'+str(level)+'.ome.tif'))
    elif axis.lower()=='y':
        files = glob.glob(os.path.join(directory, '*'+'0'+str(level)+'x*.ome.tif'))
    print(str(len(files)) + ' files found for cropping')
    for f in files:
        n, e = os.path.splitext(os.path.basename(f))
        if n.endswith('.ome'):
            n=n.strip('.ome')
            e = '.ome'+e
        if plugin=='skimage':
            im = io.imread(f)
        elif plugin=='aicsimageio':
            from aicsimageio.readers import OmeTiffReader
            from aicsimageio.writers import OmeTiffWriter
            im = OmeTiffReader(f)
            if not pixelSizes:
                pixelSizes = im.physical_pixel_sizes
            elif pixelSizes:
                from aicsimageio.types import PhysicalPixelSizes as ps
                pixelSizes = ps(Z=pixelSizes[0], Y=pixelSizes[1], X=pixelSizes[2])
            im = im.get_image_data('ZYX')
        if side.lower()=='left':
            im = im[:,:,cropDim:]
        elif side.lower()=='right':
            im = im[:,:,:im.shape[2]-cropDim]
        elif side.lower()=='top':
            im = im[:,cropDim:,:]
        elif side.lower()=='bottom':
            im = im[:,:im.shape[1]-cropDim,:]
        else:
            raise ValueError("side argument must be one of 'left', 'right', 'top', 'bottom'")
        if plugin=='skimage':
            io.imsave(os.path.join(saveDir, n+'_cropped'+e), im, imagej=True, check_contrast=False)
        elif plugin=='aicsimageio':
            OmeTiffWriter.save(data=im, 
                            uri=os.path.join(saveDir, n+'_cropped'+e), 
                            physical_pixel_sizes=pixelSizes, 
                            dimension_order="ZYX")

if __name__=='__main__':
    p  = argparse.ArgumentParser()
    p.add_argument('-d', '--directory', type=str, default=os.getcwd(), help='Directory containing images to crop')
    p.add_argument('-l', '--level', type=int, default=3, help='Level of tile to crop')
    p.add_argument('-a', '--axis', type=str, default='x', help='Axis to crop from')
    p.add_argument('-c', '--cropDim', type=int, help='Number of pixels to crop out')
    p.add_argument('-s', '--side', type=str, default='left', help='Side to crop off of')
    args = p.parse_args()
    cropZeissTiles(args.directory, args.level, args.axis, args.cropDim, args.side)


