#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:08:23 2021

@author: smith
"""
    
import os
from skimage import io
import glob
import numpy as np
import argparse

directory = '/d2/studies/ClearMap/LCT_SmartSPIM/ROC32_LCT/'
dataDir = os.path.join(directory, '642nm/')
sampleName = 'ROC32_SmartSPIM_642nm'
load_pattern='642nm/080120*.tif'
zindices=(2150,3300)
saveTIF=True
crop=False
memmap=True 
#crop={'width': 5200, 'height':6900, 'x':0, 'y':0}
savePlanes=True
returnArray=False

def filesToStack(directory, load_pattern, sampleName, zindices=None, crop=False, saveTIF=True, savePlanes=True, memmap=True, returnArray=False, imagej=True):
    files = glob.glob(os.path.join(directory, load_pattern))
    files = sorted(files)
    if zindices:
        files = files[zindices[0]:zindices[1]]
        print("Processing Z planes " + str(zindices[0]) + ' to ' + str(zindices[1]))        
    im = io.imread_collection(files)
    imc = im.concatenate()
    print("Loaded " + str(imc.shape[0]) + " images with XY dimensions " + str(imc.shape[2]) + "x" + str(imc.shape[1]))
    if crop:
        if not isinstance(crop, dict):
            raise ValueError("Crop argument must be dict of {width, height, x, y} - i.e. ImageJ-style crop coordinates")
        imc = imc[:,crop['y']:crop['y']+crop['height'],crop['x']:crop['x']+crop['width']]
        print("Cropped to XY dimensions " + str(crop['width']) +'x'+ str(crop['height']), "starting at X" + str(crop['x']) + ' Y'+str(crop['y']))
    if memmap:
        arrName = os.path.join(directory, sampleName+'_Z'+str(zindices[0])+'-Z'+str(zindices[1])+'.npy')
        np.save(arrName, imc)
        if not returnArray:
            del imc
        if saveTIF:
            fname, ext = os.path.splitext(arrName)
            if memmap:
                im = np.load(arrName, allow_pickle=True, mmap_mode='r+')
            if savePlanes:
                if not os.path.exists(fname + '/'):
                    print("Making directory " + fname + "/'")
                    os.mkdir(fname + '/')
                    saveDir=fname+'/'
                if imagej:
                    im = im[np.newaxis,:,np.newaxis,:,:]
                    for i in range(im.shape[1]):
                        if zindices:
                            if zindices[0]!=0:
                                plane=str(zindices[0]+i).zfill(4)
                                bname=os.path.basename(fname)
                                io.imsave(os.path.join(saveDir, bname+'_'+plane+'.tif'), im[:,i], check_contrast=False)
                        elif not zindices:
                            plane = str(i).zfill(4)
                            io.imsave(fname+'_'+plane+'.tif', im[i], check_contrast=False, plugin='tifffile', imagej=True)
                elif not imagej:
                    for i in range(im.shape[0]):
                        if zindices:
                            if zindices[0]!=0:
                                plane=str(zindices[0]+i).zfill(4)
                                bname=os.path.basename(fname)
                                io.imsave(os.path.join(saveDir, bname+'_'+plane+'.tif'), im[i], check_contrast=False)
                        elif not zindices:
                            plane = str(i).zfill(4)
                            io.imsave(fname+'_'+plane+'.tif', im[i], check_contrast=False, plugin='tifffile', imagej=True)
            elif not savePlanes:
                if imagej:
                    """imagej argument docs:
                        If True, write an ImageJ hyperstack compatible file. 
                        This format can handle data types uint8, uint16, or float32 and data shapes up to 
                        6 dimensions in TZCYXS order. RGB images (S=3 or S=4) must be uint8. 
                        ImageJ's default byte order is big endian but this implementation uses the 
                        system's native byte order by default. ImageJ does not support BigTIFF format or 
                        LZMA compression. The ImageJ file format is undocumented.
                    """
                    im = im[np.newaxis,:,np.newaxis,:,:]
                    io.imsave(fname + '_stack.tif', im, check_contrast=False, plugin='tifffile', imagej=True)
                elif not imagej:
                    io.imsave(fname + '_stack.tif', im, check_contrast=False, plugin='tifffile', imagej=False)
    elif not memmap:
        if saveTIF:
            io.imsave(os.path.join(directory, sampleName+'.tif'), imc)           
    if returnArray:
        return imc

filesToStack(directory, load_pattern, sampleName, zindices=zindices, crop=crop, saveTIF=True, savePlanes=True, imagej=True)

if __name__=='__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-d','--directory',type=str,default=os.getcwd(),help='Directory containing image files')
    p.add_argument('-lp','--load_pattern',type=str,default='*.tif',help='Filename pattern of files to load')
    p.add_argument('-s','--sampleName',type=str,help='Name to use for saving')
    p.add_argument('-z','--zindices',type=tuple,default=None,help='Z indices to use if creating a substack')
    p.add_argument('-c', '--crop',type=dict,default=None,help='{width, height, x, y} dict to crop images to (i.e. imageJ style cropping coordinates)')
    p.add_argument('--save',type=bool,default=True,help='Boolean whether to save TIF output')
    args = p.parse_args()
    filesToArray(args.directory, args.load_pattern, args.sampleName, zindices=args.zindices, crop=args.crop, saveTIF=args.saveTIF)
    