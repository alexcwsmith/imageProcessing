#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 16:24:07 2022

@author: smith
"""

import h5py
import numpy as np
from skimage import io
import argparse
import os

def TIFtoH5(img_path, key='Data'):
    f, e = os.path.splitext(img_path)
    if f.endswith('.ome'):
        f = f.strip('.ome')
    save_path = f+'.h5'
    hf = h5py.File(save_path, 'a')
    im = io.imread(img_path)
    print("Image loaded with size " + str(im.shape))
    dset = hf.create_dataset(key, data=im)
    hf.close()
    print("Wrote h5 file with key '" + str(key) + "' to " + str(save_path))

def h5toTIF(img_path, key='Data'):
    if not img_path.endswith('.h5'):
        raise NameError("Only .h5 files supported at this time.")
    f, e = os.path.splitext(img_path)
    if f.endswith('.lux'):
        f = f.strip('.lux')
    save_path = f+'.tif'
    hf = h5py.File(img_path, 'r')
    data = hf.get(key)
    if not data:
        print("Invalid h5 key, detected keys are: " + str(hf.keys()))
    print("Image loaded with size " + str(data.shape))
    io.imsave(save_path, np.array(data), bigtiff=True, check_contrast=False)
    print("Wrote tif file to " + str(save_path))
    
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--img_path',type=str,default=None,help='Full path to image for conversion')
    p.add_argument('-k', '--key', type=str,default='Data',help="Key for accessing data in h5 file. Default 'Data'")
    args = p.parse_args()
    if args.img_path==None:
        raise NameError("img_path argument must be full path to image")
    if args.img_path.endswith('.h5'):
        h5toTIF(args.img_path, key=args.key)
    elif args.img_path.endswith('.tif'):
        TIFtoH5(args.img_path, key=args.key)
