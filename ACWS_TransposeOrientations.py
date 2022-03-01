#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 11:57:33 2022

@author: smith
"""

import os
import numpy as np
from skimage import io
import argparse

def horizontalToCoronal(img_path):
    f, e = os.path.splitext(img_path)
    if f.endswith('.ome'):
        f = f.strip('.ome')
    save_path = f+'_Coronal.tif'
    im = io.imread(img_path)
    print("Image loaded with size " + str(im.shape))
    im_corr = im.transpose(1,2,0)
    io.imsave(save_path, im_corr, check_contrast=False)

def horizontalToSagittal(img_path):
    f, e = os.path.splitext(img_path)
    if f.endswith('.ome'):
        f = f.strip('.ome')
    save_path = f+'_Sagittal.tif'
    im = io.imread(img_path)
    print("Image loaded with size " + str(im.shape))
    im_sag = im.transpose(2,0,1)
    io.imsave(save_path, im_sag, check_contrast=False)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--img_path',type=str,default=None,help='Full path to image for resampling')
    p.add_argument('-s', '--source',type=str,default='horizontal',help='Orientation of source image; default horizontal.')
    p.add_argument('-t', '--target',type=str,default='coronal',help='Target orientation. Default coronal.')
    args = p.parse_args()
    args.source = args.source.lower()
    args.target = args.target.lower()
    if args.img_path==None:
        raise NameError("img_path argument must be full path to image")
    if not args.source=='horizontal':
        raise NameError("Sorry, at this time this script only supports horizontal source images")
    if args.target=='coronal':
        horizontalToCoronal(args.img_path)
    elif args.target=='sagittal':
        horizontalToSagittal(args.img_path)
