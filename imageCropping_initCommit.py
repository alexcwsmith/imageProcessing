#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:52:40 2019

@author: smith
"""


import os
import skimage
import glob
from skimage.external import tifffile

sampleID = 'N1_subset_test_auto'

path = '/d2/studies/ClearMap/Molly/N1_subset_test/auto/'
resultpath = os.path.join(path, 'skimage_rescaled/')
if not os.path.exists(resultpath):
    os.mkdir(resultpath)    
dirs = os.listdir(path)

def cropAndResample():
    for item in dirs:
        fullpath = os.path.join(path,item)
        if os.path.isfile(fullpath):
            im = skimage.io.imread(fullpath)
            imCrop = im[32:8736,833:6959] #replaced 1008:7040, [32:8736,1008:6784]
            filename = os.path.basename(fullpath)
            f, e = os.path.splitext(filename)
            rescaled = skimage.transform.rescale(imCrop, 0.5, anti_aliasing=False)
            rescaled16 = skimage.img_as_uint(rescaled)
            skimage.io.imsave(os.path.join(resultpath, f + '_skimage_rescaled.tif'), rescaled16)
    with tifffile.TiffWriter(path + sampleID + '.ome.tif', bigtiff=True, imagej=True) as stack:
        for filename in sorted(glob.glob(os.path.join(resultpath, '*.tif'))):
            stack.save(tifffile.imread(filename), photometric='minisblack')

cropAndResample()

sampleID = 'N1_subset_test_cfos'
path = '/d2/studies/ClearMap/Molly/N1_subset_test/cfos/'
resultpath = os.path.join(path, 'skimage_rescaled/')
if not os.path.exists(resultpath):
    os.mkdir(resultpath)    
dirs = os.listdir(path)

cropAndResample()

