#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 15:08:02 2020

@author: smith

This script finds the edges of the TIF stack and crops them to be +/- 50 pixels from those edges.
It writes a .xlsx file with the new coordinates in ImageJ format. While it is possible to change the threshold
values established here to use this script to crop signal channels, because of greater variability in signal 
than autofluorescence, it is best to use this script on the autofluorescence file, then manually crop the 
data image stack with info from the .xlsx output (using Imagej --> edit --> selection --> specify)
"""

import os
import numpy as np
import pandas as pd
import skimage


sampleID = 'ROC_17'
baseDirectory = '/d2/studies/ClearMap/ROC_iDISCO/' + sampleID
data = skimage.io.imread(os.path.join(baseDirectory, sampleID + '_auto_stack.ome.tif'))

X1 = 500
X2 = 1500
XminData = data[600:700,1250:1550,:X1] #note dimension indices are [z,y,x], these need to contain the widest part of your brain, but should be limited to reduce computation
XmaxData = data[600:700,1250:1550,X2:]
Xmin = np.nonzero(XminData[2]>2000)[1]
Xmin = (Xmin.min()-50)
Xmax = np.nonzero(XmaxData[2]<200)[1]
Xmax = ((Xmax.min())+50)
Xmax = Xmax + X2
Y1 = 500
Y2 = 2000
YminData = data[600:700,:Y1,850:1250] #again these dimension indices are [z,y,x], and should contain the widest part of your brain.
YmaxData = data[600:700,Y2:,850:1250]

Ymin = np.nonzero(YminData[1]>2000)[0]
Ymin = (Ymin.min()-50)
Ymax = np.nonzero(YmaxData[1]<200)[0]
Ymax = ((Ymax.min())+50)
Ymax = Ymax + Y2

Width = (Xmax-Xmin).astype(int)
Height = Ymax-Ymin.astype(int)

if Width < 1500:
    raise ValueError('Small width dimension detected, proceed with caution.')
if Height < 1500:
    raise ValueError('Small height dimension detected, proceed with caution.')
else: print(Width, Height)

if Height > 2560:
    YMin = 0
    Ymax = 2560

Coords = pd.DataFrame(data=[Width,Height,Xmin,Ymin], index=['Width', 'Height', 'Xmin', 'Ymin'], columns=['Coordinate(pixel)'])
Coords.to_excel(os.path.join(baseDirectory, 'Cropping_coordinates.xlsx'))

dataCrop = data[:data.shape[0],Ymin:Ymax,Xmin:Xmax,]
skimage.io.imsave(os.path.join(baseDirectory, sampleID + '_auto_cropped.ome.tif'), dataCrop)

sampleID = 'ROC_17'
path = '/d2/studies/ClearMap/ROC_iDISCO/' + sampleID + '/SPmR/'
resultpath = os.path.join(path, 'SPmR_cropped/')
if not os.path.exists(resultpath):
    os.mkdir(resultpath)    
dirs = os.listdir(path)

def cropAndResample():
    for item in dirs:
        fullpath = os.path.join(path,item)
        if os.path.isfile(fullpath):
            im = skimage.io.imread(fullpath)
            imCrop = im[Ymin:Ymax,Xmin:Xmax] #replaced 1008:7040, [32:8736,1008:6784]
            filename = os.path.basename(fullpath)
            f, e = os.path.splitext(filename)
#            rescaled = skimage.transform.rescale(imCrop, 0.5, anti_aliasing=False)
#            rescaled16 = skimage.img_as_uint(rescaled)
            skimage.io.imsave(os.path.join(resultpath, f + '_cropped.tif'), imCrop)
#    with tifffile.TiffWriter(path + sampleID + '.ome.tif', bigtiff=True, imagej=True) as stack:
#        for filename in sorted(glob.glob(os.path.join(resultpath, '*.tif'))):
#            stack.save(tifffile.imread(filename), photometric='minisblack')

cropAndResample()








