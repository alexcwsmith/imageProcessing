#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 14:14:35 2020

@author: smith
"""

import os
from skimage import io
from skimage import morphology

sampleName = 'F1_RT'
directory = '/d2/studies/ClearMap/FosTRAP_ChR2/' + sampleName + '/seg-ChR2/'
images = []
files = os.listdir(directory)
for file in files:
    if file.endswith('.tif'):
        fullpath = os.path.join(directory, file)
        images.append(fullpath)
images = sorted(images)
       
data= []
for image in images:
    im = io.imread(image)
    data.append(im)

arr = io.concatenate_images(data)

#Threshold & Skeletonize Images. This creates a weighted skeleton that retains model confidence information.
#Deleting variables at the end of each step is a memory-saving process that may not be necessary.
im2 = (arr > 0.2).astype('float32')
skel2 = morphology.skeletonize_3d(im2).astype('float32') * 0.2

im3 = (arr > 0.3).astype('float32')
skel3 = morphology.skeletonize_3d(im3).astype('float32') * 0.3
skel2_3 = skel2 + skel3
del im2; del skel2; del im3; del skel3

im4 = (arr > 0.4).astype('float32')
skel4 = morphology.skeletonize_3d(im4).astype('float32') * 0.4
skel2_4 = skel2_3 + skel4
del im4; del skel4; del skel2_3

im5 = (arr > 0.5).astype('float32')
skel5 = morphology.skeletonize_3d(im5).astype('float32') * 0.5
skel2_5 = skel2_4 + skel5
del im5; del skel5; del skel2_4

im6 = (arr > 0.6).astype('float32')
skel6 = morphology.skeletonize_3d(im6).astype('float32') * 0.6
skel2_6 = skel2_5 + skel6
del im6; del skel6; del skel2_5

im7 = (arr > 0.7).astype('float32')
skel7 = morphology.skeletonize_3d(im7).astype('float32') * 0.7
skel2_7 = skel2_6 + skel7
del im7; del skel7; del skel2_6

im8 = (arr > 0.8).astype('float32')
skel8 = morphology.skeletonize_3d(im8).astype('float32') * 0.8
skel2_8 = skel2_7 + skel8
del im8; del skel8; del skel2_7

im9 = (arr > 0.9).astype('float32')
skel9 = morphology.skeletonize_3d(im9).astype('float32') * 0.9
skel2_9 = skel2_8 + skel9
del im9; del skel9; del skel2_8

io.imsave(os.path.join(directory, '../' + sampleName + '_TrailMap_ThresholdedSkeleton.tif'), skel2_9)



