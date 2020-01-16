#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 10:36:50 2020

@author: smith



This file writes a TIF stack containing cell center coordinates transformed to the atlas, and then
isolates a region, and splits into two hemispheres. The result can be used as input for the clearMapSubregionParser script."
"""


import ClearMap.IO.IO as io
import os
import ClearMap.Visualization.Plot as plt
import ClearMap.Analysis.Label as lbl
import numpy as np

samples = ['IA1_RB', 'IA1_LT', 'IA1_LB', 
           'IA2_NP', 'IA2_RT', 'IA2_LT', 'IA2_LB']

batchDirectory = '/d2/studies/ClearMap/IA_iDISCO/batchFiles/'
batchList = []

for mouse in samples:
    paramFile = os.path.join(batchDirectory, 'parameter_file_' + mouse + '.py')
    batchList.append(paramFile)
    execfile(paramFile)

    baseDirectory = BaseDirectory

#baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/IA1_RB/'
    sampleName = mouse
    region = 'Caudoputamen'


    points = io.readPoints(TransformedCellsFile)
    data = plt.overlayPoints(AnnotationFile, points.astype(int), pointColor = None)
    data = data[:,:,:,1:]
    io.writeData(os.path.join(BaseDirectory, sampleName + '_Points_Transformed.tif'), data)

#If you are using the same annotation file for every sample, comment out lines 45-54 to drastically reduce run time
    label = io.readData(AnnotationFile)
    label = label.astype('int32')
    labelids = np.unique(label)

    outside = np.zeros(label.shape, dtype = bool);


    for l in labelids:
        if not (lbl.labelAtLevel(l, 6) == 672):
            outside = np.logical_or(outside, label == l);

#DP = 814 (level 6)
#MHb = 483 (level 7)
#Caudoputamen = 672 (level 6)
#Accumbens = 56
#CA3 = 463 (Level 8)
#Prelimbic = 972 (Level 6)

    heatmap = io.readData(os.path.join(baseDirectory, sampleName + '_Points_Transformed.tif'))
    heatmap[outside] = 0;


    Xmin = np.amin(np.nonzero(heatmap)[0])
    Xmax = np.amax(np.nonzero(heatmap)[0])
    Ymin = np.amin(np.nonzero(heatmap)[1])
    Ymax = np.amax(np.nonzero(heatmap)[1])

    heatmap_left = heatmap[Xmin-10:heatmap.shape[0]/2,Ymin-10:Ymax+10,:]
    heatmap_right = heatmap[heatmap.shape[0]/2:Xmax+10,Ymin-10:Ymax+10:,:]

    io.writeData(os.path.join(baseDirectory, sampleName + '_' + region + '_isolated_points_left.tif'), heatmap_left)
    io.writeData(os.path.join(baseDirectory, sampleName + '_' + region + '_isolated_points_right.tif'), heatmap_right)
