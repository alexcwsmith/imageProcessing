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


execfile('/d2/studies/ClearMap/IA_iDISCO/IA2_RB/parameter_file_IA2_RB.py')

baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/IA2_RB/'
sampleName = 'IA2_RB'
region = 'Caudoputamen'


points = io.readPoints(TransformedCellsFile)
data = plt.overlayPoints(AnnotationFile, points.astype(int), pointColor = None)
data = data[:,:,:,1:]
io.writeData(os.path.join(BaseDirectory, 'IA2_RB_Points_Transformed.tif'), data)

label = io.readData(AnnotationFile)
label = label.astype('int32')
labelids = np.unique(label)

outside = np.zeros(label.shape, dtype = bool);


"""
In order to find out the level to use, in console input:
>>> lbl.labelAtLevel(r, n)
where r is region ID, and n is level (usually start at 5), if the output is not the
region ID, increase n.
"""
for l in labelids:
    if not (lbl.labelAtLevel(l, 6) == 672):
       outside = np.logical_or(outside, label == l);

#DP = 814 (level 6)
#MHb = 483 (level 7)
#Caudoputamen = 672 (level 6)
#Accumbens = 56
#CA3 = 463 (Level 8)
#Prelimbic = 972 (Level 6)

heatmap = io.readData(os.path.join(baseDirectory, 'Points_Transformed.tif'))
heatmap[outside] = 0;


Xmin = np.amin(np.nonzero(heatmap)[0])
Xmax = np.amax(np.nonzero(heatmap)[0])
Ymin = np.amin(np.nonzero(heatmap)[1])
Ymax = np.amax(np.nonzero(heatmap)[1])

heatmap_left = heatmap[Xmin-10:heatmap.shape[0]/2,Ymin-10:Ymax+10,:]
heatmap_right = heatmap[heatmap.shape[0]/2:Xmax+10,Ymin-10:Ymax+10:,:]

io.writeData(os.path.join(baseDirectory, sampleName + '_' + region + '_isolated_points_left.tif'), heatmap_left)
io.writeData(os.path.join(baseDirectory, sampleName + '_' + region + '_isolated_points_right.tif'), heatmap_right)
