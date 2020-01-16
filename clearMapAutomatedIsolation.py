#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 16:43:27 2018

@author: smith

This script performs automated isolation of regions from transformed ClearMap data.
"""

#automated isolation

#import ClearMap.Analysis.Statistics as stat
import ClearMap.Analysis.Label as lbl
import ClearMap.IO.IO as io
import numpy, os
import ClearMap.Alignment.Resampling as rsp


AnnotationFile = '/d2/studies/ClearMap/IA_iDISCO/Annotations_Horizontal_8-302_Full.tif'

label = io.readData(AnnotationFile)
label = label.astype('int32')
labelids = numpy.unique(label)

outside = numpy.zeros(label.shape, dtype = bool);


"""
In order to find out the level to use, in console input:
>>> lbl.labelAtLevel(r, n)
where r is region ID, and n is level (usually start at 5), if the output is not the
region ID, increase n.
"""
for l in labelids:
    if not (lbl.labelAtLevel(l, 6) == 672):
       outside = numpy.logical_or(outside, label == l);

#DP = 814 (level 6)
#MHb = 483 (level 7)
#Caudoputamen = 672 (level 6)
#Accumbens = 56
#CA3 = 463 (Level 8)
#Prelimbic = 972 (Level 6)

samples = ['IA1_RT', 'IA1_RB', 'IA1_LT', 'IA1_LB', 
           'IA2_NP', 'IA2_RT', 'IA2_RB', 'IA2_LT', 'IA2_LB']

for mouse in samples:
    baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/' + mouse
    heatmap = io.readData(os.path.join(baseDirectory, 'cells_heatmap_vox15.tif'))
    if not os.path.exists(heatmap):
        raise ValueError('Data does not exist. Check naming convention')



for mouse in samples:
    baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/' + mouse
    sampleName = mouse
    region = 'Background'
    
    
    heatmap = io.readData(os.path.join(baseDirectory, 'cells_heatmap_vox15.tif'))
    heatmap[outside] = 0;
        
    io.writeData(os.path.join(baseDirectory, sampleName + '_' + region + '_isolated.tif'), heatmap)
    
    #convert sagittal to coronal
    io.writeData(os.path.join(baseDirectory, sampleName + '_' + region + '_isolated_coronal.tif'), rsp.sagittalToCoronalData(os.path.join(baseDirectory, sampleName + '_' + region + '_isolated.tif')));
        
        
