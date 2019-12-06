#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 16:43:27 2018

@author: smith
"""

#automated isolation

#import ClearMap.Analysis.Statistics as stat
import ClearMap.Analysis.Label as lbl
import ClearMap.IO.IO as io
import numpy, os
import ClearMap.Alignment.Resampling as rsp


baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/IA1_RT/'
sampleName = 'IA1_RT'
AnnotationFile = '/d2/studies/ClearMap/IA_iDISCO/IA1_RT/IA1_RT_annotations_horizontal_NEW.ome.tif'
region = 'Caudoputamen'



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
    if not (lbl.labelAtLevel(l, 7) == 672):
       outside = numpy.logical_or(outside, label == l);

#DP = 814 (level 6)
#MHb = 483 (level 7)
#Caudoputamen = 672 (level 7)
#Accumbens = 56
#CA3 = 463 (Level 8)
#Prelimbic = 972 (Level 6)

heatmap = io.readData(os.path.join(baseDirectory, 'cells_heatmap.tif'))
heatmap[outside] = 0;
    
io.writeData(os.path.join(baseDirectory, sampleName + '_' + region + '_isolated.tif'), heatmap)

#convert sagittal to coronal
io.writeData(os.path.join(baseDirectory, sampleName + '_' + region + '_isolated_coronal.tif'), rsp.sagittalToCoronalData(os.path.join(baseDirectory, sampleName + '_' + region + '_isolated.tif')));
    

