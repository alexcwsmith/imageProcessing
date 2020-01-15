#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 15:08:02 2020

@author: smith

This script finds the edges of the TIF stack and crops them to be +/- 50 pixels from those edges.


"""

import os
import numpy as np
import skimage

sampleID = 'ROC_13'
baseDirectory = '/d2/studies/ClearMap/ROC_iDISCO/' + sampleID
data = skimage.io.imread(os.path.join(baseDirectory, sampleID + '_auto_stack.ome.tif'))

X1 = 500
X2 = 1500
XminData = data[600:700,1250:1550,:X1] #note dimension indices are [z,y,x], these need to contain the widest part of your brain
XmaxData = data[600:700,1250:1550,X2:]
Xmin = np.nonzero(XminData[2]>2000)[1]
Xmin = (Xmin.min()-50)
Xmax = np.nonzero(XmaxData[2]<500)[1]
Xmax = ((Xmax.min())+50)
Xmax = Xmax + X2
Y1 = 500
Y2 = 2000
YminData = data[600:700,:Y1,850:1250] #again these dimension indices are [z,y,x], and should contain the widest part of your brain.
YmaxData = data[600:700,Y2:,850:1250]

Ymin = np.nonzero(YminData[1]>2000)[0]
Ymin = (Ymin.min()-50)
Ymax = np.nonzero(YmaxData[1]<500)[0]
Ymax = ((Ymax.min())+50)
Ymax = Ymax + Y2


dataCrop = data[Xmin:Xmax,Ymin:Ymax]
skimage.io.imsave(os.path.join(baseDirectory, sampleID + '_auto_cropped.ome.tif'), dataCrop)









