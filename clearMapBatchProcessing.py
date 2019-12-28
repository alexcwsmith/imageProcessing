#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 14:18:16 2019

@author: smith
"""

import pandas as pd
import numpy as np
import os
import glob
import skimage

samples = ['IA1_RT', 'IA1_RB', 'IA1_LT', 'IA1_LB', 
           'IA2_NP', 'IA2_RT', 'IA2_RB', 'IA2_LT', 'IA2_LB']

batchDirectory = '/d2/studies/ClearMap/IA_iDISCO/batchFiles/'
batchList = []

for mouse in samples:
    paramFile = os.path.join(batchDirectory, 'parameter_file_' + mouse + '.py')
    paramList.append(paramFile)
    execfile(paramFile)
    
#Loading the results of detectCells func
    points, intensities = io.readPoints(ImageProcessingParameter["sink"]);
#Filtering (here by voxel size)
    points, intensities = thresholdPoints(points, intensities, threshold = (4,200), row = (3,3));
    io.writePoints(FilteredCellsFile, (points, intensities));
#Transform points to atlas
    points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
    points = resamplePoints(**CorrectionResamplingPointsParameter);
    points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
    CorrectionResamplingPointsInverseParameter["pointSource"] = points;
    points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
    RegistrationResamplingPointParameter["pointSource"] = points;
    points = resamplePoints(**RegistrationResamplingPointParameter);
    points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
    io.writePoints(TransformedCellsFile, points);



#Load Voxelized Data
    VoxelizationFile = os.path.join(BaseDirectory, 'points_voxelized.tif');

# Parameter to calculate the density of the voxelization
    voxelizeParameter = {
    #Method to voxelize
        "method" : 'Spherical', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
        "size" : (25,25,25),  

    # Voxelization weigths (e/g intensities)
        "weights" : None
        };

    points = io.readPoints(TransformedCellsFile)
    intensities = io.readPoints(FilteredCellsFile[1])

    #Without weights:
    vox = voxelize(points, AtlasFile, **voxelizeParameter);
    if not isinstance(vox, basestring):
        io.writeData(os.path.join(batchDirectory, 'cells_heatmap_vox25' + '_' + mouse + '.tif'), vox.astype('int32'));
#

#With weigths from the intensity file (here raw intensity):
    voxelizeParameter["weights"] = intensities[:,0].astype(float);
    vox = voxelize(points, AtlasFile, **voxelizeParameter);
    if not isinstance(vox, basestring):
        io.writeData(os.path.join(batchDirectory, 'cells_heatmap_weighted_vox25' + '_' + mouse + '.tif'), vox.astype('int32'));
#


heatmapList = ['/d2/studies/ClearMap/IA_iDISCO/parameterFiles/cells_heatmap_vox25_IA2_LB.tif',
               '/d2/studies/ClearMap/IA_iDISCO/parameterFiles/cells_heatmap_vox25_IA2_LT.tif',
               '/d2/studies/ClearMap/IA_iDISCO/parameterFiles/cells_heatmap_vox25_IA2_RB.tif',
               '/d2/studies/ClearMap/IA_iDISCO/parameterFiles/cells_heatmap_vox25_IA2_RT.tif',
               '/d2/studies/ClearMap/IA_iDISCO/parameterFiles/cells_heatmap_vox25_IA2_NP.tif',
               '/d2/studies/ClearMap/IA_iDISCO/parameterFiles/cells_heatmap_vox25_IA1_LB.tif',
               '/d2/studies/ClearMap/IA_iDISCO/parameterFiles/cells_heatmap_vox25_IA1_LT.tif',
               '/d2/studies/ClearMap/IA_iDISCO/parameterFiles/cells_heatmap_vox25_IA1_RB.tif',
               '/d2/studies/ClearMap/IA_iDISCO/parameterFiles/cells_heatmap_vox25_IA1_RT.tif']

names = ['IA2_LB', 'IA2_LT', 'IA2_RB', 'IA2_RT', 'IA2_NP', 'IA1_LB', 'IA1_LT', 'IA1_RB', 'IA1_RT']

for heatmap in heatmapList:
    filename = os.path.basename(heatmap)
    f, e = os.path.splitext(filename)
    heatmap = skimage.io.imread(heatmap)
    heatmap.resize(250,528,456)
    skimage.color.rgb2gray(heatmap)
    heatmap = skimage.img_as_float32(heatmap)
    skimage.io.imsave(os.path.join(batchDirectory, f + '_resized.tif'), heatmap)

    
