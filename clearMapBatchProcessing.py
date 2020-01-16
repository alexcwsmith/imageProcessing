#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 14:18:16 2019

@author: smith

This script does all of the steps in the ClearMap process template, except for cell detection, for all samples in a batch.

In order to run this script, copy/paste parameter files for each mouse into a new directory called 'batchDirectory'.


"""


import os
import skimage

samples = ['IA1_RT', 'IA1_RB', 'IA1_LT', 'IA1_LB', 
           'IA2_NP', 'IA2_RT', 'IA2_RB', 'IA2_LT', 'IA2_LB']

batchDirectory = '/d2/studies/ClearMap/IA_iDISCO/analysisFiles/'
batchList = []

for mouse in samples:
    paramFile = os.path.join(batchDirectory, 'parameter_file_' + mouse + '.py')
    batchList.append(paramFile)
    execfile(paramFile)
    
    
#####THIS SECTION CAN BE RUN AFTER EXECFILE ABOVE IN ORDER TO OVERRIDE THE INDIVIDUAL PARAMETER FILES AND TEST NEW SETTINGS FOR ALL SAMPLES AT ONCE
    #UNCOMMENT TO USE
#    detectCellShapeParameter = {
#    "threshold" : 125,     # (float or None)      threshold to determine mask. Pixels below this are background if None no mask is generated
#    "save"      : None, #os.path.join(BaseDirectory, 'cellShape/cellShape\d{4}.ome.tif'), # (str or None)        file name to save result of this operation if None dont save to file 
#    "verbose"   : True      # (bool or int)        print / plot information about this step if None take intensities at the given pixels
#    }
    

#    removeBackgroundParameter = {
#        "size"    : (8,8),  # size in pixels (x,y) for the structure element of the morphological opening
#        "save"    : None, #os.path.join(BaseDirectory, 'Background/background\d{4}.ome.tif'), # file name to save result of this operation
#        "verbose" : True  # print / plot information about this step       
#    }
    
    #Difference of Gaussians filter: to enhance the edges. Useful if the objects have a non smooth texture (eg: amyloid deposits)
#    filterDoGParameter = {
#        "size"    : (4,4,4),        # (tuple or None)      size for the DoG filter in pixels (x,y,z) if None, do not correct for any background
#        "sigma"   : None,        # (tuple or None)      std of outer Gaussian, if None automatically determined from size
#        "sigma2"  : None,        # (tuple or None)      std of inner Gaussian, if None automatically determined from size
#        "save"    : None, #os.path.join(BaseDirectory, 'DoG/DoG\d{4}.ome.tif'),        # (str or None)        file name to save result of this operation if None dont save to file 
#        "verbose" : True         # (bool or int)        print / plot information about this step
#    }
    
#
####Resampling Data:
    resampleData(**CorrectionResamplingParameterCfos);
    resampleData(**CorrectionResamplingParameterAutoFluo);
    
    #Downsampling for alignment to the Atlas:
    resampleData(**RegistrationResamplingParameter);


    #Alignment operations:
    ######################
    #correction between channels:
    resultDirectory  = alignData(**CorrectionAlignmentParameter);
    
    #alignment to the Atlas:
    resultDirectory  = alignData(**RegistrationAlignmentParameter);


    
#Loading the results of detectCells func
    points, intensities = io.readPoints(ImageProcessingParameter["sink"]);
#Filtering (here by voxel size)
    points, intensities = thresholdPoints(points, intensities, threshold = (10,200), row = (3,3));
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
        "size" : (15,15,15),  

    # Voxelization weigths (e/g intensities)
        "weights" : None
        };

    points = io.readPoints(TransformedCellsFile)
    intensities = io.readPoints(FilteredCellsFile[1])

    #Without weights:
    vox = voxelize(points, AtlasFile, **voxelizeParameter);
    if not isinstance(vox, basestring):
        io.writeData(os.path.join(batchDirectory, 'cells_heatmap_vox15' + '_' + mouse + '.tif'), vox.astype('int32'));


#With weigths from the intensity file (here raw intensity):
#    voxelizeParameter["weights"] = intensities[:,0].astype(float);
#    vox = voxelize(points, AtlasFile, **voxelizeParameter);
#    if not isinstance(vox, basestring):
#        io.writeData(os.path.join(batchDirectory, 'cells_heatmap_weighted_vox15' + '_' + mouse + '.tif'), vox.astype('int32'));
#        
#    ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
#    table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
#    table["id"] = ids;
#    table["counts"] = counts;
#    table["name"] = labelToName(ids);
#    io.writeTable(os.path.join(batchDirectory, 'Annotated_counts_intensities' + '_' + mouse + '.csv'), table);

#Without weigths (pure cell number):
    ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
    table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
    table["id"] = ids;
    table["counts"] = counts;
    table["name"] = labelToName(ids);
    io.writeTable(os.path.join(batchDirectory, 'Annotated_counts' + '_' + mouse + '.csv'), table);


#

#### The below code resizes heatmaps so they are all identical, in case you did not use the same annotation file for all samples.
#### This code should be used with caution and outputs should be manually checked for correctness.


heatmapList = ['/d2/studies/ClearMap/IA_iDISCO/batchDirectory/cells_heatmap_vox25_IA2_LB.tif',
               '/d2/studies/ClearMap/IA_iDISCO/batchDirectory/cells_heatmap_vox25_IA2_LT.tif',
               '/d2/studies/ClearMap/IA_iDISCO/batchDirectory/cells_heatmap_vox25_IA2_RB.tif',
               '/d2/studies/ClearMap/IA_iDISCO/batchDirectory/cells_heatmap_vox25_IA2_RT.tif',
               '/d2/studies/ClearMap/IA_iDISCO/batchDirectory/cells_heatmap_vox25_IA2_NP.tif',
               '/d2/studies/ClearMap/IA_iDISCO/batchDirectory/cells_heatmap_vox25_IA1_LB.tif',
               '/d2/studies/ClearMap/IA_iDISCO/batchDirectory/cells_heatmap_vox25_IA1_LT.tif',
               '/d2/studies/ClearMap/IA_iDISCO/batchDirectory/cells_heatmap_vox25_IA1_RB.tif',
               '/d2/studies/ClearMap/IA_iDISCO/batchDirectory/cells_heatmap_vox25_IA1_RT.tif']

names = ['IA2_LB', 'IA2_LT', 'IA2_RB', 'IA2_RT', 'IA2_NP', 'IA1_LB', 'IA1_LT', 'IA1_RB', 'IA1_RT']

for heatmap in heatmapList:
    filename = os.path.basename(heatmap)
    f, e = os.path.splitext(filename)
    heatmap = skimage.io.imread(heatmap)
    heatmap.resize(250,528,456)
    skimage.color.rgb2gray(heatmap)
    heatmap = skimage.img_as_float32(heatmap)
    skimage.io.imsave(os.path.join(batchDirectory, f + '_resized.tif'), heatmap)

    
