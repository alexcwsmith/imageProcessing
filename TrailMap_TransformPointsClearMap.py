#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 16:02:33 2020

@author: smith
"""

import imageio
import numpy as np
import os

sampleName = 'F1_RT'
BaseDirectory = '/d2/studies/ClearMap/FosTRAP_ChR2/' + sampleName
execfile('/d2/studies/ClearMap/FosTRAP_ChR2/' + sampleName + '/parameter_file_TrailMap_' + sampleName + '.py')


arr = imageio.mimread('/d2/studies/ClearMap/FosTRAP_ChR2/' + sampleName + '/' + sampleName + '_TrailMap_ThresholdedSkeleton.tif', memtest=False)
coords = np.nonzero(arr)
points = np.transpose(coords)
points = np.flip(points, axis=1)
np.save(os.path.join(BaseDirectory, 'cells.npy'), points)
del arr
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
io.writePoints(TransformedCellsFile, points);


# Heat map generation
#####################
points = io.readPoints(TransformedCellsFile)

#Without weigths:
vox = voxelize(points, AtlasFile, **voxelizeParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_' + sampleName + '_TrailMap.tif'), vox.astype('int32'));
#


ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None, collapse = True);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
io.writeTable(os.path.join(BaseDirectory, 'Annotated_counts_collapse_' + sampleName + '_TrailMap.csv'), table);

#Without weigths (pure cell number):
ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None, collapse = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
io.writeTable(os.path.join(BaseDirectory, 'Annotated_counts_' + sampleName + '_TrailMap.csv'), table);

