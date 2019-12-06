#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  2 03:39:40 2018

@author: smith
"""

import os
import ClearMap.IO as io
import ClearMap.Settings as settings
import ClearMap.Visualization.Plot as plt
import ClearMap.ImageProcessing.BackgroundRemoval as bgr
import matplotlib as mplt



BaseDirectory = '/d2/studies/ClearMap/Opto_iDISCO/DP9_chR2_1-3x/'
filename = os.path.join(BaseDirectory, 'DP9_cfos_stack.ome.tif');


data = io.readData(filename, z = (650, 750));
# size restriction 
fileRange = 'x = (550, 650), y = (600, 700), z = (25,34)' 

""" note that the Z coordinate in fileRange variable are relative to the planes imported in the data variable.
So in this case, you would be imaging planes 675-684.
"""

plt.plotTiling(data, inverse = True, x = (550, 650), y = (600, 700), z = (25,34)) #
# background subtraction
dataBGR = bgr.removeBackground(data.astype('float'), size=(7,7), verbose = False, save = '7_background.tif');
dataBGR_write = plt.plotTiling(dataBGR, inverse = True, x = (550, 650), y = (600, 700), z = (25,34));
dataBGR_write = plt.overlayPoints(dataBGR, fileRange);
mplt.pyplot.savefig('dataBGR_write.tif')

io.writeData(os.path.join(BaseDirectory, 'background_15.tif'), dataBGR_write);
io.writeData(os.path.join(BaseDirectory, 'cells_check.tif'), data);


pointSource= os.path.join(BaseDirectory, FilteredCellsFile[0]);
data_write = plt.overlayPoints(filename, dataBGR_write, fileRange, pointColor = None);
io.writeData(os.path.join(BaseDirectory, 'cells_check.tif'), data);


#DoG Filter
from ClearMap.ImageProcessing.Filter.DoGFilter import filterDoG
dataDoG = filterDoG(dataBGR, size=(7,7,9), verbose = False);
plt.plotTiling(dataDoG, inverse = True, x = (600, 700), y = (600, 700), z = (1, 15));

#Find Extended Maxima
from ClearMap.ImageProcessing.MaximaDetection import findExtendedMaxima
dataMax = findExtendedMaxima(dataDoG, 
hMax = None, 
verbose = True, 
threshold = 500);
plt.plotOverlayLabel( dataDoG / dataDoG.max(), dataMax.astype('int'), x = (1750, 1850), y = (1050, 1150), z = (1, 9))




#Find Peak Intensity
findIntensityParameter = {
    "method" : 'Max',       # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
    "size"   : (7,7,5)      # (tuple)             size of the search box on which to perform the *method*
}
#Cell Shape Parameters
detectCellShapeParameter = {
    "threshold" : 500,     # (float or None)      threshold to determine mask. Pixels below this are background if None no mask is generated
    "save"      : None, # (str or None)        file name to save result of this operation if None dont save to file 
    "verbose"   : True      # (bool or int)        print / plot information about this step if None take intensities at the given pixels
}

# result files for cell coordinates (csv, vtk or ims)
SpotDetectionParameter = {
    "source" : cFosFile,
    "sink"   : (os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')),
    "detectSpotsParameter" : detectSpotsParameter
};
SpotDetectionParameter = joinParameter(SpotDetectionParameter, cFosFileRange)

ImageProcessingParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter);
#Cell detection:
################
detectCells(**ImageProcessingParameter);

#Thresholding: the threshold parameter is either intensity or size in voxel, depending on the chosen "row"
#row = (0,0) : peak intensity from the raw data
#row = (1,1) : peak intensity from the DoG filtered data
#row = (2,2) : peak intensity from the background subtracted data
#row = (3,3) : voxel size from the watershed
points, intensities = thresholdPoints(points, intensities, threshold = (20, 900), row = (3,3));
io.writePoints(FilteredCellsFile, (points, intensities));
# Transform point coordinates
#############################
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
io.writePoints(TransformedCellsFile, points);


#Detect Cell Coordinates
from ClearMap.ImageProcessing.MaximaDetection import findCenterOfMaxima
cells = findCenterOfMaxima(data, dataMax);
plt.plotOverlayPoints(data, cells, z = (1, 9))
print cells.shape

#Detect Cell Coordinates
from ClearMap.ImageProcessing.MaximaDetection import findCenterOfMaxima
cells = findCenterOfMaxima(data, dataMax);
print cells.shape

#Cell Shape Detection
from ClearMap.ImageProcessing.CellSizeDetection import detectCellShape
dataShape = detectCellShape(dataDoG, cells, threshold = 500);
plt.plotOverlayLabel(dataDoG / dataDoG.max(), dataShape, z = (1, 9))

#Loading the results:
points, intensities = io.readPoints(ImageProcessingParameter["sink"]);






