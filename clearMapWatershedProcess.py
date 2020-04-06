# -*- coding: utf-8 -*-
"""
Template to run the processing pipeline
"""
#load the parameters:
execfile('/d2/studies/ClearMap/IA_iDISCO/IA2_RT/parameter_file_IA2_RT_2.py')

#resampling operations:
#######################
#resampling for the correction of stage movements during the acquisition between channels:
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

execfile('/d2/studies/ClearMap/IA_iDISCO/IA2_RT/parameter_file_IA2_RT_2.py')

#Cell detection:
################
detectCells(**ImageProcessingParameter);

#Filtering of the detected peaks:
#################################
#Loading the results:
points, intensities = io.readPoints(ImageProcessingParameter["sink"]);

#Thresholding: the threshold parameter is either intensity or size in voxel, depending on the chosen "row"
#row = (0,0) : peak intensity from the raw data
#row = (1,1) : peak intensity from the DoG filtered data
#row = (2,2) : peak intensity from the background subtracted data
#row = (3,3) : voxel size from the watershed
points, intensities = thresholdPoints(points, intensities, threshold = (8,200), row = (3,3));
io.writePoints(FilteredCellsFile, (points, intensities));


## Check Cell detection (For the testing phase only, remove when running on the full size dataset)
#########################
import ClearMap.Visualization.Plot as plt;
pointSource= os.path.join(BaseDirectory, FilteredCellsFile[0]);
data = plt.overlayPoints(cFosFile, pointSource, pointColor = None, **cFosFileRange);
io.writeData(os.path.join(BaseDirectory, 'cells_check_striatum.tif'), data);


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




# Heat map generation
#####################
points = io.readPoints(TransformedCellsFile)
intensities = io.readPoints(FilteredCellsFile[1])

#Without weigths:
vox = voxelize(points, AtlasFile, **voxelizeParameter);
if not isinstance(vox, basestring):
  io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_vox15_2.tif'), vox.astype('int32'));
#

#With weigths from the intensity file (here raw intensity):
voxelizeParameter["weights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **voxelizeParameter);
if not isinstance(vox, basestring):
  io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted_vox15_2.tif'), vox.astype('int32'));
#




#Table generation:
##################
#With integrated weigths from the intensity file (here raw intensity):
ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
io.writeTable(os.path.join(BaseDirectory, 'Annotated_counts_intensities_2.csv'), table);

#Without weigths (pure cell number):
ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
io.writeTable(os.path.join(BaseDirectory, 'Annotated_counts_2.csv'), table);



#####################
#####################
#####################
#####################
#####################
#####################
