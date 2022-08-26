# -*- coding: utf-8 -*-
"""
Template to run the processing pipeline
"""

#load the parameters:
sampleName = 'BH2-4'
signalChannel='Npas4'
testing=False

execfile('/d2/studies/CM2/SmartSPIM/BH_Npas4_TRAP/' + sampleName + '_Stitched/clearMap_Ilastik_parameters_'+sampleName+'.py')


#resampling operations:
#######################
#resampling for the correction of stage movements during the acquisition between channels:
resampleData(**CorrectionResamplingParameterCfos)
resampleData(**CorrectionResamplingParameterAutoFluo)

#Downsampling for alignment to the Atlas:
resampleData(**RegistrationResamplingParameter)

#Alignment operations:
######################
#correction between channels:
resultDirectory  = alignData(**CorrectionAlignmentParameter)

#alignment to the Atlas:
resultDirectory  = alignData(**RegistrationAlignmentParameter)

#execfile('/d2/studies/ClearMap/ROC_iDISCO/ROC_9/parameter_file_Ilastik_ROC9.py')

#Cell detection:
################
printTime()
detectCells(**ImageProcessingParameter)
printTime()
writeRunParameters()

#Filtering of the detected peaks:
#################################
#Loading the results:

points, intensities = io.readPoints(ImageProcessingParameter["sink"]);

import matplotlib.pyplot as mplt
print(str(points.shape[0]) + " points detected")
topPerc = numpy.percentile(intensities[:,1], 99.9)
fig = mplt.figure()
vox = intensities[intensities[:,1]<topPerc][:,1]
mplt.hist(vox, bins=10)
#mplt.show()
fig.savefig(os.path.join(BaseDirectory, saveName+'_voxelSizeHistogram_'+testRegion+'.png'))


#Thresholding: the threshold parameter is either intensity or size in voxel, depending on the chosen "row"
#row = (0,0) : peak intensity from the raw data
#row = (1,1) : peak intensity from the DoG filtered data
#row = (2,2) : peak intensity from the background subtracted data
#row = (3,3) : voxel size from the watershed
print("Unfiltered # of points: " + str(points.shape[0]))
minSize=8
maxSize=100
print("Filtering out " + str(intensities[intensities[:,1]<minSize].shape[0]) + " points smaller than " + str(minSize) + " pixels")
print("Filtering out " + str(intensities[intensities[:,1]>maxSize].shape[0]) + " points larger than " + str(maxSize) + " pixels")
points, intensities = thresholdPoints(points, intensities, threshold = (minSize, maxSize), row = (1,1));
io.writePoints(FilteredCellsFile, (points, intensities));

now = dt.now().strftime('%m-%d-%y_%H:%M:%S')
params = "Min size: " + str(minSize) + "\n" + "Max size: " + str(maxSize)
with open(os.path.join(BaseDirectory, 'filterParams_'+now+'.txt'), 'w') as f:
    f.write(params)
    f.close()
    

## Check Cell detection (For the testing phase only, remove when running on the full size dataset)
#######################
if testing:
    import ClearMap.Visualization.Plot as plt;
    pointSource= os.path.join(BaseDirectory, FilteredCellsFile[0]);
    data = plt.overlayPoints(signalFile, pointSource, pointColor = None, **signalFileRange);
    io.writeData(os.path.join(BaseDirectory, sampleName+'cells_check_'+testRegion+'_'+str(minSize)+'-'+str(maxSize)+'px.tif'), data);


# Transform point coordinates
#############################
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"])
points = resamplePoints(**CorrectionResamplingPointsParameter)
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None)
CorrectionResamplingPointsInverseParameter["pointSource"] = points
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter)
RegistrationResamplingPointParameter["pointSource"] = points
points = resamplePoints(**RegistrationResamplingPointParameter)
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None)
io.writePoints(TransformedCellsFile, points)


# Heat map generation
#####################
points = io.readPoints(TransformedCellsFile)
intensities = io.readPoints(FilteredCellsFile[1])

#Without weigths:
vox = voxelize(points, AtlasFile, **voxelizeParameter)
io.writeData(os.path.join(BaseDirectory, saveName + '_cells_heatmap.tif'), vox.astype('int32'))
#
#With weigths from the intensity file (here raw intensity):
voxelizeParameter["weights"] = intensities[:,0].astype(float)
vox = voxelize(points, AtlasFile, **voxelizeParameter)
io.writeData(os.path.join(BaseDirectory, saveName + '_cells_heatmap_weighted.tif'), vox.astype('int32'))
#


###Counts analysis
regs = pd.read_csv(regionIndex, index_col=0)
vols = pd.read_csv(volumeIndex, index_col=0)


#Without weigths (pure cell number):
ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None)
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids
table["counts"] = counts
table["name"] = labelToName(ids)
io.writeTable(os.path.join(BaseDirectory, saveName + '_Annotated_counts.csv'), table)

df_counts = pd.DataFrame(table)
df_counts.set_index('id',inplace=True)
df_counts['Volume']=vols['Volume']
df_counts.dropna(how='any',axis=0,inplace=True)
df_counts['Count_Density']=df_counts['counts']/df_counts['Volume']
df_counts.sort_values(by='Count_Density',ascending=False,inplace=True)
df_counts = df_counts.loc[df_counts['Volume']!=1169].copy()
df_counts = df_counts.loc[df_counts['name']!='fiber tracts'].copy()
df_counts['Zscore']=zscore(df_counts['Count_Density'])
df_counts.to_csv(os.path.join(BaseDirectory, saveName + '_Annotated_Counts_Processed.csv'))

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None, collapse=True)
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids
table["counts"] = counts
table["name"] = labelToName(ids)
io.writeTable(os.path.join(BaseDirectory, saveName + '_Annotated_counts_collapse.csv'), table)

df_counts_collapse = pd.DataFrame(table)
df_counts_collapse.set_index('id',inplace=True)
df_counts_collapse['Volume']=vols['Volume']
df_counts_collapse.dropna(how='any',axis=0,inplace=True)
df_counts_collapse['Count_Density']=df_counts_collapse['counts']/df_counts_collapse['Volume']
df_counts_collapse.sort_values(by='Count_Density',ascending=False,inplace=True)
df_counts_collapse = df_counts_collapse.loc[df_counts_collapse['Volume']!=1169].copy()
df_counts_collapse = df_counts_collapse.loc[df_counts_collapse['name']!='fiber tracts'].copy()
df_counts_collapse['Zscore']=zscore(df_counts_collapse['Count_Density'])
df_counts_collapse.to_csv(os.path.join(BaseDirectory, saveName + '_Annotated_Counts_Collapsed_Processed.csv'))
print(str(df_counts_collapse['counts'].sum()) + " total cells")
