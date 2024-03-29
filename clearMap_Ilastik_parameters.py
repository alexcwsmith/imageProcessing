# -*- coding: utf-8 -*-
"""
Example script to set up the parameters for the image processing pipeline.

Note even though variables indicate 'signalFile', the signal channel can be whatever you are looking at.
Just too much work to update all of the variables in the package to indicate other than cFos.
"""

######################### Import modules

import os, numpy, math
os.chdir('/d1/software/ClearMap/')
import pandas as pd
from scipy.stats import zscore
from datetime import datetime as dt

import ClearMap.Settings as settings
import ClearMap.IO as io

from ClearMap.Alignment.Resampling import resampleData
from ClearMap.Alignment.Elastix import alignData, transformPoints
from ClearMap.ImageProcessing.CellDetection import detectCells
from ClearMap.Alignment.Resampling import resamplePoints, resamplePointsInverse
from ClearMap.Analysis.Label import countPointsInRegions
from ClearMap.Analysis.Voxelization import voxelize
from ClearMap.Analysis.Statistics import thresholdPoints
from ClearMap.Utils.ParameterTools import joinParameter
from ClearMap.Analysis.Label import labelToName


######################### Data parameters

#Directory to save all the results, usually containing the data for one sample
sampleName='BH2-4'
today = dt.now().strftime("%m-%d-%Y")
analysisLabel=today # this will be appended to all result file names
saveName = '_'.join([sampleName, analysisLabel])

testRegion='all'
signalChannel='tdTomato'
prefix = '_'.join([sampleName, signalChannel])

BaseDirectory = '/d2/studies/CM2/SmartSPIM/BH_Npas4_TRAP/'+sampleName+'_Stitched/'

#Data File and Reference channel File, usually as a sequence of files from the microscope
#Use \d{4} for 4 digits in the sequence for instance. As an example, if you have cfos-Z0001.ome.tif :
#os.path.join() is used to join the BaseDirectory path and the data paths:
signalFile = os.path.join(BaseDirectory, 'stacks', sampleName+'_Stitched_'+signalChannel+'.ome.tif')
AutofluoFile = os.path.join(BaseDirectory, 'stacks', sampleName+'Stitched_auto.ome.tif')

#Specify the range for the cell detection. This doesn't affect the resampling and registration operations
if testRegion=='all':   
    signalFileRange = {'x' : all, 'y' : all, 'z' : all};
elif testRegion=='pfc':
    signalFileRange = {'x' : (1225,1525), 'y' : (350,650), 'z' : (800,900)};
elif testRegion=='hpc':
    signalFileRange = {'x' : (2100,2400), 'y' : (1500,1800), 'z' : (700,800)};

#Resolution of the Raw Data (in um / pixel)OriginalResolution = (5.0, 5.0, 4.5)
OriginalResolution = (4.0, 4.0, 4.0)

#Orientation: 1,2,3 means the same orientation as the reference and atlas files.
#Flip axis with - sign (eg. (-1,2,3) flips x). 3D Rotate by swapping numbers. (eg. (2,1,3) swaps x and y)
FinalOrientation = (1,2,3)

#Resolution of the Atlas (in um/ pixel)
AtlasResolution = (20, 20, 20)

#Path to registration parameters and atlases
PathReg        = '/d2/studies/ClearMap/Parameter_files/Parameter_Files_Perens';
AtlasFile      = '/d2/studies/CM2/SmartSPIM/LSFM_Atlas_Cropped.ome.tif';
AnnotationFile = '/d2/studies/CM2/SmartSPIM/LSFM_Annotations_Cropped.ome.tif';
regionIndex = '/d2/studies/ClearMap/RegionID_Index.csv'
volumeIndex = '/d2/studies/ClearMap/RegionVolumeIndex.csv'

######################### Cell Detection Parameters using custom filters

#Spot detection method: faster, but optimised for spherical objects.
#You can also use "Ilastik" for more complex objects
ImageProcessingMethod = "Ilastik"

ilastikParameter = {
    "classifier" : "/d2/studies/CM2/SmartSPIM/BH_Npas4_TRAP/ilastik/BH_Npas4_TRAP_Classifier_Aug23.ilp",
    "classindex" : 0,
    "save"       : os.path.join(BaseDirectory, 'IlastikClasses/',sampleName+'_IlastikClasses_Z\d{4}.ome.tif'),      # (str or None)       file name to save result of this operation if None dont save to file 
    "verbose"    : True       # (bool or int)       print / plot information about this step
}


#If the maximum instensity is not at the gravity center of the object, look for a peak intensity around the center of mass. 
findCellIntensityParameter = {
    "method" : 'Max',       # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
    "verbose": True         # (bool or int)       print / plot information about this step
}


## Parameters for cell detection using spot detection algorithm 
detectCellsParameter = {
    "classifyCellsParameter"  : ilastikParameter,
    "findCellIntensityParameter"  : findCellIntensityParameter
}





#################### Heat map generation

##Voxelization: file name for the output:
VoxelizationFile = os.path.join(BaseDirectory, sampleName+'_points_voxelized.tif')

# Parameter to calculate the density of the voxelization
voxelizeParameter = {
    #Method to voxelize
    "method" : 'Spherical', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
    "size" : (15,15,15),  

    # Voxelization weigths (e/g intensities)
    "weights" : None
}


############################ Config parameters

#Processes to use for Resampling (usually twice the number of physical processors)
ResamplingParameter = {
    "processes": 18
}


#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes. Be careful of the memory footprint of each process!
    "processes" : 8,
   
    #chunk sizes: number of planes processed at once
    "chunkSizeMax" : 2000,
    "chunkSizeMin" : 200,
    "chunkOverlap" : 20,

    #optimize chunk size and number to number of processes to limit the number of cycles
    "chunkOptimization" : True,
    
    #increase chunk size for optimization (True, False or all = automatic)
    "chunkOptimizationSize" : all,
   
    "processMethod" : "parallel"
   }






######################## Run Parameters, usually you don't need to change those


### Resample Fluorescent and CFos images
# Autofluorescent cFos resampling for aquisition correction

ResolutionAffinesignalAutoFluo =  (16, 16, 16)

CorrectionResamplingParametersignal = ResamplingParameter.copy()

CorrectionResamplingParametersignal["source"] = signalFile
CorrectionResamplingParametersignal["sink"]   = os.path.join(BaseDirectory, sampleName+'_'+signalChannel+'_resampled.tif')
    
CorrectionResamplingParametersignal["resolutionSource"] = OriginalResolution
CorrectionResamplingParametersignal["resolutionSink"]   = ResolutionAffinesignalAutoFluo

CorrectionResamplingParametersignal["orientation"] = FinalOrientation
   
   
   
#Files for Auto-fluorescence for acquisition movements correction
CorrectionResamplingParameterAutoFluo = CorrectionResamplingParametersignal.copy();
CorrectionResamplingParameterAutoFluo["source"] = AutofluoFile;
CorrectionResamplingParameterAutoFluo["sink"]   = os.path.join(BaseDirectory, sampleName+'_autofluo_for_signal_resampled.tif');
   
#Files for Auto-fluorescence (Atlas Registration)
RegistrationResamplingParameter = CorrectionResamplingParameterAutoFluo.copy();
RegistrationResamplingParameter["sink"]            =  os.path.join(BaseDirectory, sampleName+'_autofluo_resampled.tif');
RegistrationResamplingParameter["resolutionSink"]  = AtlasResolution;
   

### Align cFos and Autofluo

CorrectionAlignmentParameter = {            
    #moving and reference images
    "movingImage" : os.path.join(BaseDirectory, sampleName+'_autofluo_for_signal_resampled.tif'),
    "fixedImage"  : os.path.join(BaseDirectory, sampleName+'_'+signalChannel+'_resampled.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(PathReg, 'Par0000affine_acquisition.txt'),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_signal_to_auto')
    }
  

### Align Autofluo and Atlas

CorrectionAlignmentParameter = {            
    #moving and reference images
    "movingImage" : os.path.join(BaseDirectory, sampleName+'_autofluo_for_signal_resampled.tif'),
    "fixedImage"  : os.path.join(BaseDirectory, sampleName+'_'+signalChannel+'_resampled.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(PathReg, 'Par0000affine_acquisition.txt'),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_signal_to_auto')
    }; 
  

### Align Autofluo and Atlas

#directory of the alignment result
RegistrationAlignmentParameter = CorrectionAlignmentParameter.copy();

RegistrationAlignmentParameter["resultDirectory"] = os.path.join(BaseDirectory, 'elastix_auto_to_atlas');
    
#moving and reference images
RegistrationAlignmentParameter["movingImage"]  = AtlasFile;
RegistrationAlignmentParameter["fixedImage"]   = os.path.join(BaseDirectory, sampleName+'_autofluo_resampled.tif');

#elastix parameter files for alignment
RegistrationAlignmentParameter["affineParameterFile"]  = os.path.join(PathReg, 'Par0000affine.txt');
RegistrationAlignmentParameter["bSplineParameterFile"] = os.path.join(PathReg, 'Par0000bspline.txt');



# result files for cell coordinates (csv, vtk or ims)
CellDetectionParameter = {
    "method" : ImageProcessingMethod,
    "source" : signalFile,
    "sink"   : (os.path.join(BaseDirectory, sampleName+'_cells-allpoints_'+today+'.npy'),  os.path.join(BaseDirectory,  sampleName+'_intensities-allpoints_'+today+'.npy')),
}
CellDetectionParameter = joinParameter(CellDetectionParameter, detectCellsParameter)
CellDetectionParameter = joinParameter(CellDetectionParameter, signalFileRange)

ImageProcessingParameter = joinParameter(StackProcessingParameter, CellDetectionParameter)

FilteredCellsFile = (os.path.join(BaseDirectory, sampleName+'_cells_'+today+'.npy'), os.path.join(BaseDirectory,  sampleName+'_intensities'+today+'.npy'))

TransformedCellsFile = os.path.join(BaseDirectory, sampleName+'_cells_transformed_to_Atlas_'+today+'.npy')

### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to referenece image size

CorrectionResamplingPointsParameter = CorrectionResamplingParametersignal.copy();
CorrectionResamplingPointsParameter["pointSource"] = os.path.join(BaseDirectory, sampleName+'_cells_'+today+'.npy');
CorrectionResamplingPointsParameter["dataSizeSource"] = signalFile;
CorrectionResamplingPointsParameter["pointSink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSizeSource"] = signalFile;
CorrectionResamplingPointsInverseParameter["pointSink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSizeSource"] = signalFile;
RegistrationResamplingPointParameter["pointSink"]  = None;

def printTime():
    print(dt.now().strftime('%H:%M:%S'))
    return(dt.now().strftime('%H:%M:%S'))

def writeRunParameters():
    date = dt.now().strftime("%Y-%m-%d")
    time = dt.now().strftime('%H:%M:%S')
    if ImageProcessingMethod=='ilastik':
        writeString = str(ilastikParameter)
    params = str(date + '\n' + time + '\n Sample: ' + str(sampleName) + '\n ImageProcessingMethod: ' + str(ImageProcessingMethod) + '\n Classifier: ' + str(ilastikParameter['classifier']))
    with open(os.path.join(BaseDirectory, 'runParams_'+date+'_'+time), 'w') as f:
        f.write(params)
        f.close()


