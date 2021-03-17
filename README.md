# imageProcessing
Image Processing Scripts

Primarily pre- and post-processing of TIF stacks for ClearMap. 

Full ClearMap repository found at: https://github.com/ChristophKirst/ClearMap

Email me at alexander.smith@mssm.edu with any questions.

Currently, my analysis pipeline looks like this:

1) clearMapStackCropping.py script - Takes an uncropped .ome.tif stack of autofluorescence data, detects brain edges, and crops +/- 50 pixels from edges of both autofluorescence and signal files.

2) clearMapWatershedParameters.py - Set the parameters for cell detection (as described in methods section).

3) clearMapWatershedProcess.py - Performs resampling, alignment to atlas, cell detection, voxelization, and region segmentation.

4) clearMapBatchProcessing.py script runs all steps from the process template other than cell detection (i.e. downsampling, registration, transformation, voxeliation, region counts) en masse for all animals at once. This saves a lot of time.

5) clearMapAnalyzeTemplate.py - Once data is generated for each animal, this script runs the groupwise comparisons.

6) clearMapAdvancedAnalysis_fromLauraDeNardo.py is a script that I was given from Laura DeNardo, as the name implies, and was used in her recent Nature Neuroscience paper. This script does some fun things like clustering, tSNE, PCA, LDA, correlational analysis, etc.
