# imageProcessing
Image Processing Scripts

Primarily pre- and post-processing of TIF stacks for ClearMap. 

Full ClearMap repository found at: https://github.com/ChristophKirst/ClearMap

Email me at alexander.smith@mssm.edu with any questions.

Currently, my analysis pipeline looks like this:

1) clearMapStackCropping.py script - Currently I still open the autofluorescence files manually in ImageJ and save them as one large TIF stack. This basically ensures that I actually look at each brain to make sure nothing is super wrong with the imaging. I do check all of the signal channel images while tuning parameters later, but for the most part run the analysis blind, and use the same parameters for all brains, of course.

2) clearMap_Ilastik_parameters.py for machine learning classification, using a classifier that has previously been trained using Ilastik (ilastik.org) or clearMapWatershedParameters.py for detection of cells based on size/shape/intensity (better for regular circular c-Fos staining etc). Tune the parameters based on intensity values you check in the signal channel via ImageJ.

3) clearMap_Ilastik_process or clearMapWatershedProcess.py - Currently, the only function I use in these scripts is the detectCells(**ImageProcessingParameter) line. That generates a NumPy array with data for each animal, which is computationally intensive and has to be done for one animal at a time. After those arrays have been generated, I process them with:

4) clearMapBatchProcessing.py script runs all of the other steps from the process template (downsampling, registration, transformation, voxeliation, region counts) en masse for all animals at once. This saves a lot of time.

5) clearMapAnalyzeTemplate.py is pretty self-explanatory, once you have generated data for each animal, this script runs the groupwise comparisons.

6) clearMapAdvancedAnalysis_fromLauraDeNardo.py is a script that I was given from Laura DeNardo, as the name implies, and was used in her recent Nature Neuroscience paper. This script does some fun things like clustering, tSNE, PCA, LDA, correlational analysis, etc.
