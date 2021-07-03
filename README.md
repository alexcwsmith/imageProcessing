# imageProcessing
Image Processing Scripts

Primarily pre- and post-processing of TIF stacks for ClearMap. 

Full ClearMap repository found at: https://github.com/ChristophKirst/ClearMap

Email me at alexander.smith@mssm.edu with any questions.

Currently, my analysis pipeline looks like this:

0) Optional - use clearMapStackCropping.py to detect brain edges & crop +/- 50 pixels, which may optimize atlas registration. Alternatively manually crop in ImageJ.

1) clearMapWatershedParameters.py - Set the parameters for cell detection (as described in methods section).

2) clearMapWatershedProcess.py - Performs resampling, alignment to atlas, cell detection, voxelization, and region segmentation.

3) clearMapBatchProcessing.py script runs all steps from the process template other than cell detection (i.e. downsampling, registration, transformation, voxeliation, region counts) en masse for all animals at once. This saves a lot of time.

4) clearMapAnalyzeTemplate.py - Once data is generated for each animal, this script runs the groupwise comparisons.

5) clearMapAdvancedAnalysis.py - Runs clustering & correlational analyses (thank you to Laura DeNardo for providing a script that served as the original backbone of this).
