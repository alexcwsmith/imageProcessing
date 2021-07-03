#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 17:30:50 2019

@author: smith

"""

#import ClearMap.IO.IO as io #Note for Nature Communications I revised this code to use scikit-image IO to prevent need to install the full ClearMap package
from skimage import io
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse

# directory = '/d2/studies/ClearMap/IA_iDISCO/NatComm_Supplement/']

# samples = ['IA1_LB', 'IA1_LT', 'IA1_NP', 'IA1_RB', 'IA1_RT', 
#            'IA2_LB', 'IA2_LT', 'IA2_NP', 'IA2_RB', 'IA2_RT']

def parseSubregions(directory, hemisphere, samples, save=False):
    """Docstring:
        This function works with cell center coordinates for an isolated Allen Brain Atlas (ABA) 
	parent region that are created ClearMap overlayPoints & region isolation functions (Renier et al., 2016). 
        The function of this script is to split points into sub-groups along A/P or M/L axes, 
        allowing quantification of points in distinct subregions.

        For the purposes of Smith et al., currently under review at Nature Communications, we focus on the
        Dorsal Striatum (DS; i.e. Caudoputamen in ABA, region 672). We first split the DS along the A/P axis
        into two segments. The posterior segment is the tail of the caudate, and is not examined here (as is common in the field).
        Because the borders of anterior segment of the DS has significant curvature along the AP axis (moving medial as it moves anterior),
        split this segment into 6 bins in the A/P axis. The first two bins are considered anterior, the middle two are considered middle, 
	and final two are considered posterior. For each of these 3 A/P bins, we then split into two bins along the X axis in order to 
	quantify cells in the medial/lateral subdivisions.
        
        Parameters
        ----------
        directory: string, path to parent directory containing subdirectories for each sample
        samples: list of strings, sample names (needs to be consistent with folder and file names within directory)
        hemisphere: string, 'left' or 'right'
        
    """
    if hemisphere != 'right' and hemisphere != 'left':
        raise NameError("Invalid hemisphere name, must be 'left' or 'right'")
    if type(samples) != list:
        samples = [samples]
    
    dataList = []

    for sampleName in samples:
        baseDirectory = os.path.join(directory, sampleName)        
        ##Import previously generated points data isolated from the Caudoputamen
  #      data = io.readData(os.path.join(baseDirectory, sampleName + '_Caudoputamen' + '_isolated_points' + hemisphere + '.tif'))
        data = io.imread(os.path.join(baseDirectory, sampleName + '_Caudoputamen' + '_isolated_points_' + hemisphere + '.tif'))
        data = data.transpose(2,1,0,3) #not necessary if using ClearMap IO
        points = np.nonzero(data)[:3]
        dfPoints = pd.DataFrame(points, index=['x', 'y', 'z']).T
        dfPoints.rename(columns={0: "x", 1: "y", 2: "z"})
        
        ###View Range:
        x_range=dfPoints.x.max() - dfPoints.x.min()
        y_range=dfPoints.y.max() - dfPoints.y.min()
        z_range=dfPoints.z.max() - dfPoints.z.min()
    #    print(x_range, y_range, z_range)
        
        #Bin Y axis into 2 bins - the back half of the caudoputamen in the ABA is the tail of the caudate, and is not examined here (as it is typically not in literature)
        dfPoints['y_bins']=pd.cut(dfPoints['y'], bins=2)
        dfPoints['y_bins'].value_counts()
        dfPoints_counts_y = dfPoints['y_bins'].value_counts()
        dfPoints_sorted_y = dfPoints_counts_y.sort_index()
        
        figY = plt.figure(figsize=(20,10))
        figY = dfPoints_sorted_y.sort_index().plot.bar(figsize=(20,10), label=sampleName)
        
        #2nd Iteration - Splits anterior half of striatum into 6 subregions.
        #Splitting into more bins helps accurately parse medial/lateral
            #divisions along the irregular shape of the anterior/posterior axis.
        firstHalf = dfPoints_sorted_y[0]
        dfPoints_s2 = dfPoints.sort_values(by=['y'])[:firstHalf]
        dfPoints_s2['y_iter2'] = pd.cut(dfPoints_s2['y'], bins=6)
        dfPoints_s2['y_iter2'].value_counts()
        dfPoints_s2.sort_values('y_iter2')
        dfPoints_s2_count = dfPoints_s2['y_iter2'].value_counts()
        dfPoints_s2_count.sort_index(inplace=True)
        
        ant = dfPoints_s2_count[0:2]
        mid = dfPoints_s2_count[2:4]
        post = dfPoints_s2_count[4:6]
        
        ant0 = dfPoints_s2_count[0]
        ant1 = dfPoints_s2_count[1]
        mid0 = dfPoints_s2_count[2]
        mid1 = dfPoints_s2_count[3]
        post0 = dfPoints_s2_count[4]
        post1 = dfPoints_s2_count[5]
        postA = ant0 + ant1 #starting position for everything posterior to anterior subregion.
        postM = postA+mid0+mid1 #starting position for everything posterior to middle subregion.
        
        #Check that bins were split correctly: 
        if firstHalf != postM+post0+post1:
            raise ValueError('Variables Not Equal!')
        else: print('Validated correct split')
        
        #Split each of the 6 subregions into medial/lateral and count cells:
        
        dfPoints_ant0 = dfPoints_s2.sort_values(by=['y_iter2'])[:ant0]
        dfPoints_ant0['x_bins0'] = pd.cut(dfPoints_ant0['x'], bins=2)
        dfPoints_ant0_count = dfPoints_ant0['x_bins0'].value_counts()
        dfPoints_ant0_count.sort_index(inplace=True)
        
        dfPoints_ant1 = dfPoints_s2.sort_values(by=['y_iter2'])[ant0:ant0+ant1]
        dfPoints_ant1['x_bins1'] = pd.cut(dfPoints_ant1['x'], bins=2)
        dfPoints_ant1_count = dfPoints_ant1['x_bins1'].value_counts()
        dfPoints_ant1_count.sort_index(inplace=True)
        
        postA = ant0 + ant1
        
        dfPoints_mid0 = dfPoints_s2.sort_values(by=['y_iter2'])[postA:postA+mid0]
        dfPoints_mid0['x_bins2'] = pd.cut(dfPoints_mid0['x'], bins=2)
        dfPoints_mid0_count = dfPoints_mid0['x_bins2'].value_counts()
        dfPoints_mid0_count.sort_index(inplace=True)
        
        
        dfPoints_mid1 = dfPoints_s2.sort_values(by=['y_iter2'])[postA:postA+mid0+mid1]
        dfPoints_mid1['x_bins0'] = pd.cut(dfPoints_mid1['x'], bins=2)
        dfPoints_mid1_count = dfPoints_mid1['x_bins0'].value_counts()
        dfPoints_mid1_count.sort_index(inplace=True)
        
        postM = postA + mid0 + mid1
        
        dfPoints_post0 = dfPoints_s2.sort_values(by=['y_iter2'])[postM:postM+post0]
        dfPoints_post0['x_bins1'] = pd.cut(dfPoints_post0['x'], bins=2)
        dfPoints_post0_count = dfPoints_post0['x_bins1'].value_counts()
        dfPoints_post0_count.sort_index(inplace=True)
        
        dfPoints_post1 = dfPoints_s2.sort_values(by=['y_iter2'])[postM+post0:postM+post0+post1]
        dfPoints_post1['x_bins2'] = pd.cut(dfPoints_post1['x'], bins=2)
        dfPoints_post1_count = dfPoints_post1['x_bins2'].value_counts()
        dfPoints_post1_count.sort_index(inplace=True)
        
        
        #Concatenate the bins into subdivisions:
        if hemisphere == 'left':
            aDLS_combined = dfPoints_ant0_count[0] + dfPoints_ant1_count[0]
            mDLS_combined = dfPoints_mid0_count[0] + dfPoints_mid1_count[0]
            pDLS_combined = dfPoints_post0_count[0] + dfPoints_post1_count[0]
            aDMS_combined = dfPoints_ant0_count[1] + dfPoints_ant1_count[1]
            mDMS_combined = dfPoints_mid0_count[1] + dfPoints_mid1_count[1]
            pDMS_combined = dfPoints_post0_count[1] + dfPoints_post1_count[1]
        else:
            aDLS_combined = dfPoints_ant0_count[1] + dfPoints_ant1_count[1]
            mDLS_combined = dfPoints_mid0_count[1] + dfPoints_mid1_count[1]
            pDLS_combined = dfPoints_post0_count[1] + dfPoints_post1_count[1]
            aDMS_combined = dfPoints_ant0_count[0] + dfPoints_ant1_count[0]
            mDMS_combined = dfPoints_mid0_count[0] + dfPoints_mid1_count[0]
            pDMS_combined = dfPoints_post0_count[0] + dfPoints_post1_count[0]
        
        mouse_data = [sampleName, aDLS_combined, mDLS_combined, pDLS_combined, aDMS_combined, mDMS_combined, pDMS_combined]
        #dataList.append(mouse)
        dataList.append(mouse_data)
    
    allData = pd.DataFrame(data=dataList, columns=['mouse', 'aDLS', 'mDLS', 'pDLS', 'aDMS', 'mDMS', 'pDMS'])  
    if save:
        allData.to_csv(os.path.join(directory, sampleName + '_Striatum_Subregion_Counts_' + hemisphere+ '.csv'))
    print(allData)
    return(allData)

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', default='./')
    parser.add_argument('--samples')
    parser.add_argument('--hemisphere', default='left')
    parser.add_argument('--save', default=False)
    args = parser.parse_args()

    parseSubregions(args.directory, args.hemisphere, args.samples, args.save)
