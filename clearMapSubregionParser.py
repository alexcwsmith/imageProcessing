#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 17:30:50 2019

@author: smith

This file works with isolated points for a subregion that are created by the plotTransformedPointCenters 
script. The function of this script is to split points into sub-groups, which can be done iteratively to
improve performance and accuracy for irregular shaped structures.
"""

import ClearMap.IO.IO as io
import numpy as np
import pandas as pd
import os


samples = ['IA1_RT', 'IA1_RB', 'IA1_LT', 'IA1_LB', 
           'IA2_NP', 'IA2_RT', 'IA2_RB', 'IA2_LT', 'IA2_LB']

dataList = []

for mouse in samples:
    sampleName = mouse
    baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/' + sampleName
    
    ##IMPORT PREVIOUSLY PRE-PROCESSED HeatMap Data
    hemisphere = '_left'
    data = io.readData(os.path.join(baseDirectory, sampleName + '_Caudoputamen' + '_isolated_points' + hemisphere + '.tif'))
    points = np.nonzero(data)[:3]
    dfPoints = pd.DataFrame(points, index=['x', 'y', 'z']).T
    dfPoints.rename(columns={0: "x", 1: "y", 2: "z"})
    
    ###View Range:
    x_range=dfPoints.x.max() - dfPoints.x.min()
    y_range=dfPoints.y.max() - dfPoints.y.min()
    z_range=dfPoints.z.max() - dfPoints.z.min()
    print(x_range, y_range, z_range)
    
    #Bin Y axis
    dfPoints['y_bins']=pd.cut(dfPoints['y'], bins=2)
    dfPoints['y_bins'].value_counts()
    dfPoints_counts_y = dfPoints['y_bins'].value_counts()
    dfPoints_sorted_y = dfPoints_counts_y.sort_index()
    dfPoints_sorted_y.to_excel(os.path.join(baseDirectory, 'y_binned_striatum_2' + hemisphere + '.xlsx'))
    
    figY = dfPoints_sorted_y.sort_index().plot.bar(figsize=(20,10), label='mouse')
    figY.savefig(os.path.join(baseDirectory, 'Figy.png'))
    print(dfPoints_sorted_y)
    
    #2nd Iteration - Splits anterior half of striatum into 3 subregions
    firstHalf = dfPoints_sorted_y[0]
    dfPoints_s2 = dfPoints.sort_values(by=['y'])[:firstHalf]
    dfPoints_s2['y_iter2'] = pd.cut(dfPoints_s2['y'], bins=6)
    dfPoints_s2['y_iter2'].value_counts()
    dfPoints_s2.sort_values('y_iter2')
    dfPoints_s2_count = dfPoints_s2['y_iter2'].value_counts()
    dfPoints_s2_count.sort_index()
    dfPoints_s2_count.to_excel(os.path.join(baseDirectory, 'y_binned_striatum_iter2-split6' + hemisphere + '.xlsx'))
    dfPoints_s2_count.sort_index()
    figY2 = dfPoints_s2_count.sort_index().plot.line(figsize=(20,10), label = 'mouse')
    figY2.savefig(os.path.join(baseDirectory, 'Figy2.png'))
    
    
    ant0 = dfPoints_s2_count.sort_index()[0]
    ant1 = dfPoints_s2_count.sort_index()[1]
    mid0 = dfPoints_s2_count.sort_index()[2]
    mid1 = dfPoints_s2_count.sort_index()[3]
    post0 = dfPoints_s2_count.sort_index()[4]
    post1 = dfPoints_s2_count.sort_index()[5]
    postA = ant0 + ant1
    postM = postA+mid0+mid1
    
    #Check that bins were split correctly: 
    if firstHalf != postM+post0+post1:
        raise ValueError('Variables Not Equal!')
    else: print('Values Are Equal')
    
    #Split each of the 6 subregions into medial/lateral and count cells:
    
    dfPoints_ant0 = dfPoints_s2.sort_values(by=['y_iter2'])[:ant0]
    dfPoints_ant0['x_bins0'] = pd.cut(dfPoints_ant0['x'], bins=2)
    dfPoints_ant0_count = dfPoints_ant0['x_bins0'].value_counts()
    dfPoints_ant0_count.sort_index()
    
    dfPoints_ant1 = dfPoints_s2.sort_values(by=['y_iter2'])[ant0:ant0+ant1]
    dfPoints_ant1['x_bins1'] = pd.cut(dfPoints_ant1['x'], bins=2)
    dfPoints_ant1_count = dfPoints_ant1['x_bins1'].value_counts()
    dfPoints_ant1_count.sort_index()
    
    postA = ant0 + ant1
    
    dfPoints_mid0 = dfPoints_s2.sort_values(by=['y_iter2'])[postA:postA+mid0]
    dfPoints_mid0['x_bins2'] = pd.cut(dfPoints_mid0['x'], bins=2)
    dfPoints_mid0_count = dfPoints_mid0['x_bins2'].value_counts()
    dfPoints_mid0_count.sort_index()
    
    
    dfPoints_mid1 = dfPoints_s2.sort_values(by=['y_iter2'])[postA:postA+mid0+mid1]
    dfPoints_mid1['x_bins0'] = pd.cut(dfPoints_mid1['x'], bins=2)
    dfPoints_mid1_count = dfPoints_mid1['x_bins0'].value_counts()
    dfPoints_mid1_count.sort_index()
    
    postM = postA + mid0 + mid1
    
    dfPoints_post0 = dfPoints_s2.sort_values(by=['y_iter2'])[postM:postM+post0]
    dfPoints_post0['x_bins1'] = pd.cut(dfPoints_post0['x'], bins=2)
    dfPoints_post0_count = dfPoints_post0['x_bins1'].value_counts()
    dfPoints_post0_count.sort_index()
    
    dfPoints_post1 = dfPoints_s2.sort_values(by=['y_iter2'])[postM+post0:postM+post0+post1]
    dfPoints_post1['x_bins2'] = pd.cut(dfPoints_post1['x'], bins=2)
    dfPoints_post1_count = dfPoints_post1['x_bins2'].value_counts()
    dfPoints_post1_count.sort_index()
    
    if hemisphere == '_left':
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
    
    mouse_data = [mouse, aDLS_combined, mDLS_combined, pDLS_combined, aDMS_combined, mDMS_combined, pDMS_combined]
    #dataList.append(mouse)
    dataList.append(mouse_data)

allData = pd.DataFrame(data=dataList, columns=['mouse', 'aDLS', 'mDLS', 'pDLS', 'aDMS', 'mDMS', 'pDMS'])    
allData.to_excel('/d2/studies/ClearMap/IA_iDISCO/Striatum_Subregion_Counts' + hemisphere+ '_3bins.xlsx')

