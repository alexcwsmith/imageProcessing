#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 17:30:50 2019

@author: smith
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


#sampleName = 'IA2_LT'
#baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/' + sampleName
#sampleName = mouse
##IMPORT PREVIOUSLY PRE-PROCESSED POINTS DATA
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
    dfPoints_sorted_y.to_excel(os.path.join(baseDirectory, 'y_binned_striatum_Jan11' + hemisphere + '.xlsx'))
    
    #figY = dfPoints_sorted_y.plot.bar(figsize=(20,10))
    print(dfPoints_sorted_y)
    
    #2nd Iteration - Splits anterior half of striatum into 3 subregions
    firstHalf = dfPoints_sorted_y[0]
    dfPoints_ant = dfPoints.sort_values(by=['y'])[:firstHalf]
    dfPoints_ant['y_iter2'] = pd.cut(dfPoints_ant['y'], bins=3)
    dfPoints_ant['y_iter2'].value_counts()
    dfPoints_ant.sort_values('y_iter2')
    dfPoints_ant_count = dfPoints_ant['y_iter2'].value_counts()
    dfPoints_ant_count.sort_index()
    
    ant0 = dfPoints_ant_count.sort_index()[0]
    ant1 = dfPoints_ant_count.sort_index()[1]
    ant2 = dfPoints_ant_count.sort_index()[2:]
    ant2 = ant2[0]
    
    #Check that bins were split correctly: 
    if firstHalf != ant0+ant1+ant2:
        raise ValueError('Variables Not Equal!')
    else: print('Values Are Equal')
    
    #Split each of the 3 subregions into medial/lateral and count cells:  
    dfPoints_ant0 = dfPoints_ant.sort_values(by=['y_iter2'])[:ant0]
    dfPoints_ant0['x_bins0'] = pd.cut(dfPoints_ant0['x'], bins=2)
    dfPoints_ant0_count = dfPoints_ant0['x_bins0'].value_counts()
    dfPoints_ant0_count.sort_index()
    
    dfPoints_ant1 = dfPoints_ant.sort_values(by=['y_iter2'])[ant0:ant0+ant1]
    dfPoints_ant1['x_bins1'] = pd.cut(dfPoints_ant1['x'], bins=2)
    dfPoints_ant1_count = dfPoints_ant1['x_bins1'].value_counts()
    dfPoints_ant1_count.sort_index()
    
    dfPoints_ant2 = dfPoints_ant.sort_values(by=['y_iter2'])[ant0+ant1:]
    dfPoints_ant2['x_bins2'] = pd.cut(dfPoints_ant2['x'], bins=2)
    dfPoints_ant2_count = dfPoints_ant2['x_bins2'].value_counts()
    dfPoints_ant2_count.sort_index()
    
    
    if hemisphere == '_left':
        aDLS_combined = dfPoints_ant0_count[0]
        mDLS_combined = dfPoints_ant1_count[0]
        pDLS_combined = dfPoints_ant2_count[0]
        aDMS_combined = dfPoints_ant0_count[1]
        mDMS_combined = dfPoints_ant1_count[1]
        pDMS_combined = dfPoints_ant2_count[1]
    else:
        aDLS_combined = dfPoints_ant0_count[1]
        mDLS_combined = dfPoints_ant1_count[1]
        pDLS_combined = dfPoints_ant2_count[1]
        aDMS_combined = dfPoints_ant0_count[0]
        mDMS_combined = dfPoints_ant1_count[0]
        pDMS_combined = dfPoints_ant2_count[0]
    #    
    mouse_data = [sampleName, hemisphere, aDLS_combined, mDLS_combined, pDLS_combined, aDMS_combined, mDMS_combined, pDMS_combined]

    dataList.append(mouse_data)

allData = pd.DataFrame(data=dataList, columns=['mouse', 'hemisphere', 'aDLS', 'mDLS', 'pDLS', 'aDMS', 'mDMS', 'pDMS'])    
allData.to_excel('/d2/studies/ClearMap/IA_iDISCO/Striatum_Subregion_Counts' + hemisphere + '_3bins_Jan11.xlsx')


