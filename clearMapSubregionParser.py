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


sampleName = 'IA2_LB'
baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/' + sampleName

##IMPORT PREVIOUSLY PRE-PROCESSED HeatMap Data
hemisphere = '_Right'
heatmap = io.readData(os.path.join(baseDirectory, sampleName + '_Caudoputamen' + hemisphere + '.tif'))
data = io.readData(heatmap)
points = np.nonzero(data)
dfPoints = pd.DataFrame(points, index=['x', 'y', 'z']).T
dfPoints.rename(columns={0: "x", 1: "y", 2: "z"})

###View Range:
x_range=dfPoints.x.max() - dfPoints.x.min()
y_range=dfPoints.y.max() - dfPoints.y.min()
print(x_range, y_range)

#Bin Y axis
dfPoints['y_bins']=pd.cut(dfPoints['y'], bins=2)
dfPoints['y_bins'].value_counts()
dfPoints_counts_y = dfPoints['y_bins'].value_counts()
dfPoints_sorted_y = dfPoints_counts_y.sort_index()
dfPoints_sorted_y.to_excel(os.path.join(baseDirectory, 'y_binned_striatum_2' + hemisphere + '.xlsx'))

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

#Merge data sets & write to excel:
ant01 = pd.merge(dfPoints_ant0, dfPoints_ant1, how='outer')
ant012 = pd.merge(ant01, dfPoints_ant2, how='outer')

if hemisphere == '_Right':
    dfPoints_ant_compiled = pd.DataFrame([ant012['x_bins2'].value_counts(),ant012['x_bins1'].value_counts(),ant012['x_bins0'].value_counts()])
    dfPoints_ant_compiled.columns=['Medial', 'Lateral', 'Bin', 'Bin2']
else:
    dfPoints_ant_compiled = pd.DataFrame([ant012['x_bins0'].value_counts(),ant012['x_bins1'].value_counts(),ant012['x_bins2'].value_counts()])
    dfPoints_ant_compiled.columns=['Bin', 'Lateral', 'Medial']
    
dfPoints_ant_compiled.sort_index()

dfPoints_ant_compiled.to_excel(os.path.join(baseDirectory, 'Anterior_CellCounts' + hemisphere + '.xlsx'))


