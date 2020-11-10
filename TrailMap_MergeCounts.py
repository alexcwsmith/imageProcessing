#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 00:31:58 2019

@author: smith
"""

import pandas as pd
import os
import numpy as np
from scipy.stats import ttest_ind

#Set i/o directory
baseDirectory = '/d2/studies/ClearMap/FosTRAP_ChR2/batchProcessing/Counts_full/'

#Get Allen Brain region statistics
regionIndex = pd.read_csv(os.path.join(baseDirectory, 'RegionID_Index.csv'), index_col=0, encoding='ISO-8859-1')
regionVolumes = pd.read_csv(os.path.join(baseDirectory, 'RegionVolumeIndex.csv'), index_col=0, encoding='ISO-8859-1')

#Load data for each group
samples_veh = ['F1_NP', 'F1_RB', 'F2_NP', 'F2_LB']
samples_exp = ['F1_RT', 'F2_RT', 'F2_RB', 'F2_LT']

sample = 'F1_NP'
df = pd.read_csv(os.path.join(baseDirectory, 'Annotated_counts_' + sample + '_TrailMap_Reg3.csv'), header=None)
if df.shape[1] == 5:
    df.columns=['id', 'count_' + sample, 'region_' + sample, 'subregion_' + sample, 'subregion2_' + sample]
elif df.shape[1] == 6:
    df.columns=['id', 'count_' + sample, 'region_' + sample, 'subregion_' + sample, 'subregion2_' + sample, 'subregion3_' + sample]

#Merge & save data
idx = df['id'].tolist()
merged = pd.DataFrame(idx)
merged.columns=['id']
for sample in samples_exp:
    df = pd.read_csv(os.path.join(baseDirectory, 'Annotated_counts_' + sample + '_TrailMap_Reg3.csv'), header=None)
    if df.shape[1] == 5:
        df.columns=['id', 'count_' + sample, 'region_' + sample, 'subregion_' + sample, 'subregion2_' + sample]
    elif df.shape[1] == 6:
        df.columns=['id', 'count_' + sample, 'region_' + sample, 'subregion_' + sample, 'subregion2_' + sample, 'subregion3_' + sample]
    elif df.shape[1] ==4:
        df.columns=['id', 'count_' + sample, 'region_' + sample, 'subregion_' + sample]
    merged = merged.merge(df, left_on=['id'], right_on=['id'], how='outer')
    
merged.to_csv(os.path.join(baseDirectory, 'Annotated_Counts_Merged_FosTRAP.csv'))

df_exp = pd.read_csv(os.path.join(baseDirectory, 'Annotated_Counts_Merged_FosTRAP.csv'), index_col=0)
df_veh = pd.read_csv(os.path.join(baseDirectory, 'Annotated_Counts_Merged_Control.csv'), index_col=0)

#Compare groups
df_exp = df_exp.set_index('id', drop=True)
df_veh = df_veh.set_index('id', drop=True)
exp_idx = df_exp.index.tolist()
veh_idx = df_veh.index.tolist()

counts_exp = df_exp.iloc[:,[0,5,9,14]]
counts_exp['Region'] = regionIndex['name']
counts_exp['Mean_FosTRAP'] = np.mean(counts_exp, axis=1)
counts_exp['Stdev_FosTRAP'] = np.std(counts_exp, axis=1)
counts_exp['RegionVolume'] = regionVolumes['Volume']
counts_exp = counts_exp.dropna(axis=0, subset=['Region'])
counts_exp['Mean_VolumeNormalized'] = counts_exp['Mean_FosTRAP'] / counts_exp['RegionVolume']
counts_exp = counts_exp.sort_values(by='Mean_VolumeNormalized', ascending=False)
counts_exp.to_csv(os.path.join(baseDirectory, 'FosTRAP_Counts_Averages.csv'))

counts_veh = df_veh.iloc[:,[0,5,10,15]]
counts_veh['Region'] = regionIndex['name']
counts_veh['Mean_Control'] = np.mean(counts_veh,axis=1)
counts_veh['Mean_Control'] = np.std(counts_veh,axis=1)
counts_veh['RegionVolume'] = regionVolumes['Volume']
counts_veh = counts_veh.dropna(axis=0, subset=['Region'])
counts_veh['Mean_VolumeNormalized'] = counts_veh['Mean_Control'] / counts_veh['RegionVolume']
counts_veh = counts_veh.sort_values(by='Mean_VolumeNormalized', ascending=False)
counts_veh.to_csv(os.path.join(baseDirectory, 'Control_Counts_Averages.csv'))

comb_counts = pd.concat([counts_exp.iloc[:,0:4], counts_veh.iloc[:,0:4]], axis=1)
comb_counts = comb_counts.fillna(value=0)
pvals = np.array([ttest_ind(comb_counts.iloc[i, 0:4], comb_counts.iloc[i, 5:9])[1] for i in range(comb_counts.shape[0])])
pvalsDf = pd.DataFrame(pvals)
pvalsDf['id']=comb_counts.index
pvalsDf = pvalsDf.set_index('id', drop=True)
pvalsDf['Region']= regionIndex['name']
pvalsDf.columns=['pval', 'Region']
pvalsDf = pvalsDf.sort_values(by='pval', ascending=True)

pvalsDf.to_csv(os.path.join(baseDirectory, 'TrailMap_p-values.csv'))

comb_counts_vol = comb_counts.divide(regionVolumes['Volume'], axis=0)
comb_counts_vol['name'] = regionIndex['name']
comb_counts_vol = comb_counts_vol.dropna(how='any')
comb_counts_vol.to_csv(os.path.join(baseDirectory, 'Combined_Counts_VolumeNormalized.csv'))

comb_counts['name'] = regionIndex['name']
comb_counts['Volume'] = regionVolumes['Volume']
comb_counts.to_csv(os.path.join(baseDirectory, 'Combined_Counts.csv'))





