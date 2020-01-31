#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 11:55:35 2019

@author: smith
"""

import os
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_style('white')


baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/DeNardo_analysis/'
trap = pd.read_excel('/d2/studies/ClearMap/IA_iDISCO/analysisFiles/counts_table_processed.xlsx',header=0) #import whole brain cell counts
#behavior = pd.read_excel("Behavior.xlsx",header=0) #import corresponding behavioral analysis

from scipy.stats import zscore
# NOTE: Not normalizing by area volume
Alldata = np.array(trap.iloc[:,8:17])
Ydata = np.array(trap.iloc[:,8:12])
IAdata = np.array(trap.iloc[:,12:17])
dataWhich = trap.iloc[:,8:17]

All_Vol = zscore(np.array(trap.iloc[:,21:30]),axis=0)
Y_Vol = zscore(np.array(trap.iloc[:,21:25]),axis=0)
IA_Vol = zscore(np.array(trap.iloc[:,25:30]),axis=0)
dataWhich = trap.iloc[:,21:30]


Yoked_Zdata = zscore(np.array(trap.iloc[:,8:12]),axis=0) #z-score across rows
IA30_Zdata = zscore(np.array(trap.iloc[:,12:17]),axis=0) #z-score across rows

aDLS_idx = np.nonzero(trap['name']=='Anterior dorsolateral striatum')[0][0]
area_names = trap['name']
area_vols = trap['Volume']

import networkx as nx
from networkx.algorithms import community

def cluster_snn(data, max_pc=6, k=30):
    import igraph as ig
    from sklearn.neighbors import NearestNeighbors
    
    #select the algorithm that makes the most sense for your data
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='kd_tree').fit(data)
    neighbor_graph = nbrs.kneighbors_graph(data)
    g = ig.Graph()
    g = ig.GraphBase.Adjacency(neighbor_graph.toarray().tolist(), mode=ig.ADJ_UNDIRECTED)
    sim = np.array(g.similarity_jaccard())
    g = ig.GraphBase.Weighted_Adjacency(sim.tolist(), mode=ig.ADJ_UNDIRECTED)
    return np.array(g.community_multilevel(weights="weight", return_levels=False))

plots_dir = '/d2/studies/ClearMap/IA_iDISCO/DeNardo_analysis/' 
import os
from sklearn.manifold import TSNE



clu = cluster_snn(All_Vol)
clu = cluster_snn(np.corrcoef(all_fos))
#clu_z = cluster_snn(fosIA30_zscore)

cluster_areas = pd.DataFrame(data=np.array((area_names, clu)).T,columns=['Area', 'Cluster'])
#cluster_z-scores = pd.DataFrame(data=np.array((area_names, clu_z)).T, columns=['Area', 'Cluster'])
#export list of brain regions in each cluster:
cluster_areas.to_excel(os.path.join(baseDirectory, 'Combined_clusters_correlations.xlsx'))


ts = TSNE(perplexity=50,random_state=0).fit_transform(np.corrcoef(all_fos)) #select perplexity that makes the most sense for your data
plt.figure(figsize=(15,15), dpi=300)
plt.tick_params(axis='both', labelsize=24)
#can specify the color scheme by changing cmap (ex: cmap=plot.cm.Blues)
plt.scatter(ts[:,0],ts[:,1],s=100, c=clu,cmap=plt.cm.brg,lw=0,vmin=0,vmax=4)

#can enter a list of brain regions that you want to display on plots, must correspond to area_names
pl = np.argwhere(np.array(["Prelimbic" in i for i in area_names]))[0][0]
cpu = np.argwhere(np.array(["Caudoputamen" in i for i in area_names]))[0][0]
adls = np.argwhere(np.array(["Anterior dorsolateral striatum" in i for i in area_names]))[0][0]
pdms = np.argwhere(np.array(["Posterior dorsomedial striatum" in i for i in area_names]))[0][0]
bla = np.argwhere(np.array(["Basolateral amygdalar nucleus" in i for i in area_names]))[0][0]
bma = np.argwhere(np.array(["Medial amygdalar nucleus" in i for i in area_names]))[0][0]
nac = np.argwhere(np.array(["Nucleus accumbens" in i for i in area_names]))[0][0]
pld = np.argwhere(np.array(["Pallidum" in i for i in area_names]))[0][0]
ent = np.argwhere(np.array(["Entorhinal area" in i for i in area_names]))[0][0]
sma = np.argwhere(np.array(["Secondary motor area" in i for i in area_names]))[0][0]
pma = np.argwhere(np.array(["Primary motor area" in i for i in area_names]))[0][0]
hpc = np.argwhere(np.array(["Hippocampal formation" in i for i in area_names]))[0][0]
lha = np.argwhere(np.array(["Lateral hypothalamic area" in i for i in area_names]))[0][0]
cea = np.argwhere(np.array(["Central amygdalar nucleus" in i for i in area_names]))[0][0]
vta = np.argwhere(np.array(["Ventral tegmental area" in i for i in area_names]))[0][0]
sn = np.argwhere(np.array(["Substantia nigra" in i for i in area_names]))[0][0]
sub = np.argwhere(np.array(["Subiculum" in i for i in area_names]))[0][0]
presub = np.argwhere(np.array(["Presubiculum" in i for i in area_names]))[0][0]
postsub = np.argwhere(np.array(["Postsubiculum" in i for i in area_names]))[0][0]
rsc = np.argwhere(np.array(["Retrosplenial area" in i for i in area_names]))[0][0]
ca3 = np.argwhere(np.array(["Field CA3" in i for i in area_names]))[0][0]
acc = np.argwhere(np.array(["Anterior cingulate area" in i for i in area_names]))[0][0]
cpu = np.argwhere(np.array(["Caudoputamen-NOS" in i for i in area_names]))[0][0]
dg = np.argwhere(np.array(["Dentate gyrus polymorph layer" in i for i in area_names]))[0][0]

area_idx = [pl, cpu, adls, pdms, bla, bma, nac, 
            pld, ent, sma, pma, hpc, lha, cea, vta, sn, sub, presub, postsub, rsc, acc, cpu, dg]

names = ["PL", "CPU", "aDLS", "pDMS", 
         "BLA", "BMA", "NAc", "PD", "Ent", "SMA", "PMA", 
         "HPC", "LHA", "CeA", "VTA", "SN", "Sub", "PrS", "PoS", "RSC", "ACC", "CPu", "DG"]

for i, a in enumerate(area_idx):
    plt.text(ts[a,0],ts[a,1],names[i])
    matplotlib.rcParams.update({'font.size': 22})
    plt.scatter(ts[:,0],ts[:,1],s=50, c=clu,cmap=plt.cm.brg,lw=0,vmin=0,vmax=4)
plt.xlabel('t-SNE 1', fontsize=36)
plt.ylabel('t-SNE 2', fontsize=36)
plt.yticks(ticks=(np.arange(-20,21, step=5)))
plt.xticks(ticks=(np.arange(-20,21, step=5)))

#give your plot a unique name
plt.savefig(os.path.join(baseDirectory, "Combined_clusters_50p_correlations.png"), bbox_inches='tight')

baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/DeNardo_analysis/'
fos = pd.read_excel(os.path.join(baseDirectory, 'counts_table_6_processed.xlsx'))

from scipy.stats import zscore

fosCombined = np.array(fos.iloc[:,8:17])/np.tile(area_vols,[9,1]).T
fosYoked = zscore(np.array(trap.iloc[:,8:12]),axis=0)/np.tile(area_vols,[4,1]).T
fosIA30 = zscore(np.array(trap.iloc[:,12:17]), axis=0)/np.tile(area_vols,[5,1]).T


All_Zvol = zscore(np.array(trap.iloc[:,21:30]),axis=0)
Y_Zvol = zscore(np.array(trap.iloc[:,21:25]),axis=0)
IA_Zvol = zscore(np.array(trap.iloc[:,25:30]),axis=0)
dataWhich = trap.iloc[:,21:30]




from sklearn.decomposition import PCA
#fos_pc = PCA().fit(fos.iloc[:,21:30].T)
#fos_pc_proj = fos_pc.transform(fos.iloc[:,21:30].T)
fos_pc = PCA().fit(All_Vol)
fos_pc_proj = fos_pc.transform(IA_Vol)

fos_pc_projDf = pd.DataFrame(fos_pc_proj)
fos_pc_projDf.to_excel(os.path.join(baseDirectory, 'Yoked_fos_pc_proj.xlsx'))

clu = cluster_snn(fos_pc_proj)

fos_pc_proj = np.array(fos_pc_proj)
fos_pcs = pd.DataFrame(fos_pc_proj, columns=['PCProjection']) #, 'PCProjection'
cluData = pd.DataFrame(np.array((area_names, clu)).T)
cluster_areas = pd.DataFrame(data, columns=['Area', 'Cluster'])

cluData.to_excel(os.path.join(baseDirectory, 'Yoked_PCA_clusters_norm.xlsx')) 

import os
plt.figure(figsize=(5,5))
plt.axvline(0,color='k', linestyle='--')
plt.axhline(0,color='k',linestyle='--')
#plt.plot(fos_pc_proj[:7,0],fos_pc_proj[:7,1],'bo')
#plt.plot(fos_pc_proj[7:,0],fos_pc_proj[7:,1],'o',color='orange')
plt.plot(fos_pc_proj[:4,0],fos_pc_proj[:4,1],'bo')
plt.plot(fos_pc_proj[4:9,0],fos_pc_proj[4:9,1],'o',color='orange')
#plt.plot(fos_pc_proj[11:16,0],fos_pc_proj[11:16,1],'o',color='green')
#plt.plot(fos_pc_proj[17:21,0],fos_pc_proj[17:21,1],'o',color='red')
plt.xlabel('PC 1 (%0.02f\\%% variance)' % (fos_pc.explained_variance_ratio_[0]*100))
plt.ylabel('PC 2 (%0.02f\\%% variance)' % (fos_pc.explained_variance_ratio_[1]*100))
plt.savefig(os.path.join(baseDirectory, "fos_pca_experiments.pdf"))

# plot PC contributions by brain region
plt.figure(figsize=(35,50))
plt.pcolor(fos_pc.components_[0:3].T,cmap=plt.cm.bwr,vmin=-0.6,vmax=0.6)
plt.axes().set_yticks(np.arange(fos_pc.components_.shape[1]))
plt.axes().set_yticklabels(fos['name'])
plt.ylim([0, fos_pc.components_.shape[1]])
plt.colorbar()
plt.savefig(os.path.join(plots_dir, "Combined_fos_pcs_heatmap_2.pdf"))

All_Zvol = zscore(np.array(trap.iloc[:,21:30]),axis=0)
Y_Zvol = zscore(np.array(trap.iloc[:,21:25]),axis=0)
IA_Zvol = zscore(np.array(trap.iloc[:,25:30]),axis=0)
dataWhich = trap.iloc[:,21:30]


fosYoked = np.array(fos.iloc[:,21:25]).astype(float)
fosIA30 = np.array(fos.iloc[:,25:30]).astype(float)
all_fos = np.array(fos.iloc[:,21:30]).astype(float)
all_fos_zscore = zscore(all_fos, axis=0)
fosYoked_zscore = zscore(fosYoked,axis=0)
fosIA30_zscore = zscore(fosIA30, axis=0)
allfos_zscores = pd.DataFrame(all_fos_zscore)
allfos_zscores.to_excel(os.path.join(baseDirectory, 'zscores_all_volumeNormalized.xlsx'))
fosIA30_zscores=pd.DataFrame(fosIA30_zscore)
fosYoked_zscores=pd.DataFrame(fosYoked_zscore)
fosIA30_zscores.to_excel('IA30_zscores_volumeNormalized.xlsx')
fosYoked_zscores.to_excel('Yoked_zscores_volumeNormalized.xlsx')

from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
expt_pc = PCA().fit(all_fos_zscore)
expt_pc_proj = expt_pc.transform(all_fos_zscore)

labels = np.zeros((398,9))
labels[:398] = 1
expt_lda = LinearDiscriminantAnalysis().fit(fosIA30_zscore, labels)
expt_lda_proj = expt_lda.transform(fosIA30_zscore)

plt.figure(figsize=(10,10))
coef = expt_lda.coef_
plt.ylabel(pd.DataFrame(area_names[np.argsort(coef)[0]), rotation=0, fontsize=14, labelpad=20)
plt.locator_params(axis='y', tight=True)
plt.imshow(np.sort(coef).T, aspect='auto', cmap=plt.cm.bwr)
plt.savefig(os.path.join(plots_dir, 'LDA_analysis.pdf'))
area_names[np.argsort(coef)[0]]
LDA_areas = pd.DataFrame(area_names[np.argsort(coef)[0]])
LDA_areas.to_excel(os.path.join(baseDirectory, 'LDA_regions.xlsx'))

plt.imshow(np.sort(coef).T,aspect='auto', interpolation='none',cmap=plt.cm.bwr)
plt.savefig(os.path.join(plots_dir, 'LDA_analysis.pdf'))
area_names[np.argsort(coef)[0]]
LDA_areas = pd.DataFrame(area_names[np.argsort(coef)[0]])



plt.plot(expt_pc.explained_variance_ratio_*100,'ko-')
plt.xlabel('PC number')
plt.ylabel('% variance explained')
plt.savefig(os.path.join(plots_dir, "combined_pca_variance_explained.pdf"))

plt.figure(figsize=(5,5))
plt.axvline(0,color='k', linestyle='--')
plt.axhline(0,color='k',linestyle='--')
plt.plot(expt_pc_proj[:4,0],expt_pc_proj[:4,1],'bo')
plt.plot(expt_pc_proj[4:,0],expt_pc_proj[4:,1],'o',color='orange')
plt.xlabel('PC 1 (%0.02f\\%% variance)' % (expt_pc.explained_variance_ratio_[0]*100))
plt.ylabel('PC 2 (%0.02f\\%% variance)' % (expt_pc.explained_variance_ratio_[1]*100))
plt.savefig(os.path.join(plots_dir, "combined_pca_experiments.pdf"))

# li style plot
plt.figure(figsize=(30,50))
plt.pcolor(expt_pc.components_[0:200].T,cmap=plt.cm.bwr,vmin=-0.6,vmax=0.6)
plt.axes().set_yticks(np.arange(expt_pc.components_.shape[1]))
plt.axes().set_yticklabels(fos['name'])
plt.ylim([0, expt_pc.components_.shape[1]])
plt.colorbar()
plt.savefig(os.path.join(plots_dir, "pcs_heatmap.pdf"))

PCs = pd.DataFrame(data=expt_pc.components_[:3,:].T, index=fos['name'],columns=[1,2,3])
PCs.to_excel(os.path.join(plots_dir, "pcs.xlsx"))

from scipy.stats import ranksums, ttest_ind
from statsmodels.stats.multitest import multipletests

fosYoked = np.array(fos.iloc[:,21:25]).astype(float)
fosIA30 = np.array(fos.iloc[:,25:30]).astype(float)
all_fos = np.array(fos.iloc[:,21:30]).astype(float)
all_fos_zscore = zscore(all_fos, axis=0)
fosYoked_zscore = zscore(fosYoked,axis=0)
fosIA30_zscore = zscore(fosIA30, axis=0)


pvals = np.array([ttest_ind(fosYoked_zscore[i,:], fosIA30_zscore[i,:])[1] for i in range(fosYoked_zscore.shape[0])])
fdr_pass,qvals,_,_ = multipletests(pvals, method='fdr_bh',alpha=0.05)
volPdf = pd.DataFrame(pvals)
volPdf['Area'] = pd.DataFrame(area_names)
volPdf.to_excel(os.path.join(baseDirectory, 'pvalues_Zscored_VolumeNormalized.xlsx'))

pvals = np.array([ttest_ind(fosYoked[i,:], fosIA30[i,:])[1] for i in range(Yoked_Zdata.shape[0])])
fdr_pass,qvals,_,_ = multipletests(pvals, method='fdr_bh',alpha=0.05)
volPdf = pd.DataFrame(pvals)
volPdf['Area'] = pd.DataFrame(area_names)
volPdf['Latency_Correlation_Coefficient']=latency_cc
volPdf.to_excel(os.path.join(baseDirectory, 'pvalues_volumeNormalized.xlsx'))

area_names[pvals<0.05]
areaDf = pd.DataFrame(area_names[pvals<.05])
areaDf.to_excel(os.path.join(baseDirectory, 'significant_regions.xlsx'))

pvals = np.array([ttest_ind(Yoked_Zdata[i,:], IA30_Zdata[i,:])[1] for i in range(Y_Zvol.shape[0])])
fdr_pass,qvals,_,_ = multipletests(pvals, method='fdr_bh',alpha=0.05)

meanYoked = np.mean(fosYoked,axis=1) 
meanIA30 = np.mean(fosIA30,axis=1)
stdYoked = np.std(fosYoked,axis=1)
stdIA30 = np.std(fosIA30,axis=1)
all_comparisons = pd.DataFrame(np.array([pvals, qvals, meanYoked, stdYoked, meanIA30, stdIA30]).T,columns=["PVal", "QVal", "MeanYoked", "StdYoked", "MeanIA30", "StdIA30"],index=area_names)
all_comparisons.to_excel(os.path.join(baseDirectory, "fos_comparisons_descriptiveStatistics.xlsx"))

pvals = np.array([ttest_ind(Y_Zvol[i,:], IA_Zvol[i,:])[1] for i in range(Y_Zvol.shape[0])])
fdr_pass,qvals,_,_ = multipletests(pvals, method='fdr_bh',alpha=0.05)

meanYokedVol = np.mean(Y_Zvol, axis=1) 
meanIA30Vol = np.mean(IA_Zvol, axis=1)
stdYokedVol = np.std(Y_Zvol, axis=1)
stdIA30Vol = np.std(IA_Zvol, axis=1)
all_comparisons = pd.DataFrame(np.array([pvals, qvals, meanYokedVol, stdYokedVol, meanIA30Vol, stdIA30Vol]).T,columns=["PVal", "QVal", "MeanYoked", "StdYoked", "MeanIA30", "StdIA30"],index=area_names)
all_comparisons.to_excel(os.path.join(baseDirectory, "zscored_fos_comparisons_volumeNormalized.xlsx"))

cluYoked = cluster_snn(fosYoked_zscore)
cluIA30 = cluster_snn(fosIA30_zscore)

cluIA30df = pd.DataFrame(cluIA30)
cluIA30df.to_excel('IA30_zscores_clusters.xlsx')
cluYokeddf = pd.DataFrame(cluYoked)
cluYokeddf.to_excel('Yoked_zscores_clusters.xlsx')

cluIA30.astype(np.int)

area_names[cluIA30==0]
plt.figure(figsize=(10,10))
plt.imshow(area_names())



plt.figure(figsize=(25,25))
plt.imshow(np.corrcoef(fosYoked[:16]),aspect='auto', interpolation='none',cmap=plt.cm.bwr,vmin=-1,vmax=1)
matplotlib.rcParams.update({'font.size': 22})
#plt.axes().set_yticklabels(fos['name'][:])
plt.yticks(np.arange(16), (fos['name']), fontsize=28)
plt.xticks(np.arange(16), (fos['name']), fontsize=28, rotation='vertical')
#plt.axes().set_xticklabels(fos['name'][:])
plt.colorbar()
plt.savefig(os.path.join(baseDirectory, "Yoked_correlations_volumeNormalized_top16.png"), bbox_inches='tight')

YokedCorrelations = pd.DataFrame(np.corrcoef(fosYoked))
YokedCorrelations.to_excel(os.path.join(baseDirectory, 'Yoked_correlations.xlsx'))


plt.figure(figsize=(25,25))
plt.imshow(np.corrcoef(fosIA30[:16]),aspect='auto', interpolation='none',cmap=plt.cm.bwr,vmin=-1,vmax=1)
matplotlib.rcParams.update({'font.size': 22})
#plt.axes().set_yticklabels(fos['name'][:], fontsize=12)
plt.yticks(np.arange(16), (fos['name']), fontsize=28)
plt.xticks(np.arange(16), (fos['name']), fontsize=28, rotation='vertical')
#plt.axes().set_xticklabels(fos['name'][:])
plt.colorbar()
plt.savefig(os.path.join(baseDirectory, "IA30_correlations_volumeNormalized_top16.png"), bbox_inches='tight')

IA30_Correlations = pd.DataFrame(np.corrcoef(fosIA30))
IA30_Correlations.to_excel(os.path.join(baseDirectory, 'IA30_correlations.xlsx'))

plt.figure(figsize=(120,120))
plt.imshow((np.corrcoef(fosIA30[:398]) - np.corrcoef(fosYoked[:398])),aspect='auto', interpolation='none',cmap=plt.cm.bwr,vmin=-1,vmax=1)
matplotlib.rcParams.update({'font.size': 22})
#plt.axes().set_yticklabels(fos['name'][:], fontsize=12)
plt.yticks(np.arange(398), (fos['name']), fontsize=18)
plt.xticks(np.arange(398), (fos['name']), fontsize=18, rotation='vertical')
#plt.axes().set_xticklabels(fos['name'][:])
plt.colorbar()
plt.savefig(os.path.join(baseDirectory, "Delta_correlations_volumeNormalized_full.png"), bbox_inches='tight')


##CORRELATE WITH BEHAVIOR:
behavior = pd.read_excel('behavior.xlsx')
latency = behavior['Latency'].T


from sklearn.linear_model import LinearRegression
nareas = fosIA30.shape[0]
latency_cc = np.zeros((nareas,))
for i in range(nareas):
    latency_cc[i] = np.corrcoef(fosIA30[i,:], latency)[0,1]

corr_coef = pd.DataFrame(data=np.array((area_names, latency_cc)).T,columns=['Area', 'Latency'])
corr_coef['Yoked_vs_IA_pvalue']=pvals
corr_coef.to_excel('latency_cc.xlsx')

