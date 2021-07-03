# -*- coding: utf-8 -*-
"""
Example script to set up groupwise comparison of ClearMap data.

This script is lightly modified from the ClearMap package written by Christoph Kirst,
full repository found at github.com/ChristophKirst/ClearMap
"""
# import modules:

import ClearMap.Analysis.Statistics as stat
import ClearMap.Analysis.Label as lbl
import ClearMap.IO.IO as io
import ClearMap.Alignment.Resampling as rsp
import numpy, os


# Base Directory, usually where your experiment is saved:
baseDirectory = '/d2/studies/ClearMap/ROC_iDISCO/batchDirectory'
sampleName='ROC_iDISCO'


# Voxel-based statistics:
#########################

#Load the data (heat maps generated previously )
group1 = [os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC9.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC10.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC11.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC12.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC13.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC14.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC15.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC16.tif')
          ]

group2 = [os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC1.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC2.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC4.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC6.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC7.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC8.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC17.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_ROC18.tif'),
          
          ]
 

g1 = stat.readDataGroup(group1);
g2 = stat.readDataGroup(group2);



#Generated average and standard deviation maps
##############################################
g1a = numpy.mean(g1,axis = 0);
g1s = numpy.std(g1,axis = 0);

g2a = numpy.mean(g2,axis = 0);
g2s = numpy.std(g2,axis = 0);

io.writeData(os.path.join(baseDirectory, 'Group1_full_mean_horizontal.mhd'), g1a)# rsp.sagittalToCoronalData(g1a));
io.writeData(os.path.join(baseDirectory, 'Group1_full_std_horizontal.mhd'), g1s)# rsp.sagittalToCoronalData(g1s));

io.writeData(os.path.join(baseDirectory, 'Group2_full_mean_horizontal.mhd'), g2a)# rsp.sagittalToCoronalData(g2a));
io.writeData(os.path.join(baseDirectory, 'Group2_full_std_horizontal.mhd'), g2s)# rsp.sagittalToCoronalData(g2s));


#Generate the p-values map
##########################
#pcutoff: only display pixels below this level of significance
pvals, psign = stat.tTestVoxelization(g1.astype('float'), g2.astype('float'), signed = True, pcutoff = 0.05);

#color the p-values according to their sign (defined by the sign of the difference of the means between the 2 groups)
pvalsc = stat.colorPValues(pvals, psign, positive = [1,0], negative = [0,1]);
io.writeData(os.path.join(baseDirectory, sampleName + '_pvalues.tif'), pvalsc.astype('float32'))#, rsp.sagittalToCoronalData(pvalsc.astype('float32')));



#############################################################################
# Regions-based statistics:
###########################



group1 = ['/d2/studies/ClearMap/IA_iDISCO/IA1_LB/cells_transformed_to_Atlas.npy',
          '/d2/studies/ClearMap/IA_iDISCO/IA2_NP/cells_transformed_to_Atlas.npy',
          '/d2/studies/ClearMap/IA_iDISCO/IA2_LT/cells_transformed_to_Atlas.npy',
          '/d2/studies/ClearMap/IA_iDISCO/IA2_LB/cells_transformed_to_Atlas.npy'
          ]

        
group2 = ['/d2/studies/ClearMap/IA_iDISCO/IA1_RT/cells_transformed_to_Atlas.npy',
          '/d2/studies/ClearMap/IA_iDISCO/IA1_RB/cells_transformed_to_Atlas.npy',
          '/d2/studies/ClearMap/IA_iDISCO/IA1_LT/cells_transformed_to_Atlas.npy',
          '/d2/studies/ClearMap/IA_iDISCO/IA2_RT/cells_transformed_to_Atlas.npy',
          '/d2/studies/ClearMap/IA_iDISCO/IA2_RB/cells_transformed_to_Atlas.npy'
          ]

group1i = [fn.replace('cells_transformed_to_Atlas', 'intensities') for fn in group1];
group2i = [fn.replace('cells_transformed_to_Atlas', 'intensities') for fn in group2];



PathReg        = '/d1/studies/ClearMap/ClearMap_resources/Parameter_files/';
AnnotationFile = '/d2/studies/ClearMap/IA_iDISCO/Annotations_Horizontal_8-302_Full.tif';

#Calculate stats for regions. Set collapse=True to use collapsed regions (e.g. collapse cortical layers into single region)
ids, pc1, pc1i = stat.countPointsGroupInRegions(group1, intensityGroup = group1i, returnIds = True, labeledImage = AnnotationFile, returnCounts = True, collapse=None);
pc2, pc2i = stat.countPointsGroupInRegions(group2, intensityGroup = group2i, returnIds = False, labeledImage = AnnotationFile, returnCounts = True, collapse=None);


pvals, psign = stat.tTestPointsInRegions(pc1, pc2, pcutoff = None, signed = True);
pvalsi, psigni = stat.tTestPointsInRegions(pc1i, pc2i, pcutoff = None, signed = True, equal_var = True);

import ClearMap.Analysis.Tools.MultipleComparisonCorrection as FDR


iid = pvalsi < 1;

ids0 = ids[iid];
pc1i0 = pc1i[iid];
pc2i0 = pc2i[iid];
pc10 = pc1[iid];
pc20 = pc2[iid];
psigni0 = psigni[iid];
pvalsi0 = pvalsi[iid];
qvalsi0 = FDR.estimateQValues(pvalsi0);
psign0 = psign[iid];
pvals0 = pvals[iid];
qvals0 = FDR.estimateQValues(pvals0);


#make table

dtypes = [('id','int64'),('mean1','f8'),('std1','f8'),('mean2','f8'),('std2','f8'),('pvalue', 'f8'),('qvalue', 'f8'),('psign', 'int64')];
for i in range(len(group1)):
    dtypes.append(('count1_%d' % i, 'f8'));
for i in range(len(group2)):
    dtypes.append(('count2_%d' % i, 'f8'));   
dtypes.append(('name', 'a256'));

table = numpy.zeros(ids0.shape, dtype = dtypes)
table["id"] = ids0;
table["mean1"] = pc1i0.mean(axis = 1)/1000000;
table["std1"] = pc1i0.std(axis = 1)/1000000;
table["mean2"] = pc2i0.mean(axis = 1)/1000000;
table["std2"] = pc2i0.std(axis = 1)/1000000;
table["pvalue"] = pvalsi0;
table["qvalue"] = qvalsi0;

table["psign"] = psigni0;
for i in range(len(group1)):
    table["count1_%d" % i] = pc10[:,i];
for i in range(len(group2)):
    table["count2_%d" % i] = pc20[:,i];
table["name"] = lbl.labelToName(ids0);


#sort by qvalue & write intensity table
ii = numpy.argsort(pvalsi0);
tableSorted = table.copy();
tableSorted = tableSorted[ii];

with open(os.path.join(baseDirectory, 'counts-intensity_table.csv'),'w') as f:
    f.write(', '.join([str(item) for item in table.dtype.names]));
    f.write('\n');
    for sublist in tableSorted:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

#############################


#make table

dtypes = [('id','int64'),('mean1','f8'),('std1','f8'),('mean2','f8'),('std2','f8'),('pvalue', 'f8'),('qvalue', 'f8'),('psign', 'int64')];
for i in range(len(group1)):
    dtypes.append(('count1_%d' % i, 'f8'));
for i in range(len(group2)):
    dtypes.append(('count2_%d' % i, 'f8'));   
dtypes.append(('name', 'a256'));

table = numpy.zeros(ids0.shape, dtype = dtypes)
table["id"] = ids0;
table["mean1"] = pc10.mean(axis = 1);
table["std1"] = pc10.std(axis = 1);
table["mean2"] = pc20.mean(axis = 1);
table["std2"] = pc20.std(axis = 1);
table["pvalue"] = pvals0;
table["qvalue"] = qvals0;

table["psign"] = psigni0;
for i in range(len(group1)):
    table["count1_%d" % i] = pc10[:,i];
for i in range(len(group2)):
    table["count2_%d" % i] = pc20[:,i];
table["name"] = lbl.labelToName(ids0);


#sort by qvalue & write counts table
ii = numpy.argsort(pvals0);
tableSorted = table.copy();
tableSorted = tableSorted[ii];

with open(os.path.join(baseDirectory, 'counts_table.csv'),'w') as f:
    f.write(', '.join([str(item) for item in table.dtype.names]));
    f.write('\n');
    for sublist in tableSorted:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

