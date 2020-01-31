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
baseDirectory = '/d2/studies/ClearMap/IA_iDISCO/analysisFiles'



# Voxel-based statistics:
#########################

#Load the data (heat maps generated previously )
group1 = [os.path.join(baseDirectory, 'cells_heatmap_vox15_IA1_LB.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_IA2_NP.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_IA2_LT.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_IA2_LB.tif')
          ]

group2 = [os.path.join(baseDirectory, 'cells_heatmap_vox15_IA1_RT.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_IA1_RB.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_IA1_LT.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_IA2_RT.tif'),
          os.path.join(baseDirectory, 'cells_heatmap_vox15_IA2_RB.tif')
          ]
 

g1 = stat.readDataGroup(group1);
g2 = stat.readDataGroup(group2);



#Generated average and standard deviation maps
##############################################
g1a = numpy.mean(g1,axis = 0);
g1s = numpy.std(g1,axis = 0);

g2a = numpy.mean(g2,axis = 0);
g2s = numpy.std(g2,axis = 0);

io.writeData(os.path.join(baseDirectory, 'Yoked_full_mean_horizontal_vox15.mhd'), g1a)# rsp.sagittalToCoronalData(g1a));
io.writeData(os.path.join(baseDirectory, 'Yoked_full_std_horizontal_vox15.mhd'), g1s)# rsp.sagittalToCoronalData(g1s));

io.writeData(os.path.join(baseDirectory, 'IA30_full_mean_horizontal_vox15.mhd'), g2a)# rsp.sagittalToCoronalData(g2a));
io.writeData(os.path.join(baseDirectory, 'IA30_full_std_horizontal_vox15.mhd'), g2s)# rsp.sagittalToCoronalData(g2s));





#Generate the p-values map
##########################
#pcutoff: only display pixels below this level of significance
pvals, psign = stat.tTestVoxelization(g1.astype('float'), g2.astype('float'), signed = True, pcutoff = 0.05);

#color the p-values according to their sign (defined by the sign of the difference of the means between the 2 groups)
pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
io.writeData(os.path.join(baseDirectory, 'pvalues_full_vox15.tif'), pvalsc.astype('float32'))#, rsp.sagittalToCoronalData(pvalsc.astype('float32')));



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



#ids, pc1, pc1i = stat.countPointsGroupInRegions(group1, intensityGroup = group1i, withIds = True, labeledImage = lbl.DefaultLabeledImageFile, withCounts = True);
#pc2, pc2i = stat.countPointsGroupInRegions(group2, intensityGroup = group2i, withIds = False, labeledImage = lbl.DefaultLabeledImageFile, withCounts = True);

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


#sort by qvalue
ii = numpy.argsort(pvalsi0);
tableSorted = table.copy();
tableSorted = tableSorted[ii];

with open(os.path.join(baseDirectory, 'counts-intensity_table_6.csv'),'w') as f:
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


#sort by qvalue
ii = numpy.argsort(pvals0);
tableSorted = table.copy();
tableSorted = tableSorted[ii];

with open(os.path.join(baseDirectory, 'counts_table_6.csv'),'w') as f:
    f.write(', '.join([str(item) for item in table.dtype.names]));
    f.write('\n');
    for sublist in tableSorted:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


#############################################################################


baseDirectory = '/home/mtllab/Documents/bicuculline'


group2 = [#'/home/mtllab/Documents/bicuculline/2/cells_heatmap.tif',
       #   '/home/mtllab/Documents/bicuculline/3/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/bicuculline/7/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/bicuculline/6/cells_transformed_to_Atlas.npy',
       #   '/home/mtllab/Documents/bicuculline/8/cells_transformed_to_Atlas.npy',
        #  '/home/mtllab/Documents/bicuculline/9/cells_transformed_to_Atlas.npy',
       #   '/home/mtllab/Documents/bicuculline/13/cells_transformed_to_Atlas.npy']
          '/home/mtllab/Documents/bicuculline/14/cells_transformed_to_Atlas.npy']
          
group1 = ['/home/mtllab/Documents/bicuculline/10/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/bicuculline/11/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/bicuculline/12/cells_transformed_to_Atlas.npy']



group1i = [fn.replace('cells_transformed_to_Atlas', 'intensities') for fn in group1];
group2i = [fn.replace('cells_transformed_to_Atlas', 'intensities') for fn in group2];




AnnotationFile = os.path.join(PathReg, 'annotation_25_right_fullWD.tif');



#ids, pc1, pc1i = stat.countPointsGroupInRegions(group1, intensityGroup = group1i, withIds = True, labeledImage = lbl.DefaultLabeledImageFile, withCounts = True);
#pc2, pc2i = stat.countPointsGroupInRegions(group2, intensityGroup = group2i, withIds = False, labeledImage = lbl.DefaultLabeledImageFile, withCounts = True);

ids, pc1, pc1i = stat.countPointsGroupInRegions(group1, intensityGroup = group1i, returnIds = True, labeledImage = AnnotationFile, returnCounts = True, collapse=True);
pc2, pc2i = stat.countPointsGroupInRegions(group2, intensityGroup = group2i, returnIds = False, labeledImage = AnnotationFile, returnCounts = True, collapse=True);


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


#sort by qvalue
ii = numpy.argsort(pvalsi0);
tableSorted = table.copy();
tableSorted = tableSorted[ii];

with open(os.path.join(baseDirectory, 'counts-intensity_table_slow.csv'),'w') as f:
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


#sort by qvalue
ii = numpy.argsort(pvals0);
tableSorted = table.copy();
tableSorted = tableSorted[ii];

with open(os.path.join(baseDirectory, 'counts_table_slow.csv','w') as f:
    f.write(', '.join([str(item) for item in table.dtype.names]));
    f.write('\n');
    for sublist in tableSorted:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();