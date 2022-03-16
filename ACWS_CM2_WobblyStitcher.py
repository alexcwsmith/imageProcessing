#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 23:42:44 2020

@author: smith
"""
import os
os.chdir('/d1/software/ClearMap2/')
%gui qt5
from ClearMap.Environment import *  #analysis:ignore

directory = '/d1/studies/ClearMap/MOR-mCherry_4x/FMC1_ECiTest/'    
#_<Y,2>x<X,2>_Tile.ome.npy
expression_raw = 'FMC1_MOR_<Y,2>x<X,2>_Tile.ome.tif'

ws = wsp.Workspace('CellMap', directory=directory);
ws.update(raw=expression_raw)
ws.info()

io.convert_files(ws.file_list('rawlip', extension='tif'), extension='npy',
                 processes=22, verbose=True);

#Rigid Stiching
layout = stw.WobblyLayout(expression=ws.filename('raw'), tile_axes=['X', 'Y'], overlaps=(228,270));

st.align_layout_rigid_mip(layout, depth=[260, 300, None], max_shifts=[(-75,75),(-30,30),(-10,10)],
                          ranges = [None,None,None], background=(45, 200), clip=7500,
                           processes=16, verbose=True)

st.place_layout(layout, method='optimization', min_quality=-np.inf, lower_to_origin=True, verbose=True)
st.save_layout(ws.filename('layout', postfix='aligned_axis'), layout)

#Wobbly Alignment
##RESUME HERE 9/9
layout = st.load_layout(ws.filename('layout', postfix='aligned_axis'))

stw.align_layout(layout, axis_range=(None, None, 3), max_shifts=[(-75,75),(-30,30),(0,0)], axis_mip=None,
                 validate=dict(method='foreground', valid_range=(130, None), size=None),
                 prepare =dict(method='normalization', clip=None, normalize=True),
                 validate_slice=dict(method='foreground', valid_range=(45,5000), size=500),
                 prepare_slice =None,
                 find_shifts=dict(method='tracing', cutoff=3*np.sqrt(2)),
                 processes=10, verbose=True)

st.save_layout(ws.filename('layout', postfix='aligned'), layout)

layout = st.load_layout(ws.filename('layout', postfix='aligned'));

stw.place_layout(layout, min_quality = -np.inf,
                 method = 'optimization',
                 smooth = dict(method = 'window', window = 'bartlett', window_length = 100, binary = None),
                 smooth_optimized = dict(method = 'window', window = 'bartlett', window_length = 20, binary = 10),
                 fix_isolated = False, lower_to_origin = True,
                 processes = 10, verbose = True)

st.save_layout(ws.filename('layout', postfix='placed'), layout)


###Stitching
layout = st.load_layout(ws.filename('layout', postfix='placed'));

stw.stitch_layout(layout, sink = ws.filename('stitched'),
                  method = 'interpolation', processes=18, verbose=True)


resample_parameter = {
    "source_resolution" : (1.625,1.625,6.5),
    "sink_resolution"   : (6.5,6.5,6.5),
    "processes" : 12,
    "verbose" : True,
    };

res.resample(ws.filename('stitched'), sink=ws.filename('resampled'), **resample_parameter)


source = ws.source('stitched')
sink = ws.filename('stitched', extension='tif')

io.convert(source, sink, processes=18, verbose=True)

import os
from skimage import io
from skimage import restoration
from scipy.ndimage.filters import gaussian_filter
from skimage.transform import resize
directory = '/d1/studies/ClearMap/MOR-mCherry_4x/FMC1_ECiTest/'

im = io.imread(os.path.join(directory, 'FMC1_MOR_mCherry_MinSubtraction_GaussianSigma6e-1.tif'))

arrays = []
for plane in range(im.shape[0]):
    p = im[plane]
    m = im[plane].min()
    plane_corr = p-m
    arrays.append(plane_corr)
    
im_min = io.concatenate_images(arrays)
io.imsave(os.path.join(directory, 'stitched_forebrain_minSub.tif'), im_min)

im_gaus = gaussian_filter(im_min, sigma=0.6, mode='wrap')

im_gaus_med = median_filter(im_gaus, size=(11,11,3), mode='wrap')

im_resize = resize(im, output_shape=(1009,828,1736), mode='wrap')


