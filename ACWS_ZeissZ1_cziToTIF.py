#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 14:45:19 2022

@author: smith
"""

from aicsimageio.readers import CziReader
from aicsimageio.writers import OmeTiffWriter
import os
import numpy as np
import argparse
import time

tileMap4w = {0:'00x00', 1:'00x01', 2:'00x02', 3:'00x03', 4:'01x03', 5:'01x02',
           6:'01x01', 7:'01x00', 8:'02x00', 9:'02x01', 10:'02x02', 11:'02x03',
           12:'03x03', 13:'03x02', 14:'03x01', 15:'03x00', 16:'04x00', 17:'04x01',
           18:'04x02', 19:'04x03'}

tileMap3w = {0:'00x00', 1:'00x01', 2:'00x02', 3: '01x02', 4: '01x01', 5:'01x00',
            6: '02x00', 7:'02x01', 8:'02x02', 9:'03x02', 10:'03x01', 11:'03x00',
            12: '04x00', 13:'04x01', 14:'04x02', 15:'05x02', 16:'05x01', 17:'05x00'}

def parseTiles(imPath):
    """Tool for extracting TIF stacks for each tile from CZI data produced by a Zeiss Z1 light-sheet.
    
    Parameters
    ----------
    imPath : str
        Full path to .czi file to extract data from.
    """
    print("Started at " + time.ctime())
    directory = os.path.dirname(imPath)
    basename = os.path.basename(imPath)
    n, e = os.path.splitext(basename)
    im = CziReader(imPath)
    print("Loading image with dimensions " + str(im.shape))
    img = im.get_image_data()
    print(im.dims.order, im.shape)
    img = np.squeeze(img)
    if img.shape[0]>20:
        raise ValueError("Images with over 20 tiles not yet supported, add position mappings to dictionary")
    if img.shape[0] % 4!=0 and img.shape[0] %3!=0:
        raise ValueError("Image is not 3 or 4 tiles wide, does not fit position mapping dictionary")
    print(img.shape) #tiles, channels, z, y, x
    if img.shape[0]%3==0:
        tilemap=tileMap3w
    elif img.shape[0]%4==0:
        tilemap=tileMap4w
    for tile in range(img.shape[0]):
        for channel in range(img.shape[1]):
            image = img[tile, channel, :, :, :]
            
            OmeTiffWriter.save(data=image, 
                            uri=os.path.join(directory, n+ '_C'+str(channel)+'_'+tilemap[tile]+'.ome.tif'), 
                            pixels_physical_size=im.physical_pixel_sizes, 
                            dimension_order="ZYX")
    print("Finished at " + time.ctime())
     

if __name__=='__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input',type=str,default=None,help='Full path to CZI file to parse')
    args = p.parse_args()
    parseTiles(args.input)
