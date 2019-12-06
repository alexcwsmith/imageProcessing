#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 10:44:53 2019

@author: smith
"""

from PIL import Image
import numpy
import skimage
import os
from skimage import io; io.use_plugin('matplotlib')

BaseDirectory = raw_input('Enter Directory containing files: ') #'/home/smith/Desktop/'

def get_imlist(BaseDirectory):
  """  Returns a list of filenames for
    all jpg images in a directory. """

  return [os.path.join(BaseDirectory, f) for f in os.listdir(BaseDirectory) if f.endswith('.tif')]

files = get_imlist(BaseDirectory)
print files

for img in files:
    im = Image.open(img)
    stack = numpy.stack(im, axis=3).astype("uint8")
    np_im = numpy.array(stack, dtype=numpy.uint8)
    print np_im.shape
    bit_im = skimage.img_as_ubyte(np_im, force_copy=True)
    skimage.io.imsave(img + '8bit.tif', bit_im)

