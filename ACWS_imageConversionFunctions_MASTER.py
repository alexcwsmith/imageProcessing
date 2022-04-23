#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 12:36:18 2022

@author: smith
"""

import os
from skimage.io import imread, imsave, imread_collection
from skimage import io
from skimage.util import img_as_uint, img_as_float32, img_as_float64
from skimage.transform import resize
import numpy as np
import h5py
from pathos.pools import ProcessPool
import time
import glob

def horizontalToCoronal(imPath):
    """Transpose horizontal orientation to coronal

    Parameters
    ----------
    imPath : str
        Full path to image to transpose.

    Returns
    -------
    None.

    """
    f, e = os.path.splitext(imPath)
    if f.endswith('.ome'):
        f = f.strip('.ome')
        e = '.ome'+e
    save_path = f+'_Coronal'+e
    im = imread(imPath)
    print("Image loaded with size " + str(im.shape))
    im_corr = im.transpose(1,2,0)
    io.imsave(save_path, im_corr, check_contrast=False)

def horizontalToSagittal(imPath):
    """Transpose horizontal orientation to sagittal"    

    Parameters
    ----------
    imPath : str
        Full path to image to transpose.

    Returns
    -------
    None.

    """
    f, e = os.path.splitext(imPath)
    if f.endswith('.ome'):
        f = f.strip('.ome')
    save_path = f+'_Sagittal.tif'
    im = imread(imPath)
    print("Image loaded with size " + str(im.shape))
    im_sag = im.transpose(2,0,1)
    io.imsave(save_path, im_sag, check_contrast=False)

def sagittalToHorizontal(imPath):
    """Transpose sagittal orientation to horizontal"    

    Parameters
    ----------
    imPath : str
        Full path to image to transpose.

    Returns
    -------
    None.

    """

    f, e = os.path.splitext(imPath)
    if f.endswith('.ome'):
        f = f.strip('.ome')
    save_path = f+'_Horizontal.tif'
    im = imread(imPath)
    print("Image loaded with size " + str(im.shape))
    im_sag = im.transpose(1,2,0)
    io.imsave(save_path, im_sag, check_contrast=False)

def loadDir(path, ext='.tif'):
    return sorted(os.path.join(path, x) for x in os.listdir(path) if x.endswith(ext))

def changeBitDepth(image, target=np.uint16, save=True):
    """Change bit depth of image

    Parameters
    ----------
    image : np.ndarray
        Image array to convert.
    target : np.dtype, optional
        dtype to convert to. The default is np.uint16.
    save : bool, optional
        Whether to save result image. The default is True.

    Raises
    ------
    NameError
        If target is not np.uint8, np.uint16, np.float32, or np.float64

    Returns
    -------
    im : array
        Result array, returned if save=False.

    """
    im = imread(image)
    if im.dtype == target:
        print("Image is already dtype " + str(target))
        return
    else:
        if target==np.uint16:
            newType = '16bit'
            im = img_as_uint(im)
        elif target==np.float32:
            newType = '32bit'
            im = img_as_float32(im)
        elif target==np.float64:
            newType = '64bit'
            im = img_as_float64(im)
        elif target==np.uint8 or target==np.int8:
            newType = '8bit'
            im = im.astype(np.uint8)
        else:
            raise NameError("Target datatype must be np.uint8, np.uint16, np.float32, or np.float64")
    if not save:
        return im
    elif save:
        dirname = os.path.dirname(image)
        n, e = os.path.splitext(image)
        if not os.path.exists(os.path.join(dirname, newType)):
            os.mkdir(os.path.join(dirname, newType))
        imsave(os.path.join(dirname, newType, newType+'_'+n+e), im, check_contrast=False)

def _multiStackToFiles(plane, stackFile, channel, planeNum=0):
    name = os.path.basename(stackFile)
    n, e = os.path.splitext(name)
    directory = os.path.dirname(stackFile)
    
    if stackFile.endswith('.ims') or stackFile.endswith('.h5'):
        file = h5py.File(stackFile, 'r')
        data = file.get('DataSet/ResolutionLevel 0/TimePoint 0/Channel ' + str(channel)+'/Data')
        imNumPadded = str(plane).zfill(4)
        imsave(os.path.join(directory, n + '_C'+str(channel)+'/' + n + '_C' + str(channel) + '_Z' + imNumPadded + '.tif'), np.array(data[plane]), check_contrast=False)
    else:
        imNumPadded = str(planeNum).zfill(4)
        imsave(os.path.join(directory, n + '_C'+str(channel)+'/' + n + '_C' + str(channel) + '_Z' + imNumPadded + '.tif'), np.array(plane), check_contrast=False)

def multiStackToFiles(stackFile, channel, nthreads):
    """Multiprocess extraction of individual tif planes from a 3D stack.
    
    stackFile : str
        Full path to 3D stack file. Can be .tif, .ims, or .h5
    channel : str
        Channel index to extract
    nthreads : int
        Number of threads for multiprocessesing
    """
    name = os.path.basename(stackFile)
    n, e = os.path.splitext(name)
    directory = os.path.dirname(stackFile)
    if not os.path.exists(os.path.join(directory, n + '_C' + str(channel)+'/')):
        os.mkdir(os.path.join(directory, n + '_C' + str(channel)+'/'))
    print("Started at " + time.ctime())
    if stackFile.endswith('.ims') or stackFile.endswith('.h5'):
        file = h5py.File(stackFile, 'r')
        data = file.get('DataSet/ResolutionLevel 0/TimePoint 0/Channel ' + str(channel) + '/Data')
        planes = list(range(data.shape[0]))
    else:
        data = imread(stackFile)
    print("Loaded data of shape " + str(data.shape))
    pool = ProcessPool(nodes=nthreads)
    if stackFile.endswith('.ims') or stackFile.endswith('.h5'):
        pool.map(_multiStackToFiles, list(planes), list([stackFile]*len(planes)), list([channel]*len(planes)))
        pool.close()
        pool.join()
    else:        
        pool.map(_multiStackToFiles, list(data), list([stackFile]*data.shape[0]), list([channel]*data.shape[0]), list(range(data.shape[0])))
        pool.close()
        pool.join() 
    print("Finished at " + time.ctime())

def _IMStoTIF(filePath):
    name, ext = os.path.splitext(filePath)
    if ext != '.ims':
        raise TypeError('Input filePath must be .ims file')      
    file = h5py.File(filePath, 'a')
    data = file.get('DataSet').get('ResolutionLevel 0/TimePoint 0/Channel 0/Data')
    print("Saving " + name + ".tif")
    imsave(name + '.tif', np.array(data), bigtiff=True)

def multiprocessIMStoTIF(args):
    """CLI tool to extract TIF planes from Imaris .ims file.
    
    Parameters
    ----------
    args : argparse.Namespace
        Namespace consisting of args.directory, args.nthreads
    """
    print("Starting at " + str(time.ctime()))
    processes = int(args.nthreads)
    queue = loadDir(args.directory, ext='.ims')
    os.chdir(str(args.directory))
    pool = ProcessPool(nodes=processes)
    pool.map(_IMStoTIF, queue)
    pool.close()
    pool.join()
    print('Completed at ' + str(time.ctime()))


def TIFtoH5(imPath, key='Data'):
    """Convert .tif file to .h5
    
    Parameters
    ----------
    imPath : str
        Full path to .tif image.
    key : str, optional
        h5 key to save. Default is 'Data'
    """
    f, e = os.path.splitext(imPath)
    if f.endswith('.ome'):
        f = f.strip('.ome')
    save_path = f+'.h5'
    hf = h5py.File(save_path, 'a')
    im = imread(imPath)
    print("Image loaded with size " + str(im.shape))
    hf.create_dataset(key, data=im)
    hf.close()
    print("Wrote h5 file with key '" + str(key) + "' to " + str(save_path))

def h5toTIF(imPath, key='Data'):
    """Convert .h5 file to .tif
    
    Parameters
    ----------
    imPath : str
        Full path to .h5 image.
    key : str, optional
        h5 key to load. Default is 'Data'
    """

    if not imPath.endswith('.h5'):
        raise NameError("Only .h5 files supported at this time.")
    f, e = os.path.splitext(imPath)
    if f.endswith('.lux'):
        f = f.strip('.lux')
    save_path = f+'.tif'
    hf = h5py.File(imPath, 'r')
    data = hf.get(key)
    if not data:
        print("Invalid h5 key, detected keys are: " + str(hf.keys()))
    print("Image loaded with size " + str(data.shape))
    io.imsave(save_path, np.array(data), bigtiff=True, check_contrast=False)
    print("Wrote tif file to " + str(save_path))


def filesToStack(directory, load_pattern, sampleName, zindices=None, 
                 crop=False, saveTIF=True, savePlanes=True, memmap=True,
                 conserve_mem=True, returnArray=False):
    """Create 3D stack from tif series.
    
    Parameters
    ----------
    directory : str
        Directory containing .tif series with numeric labels.
    load_pattern : str
        File name pattern to retrieve.
    sampleName : str
        Name to give output file.
    zindices : tuple (ints), optional
        Indices of z-planes to retrieve, if subsampling.
    crop : bool or dict
        If cropping data, a dictionary with {'width', 'height', 'x', 'y'}, ImageJ-format cropping coordinates.
    saveTIF : bool, optional
        Whether to save the final TIF stack, if False returns as array.
    memmap : bool, optional
        Whether to use memory-mapping to limit RAM consumption. Default True.
    conserve_mem : bool, optional
        Whether to use the conserve_mem flag for skimage io. Default True
    returnArray : bool, optional
        Whether to return the final array. Default False.
    """
    files = sorted(glob.glob(os.path.join(directory, load_pattern)))
    if zindices:
        files = files[zindices[0]:zindices[1]]
        print("Processing Z planes " + str(zindices[0]) + ' to ' + str(zindices[1]))        
    im = imread_collection(files, conserve_memory=conserve_mem)
    imc = im.concatenate()
    print("Loaded " + str(imc.shape[0]) + " images with XY dimensions " + str(imc.shape[2]) + "x" + str(imc.shape[1]))
    if crop:
        if not isinstance(crop, dict):
            raise ValueError("Crop argument must be dict of {width, height, x, y} - i.e. ImageJ-style crop coordinates")
        imc = imc[:,crop['y']:crop['y']+crop['height'],crop['x']:crop['x']+crop['width']]
        print("Cropped to XY dimensions " + str(crop['width']) +'x'+ str(crop['height']), "starting at X" + str(crop['x']) + ' Y'+str(crop['y']))
    if memmap:
        if zindices:
            arrName = os.path.join(directory, sampleName+'_Z'+str(zindices[0])+'-Z'+str(zindices[1])+'.npy')
        else:
            arrName = os.path.join(directory, '../'+sampleName+'.npy')
        np.save(arrName, imc)
        if not returnArray:
            del imc
        if saveTIF:
            fname, ext = os.path.splitext(arrName)
            if memmap:
                im = np.load(arrName, allow_pickle=True, mmap_mode='r+')
            if savePlanes:
                if not os.path.exists(fname + '/'):
                    print("Making directory " + fname + "/'")
                    os.mkdir(fname + '/')
                    saveDir=fname+'/'
                for i in range(im.shape[0]):
                    if zindices:
                        if zindices[0]!=0:
                            plane=str(zindices[0]+i).zfill(4)
                            bname=os.path.basename(fname)
                            io.imsave(os.path.join(saveDir, bname+'_'+plane+'.tif'), im[i], check_contrast=False)
                    elif not zindices:
                        plane = str(i).zfill(4)
                        io.imsave(fname+'_'+plane+'.tif', im[i], check_contrast=False, plugin='tifffile', imagej=True)
            elif not savePlanes:
                """imagej argument docs:
                    If True, write an ImageJ hyperstack compatible file. 
                    This format can handle data types uint8, uint16, or float32 and data shapes up to 
                    6 dimensions in TZCYXS order. RGB images (S=3 or S=4) must be uint8. 
                    ImageJ's default byte order is big endian but this implementation uses the 
                    system's native byte order by default. ImageJ does not support BigTIFF format or 
                    LZMA compression. The ImageJ file format is undocumented.
                """
                im = im[np.newaxis,:,np.newaxis,:,:]
                io.imsave(fname + '_stack.tif', im, check_contrast=False, plugin='tifffile', imagej=True)
    elif not memmap:
        if saveTIF:
            io.imsave(os.path.join(directory, '../'+sampleName+'.tif'), imc, plugin='tifffile', imagej=True)           
    if returnArray:
        return imc

def downsampleZstack(factor, imagePath=None, im=None, name=None, save=True):
    """Downsample a Z-stack.
    
    Parameters
    ----------
    factor : flaot
        Factor by which to downsample by.
    imagePath : str, optional
        Path to image to load, default None.
    im : np.array, optional
        Array of previously loaded image. Cannot pass im and imagePath
    name : str, optional
        File name to save if passing im.
    save : bool, optional
        Whether to save result, Default True. False returns array.
    """
    if imagePath:
        directory = os.path.dirname(imagePath)
        n, e = os.path.splitext(imagePath)
        im = io.imread(imagePath)
    elif im:
        n = name
        e = '.tif'
    downZ = im.shape[0]//factor
    arrays=[]
    for i in range(downZ):
        sample = im[i:i+downZ,:,:]
        maxproj = np.max(sample, axis=0)
        arrays.append(maxproj)
    flatImage = io.concatenate_images(arrays)
    if save:
        io.imsave(os.path.join(directory, n+'_DownsampledZ'+str(factor)+e), flatImage)
    elif not save:
        return flatImage

def _downsample(imagePath, size=(12, 256,256), mode='wrap', preserve_range=True, anti_aliasing=True):
    im = io.imread(imagePath)
    directory = os.path.dirname(imagePath)
    n, e = os.path.splitext(os.path.basename(imagePath))
    imres = resize(im, size, mode=mode, preserve_range=preserve_range, anti_aliasing=anti_aliasing)
    io.imsave(os.path.join(directory, n + '_downsampled'+e), imres)

def downsample(directory, ext='.tif', nthreads=24):
    """Downsample individual planes of TIF stack with parallel processing.
    
    Parameters
    ----------
    directory : str
        Path to directory containing tif series.
    ext : str
        Extension of images to load, default '.tif'
    nthreds : int
        Number of threads for multiprocessing.
    """
    queue = loadDir(directory, ext)
    pool = ProcessPool(nodes=nthreads)
    pool.map(_downsample, queue)
    pool.close()
    pool.join()

