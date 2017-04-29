## Created 2015, Zack Gainsforth
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
from scipy.optimize import curve_fit, fmin
from QuantSpectra import Section, GaussianComponent, LinearComponent
import h5py
import QuickPlot
import ALCHEMIConfig
import shutil, os
from scipy.ndimage.interpolation import zoom

def LoadEMDFileCCD(StackName):
    f = h5py.File(StackName, 'r')
    # Access the stack as [m,n,energy]
    CCDStack = f['data']['ALCHEMI CCD']['data']

    # Load all the attributes too and print them out.
    for a in f['microscope'].attrs.keys():
        print a + ": " + str(f['microscope'].attrs[a])

    return CCDStack, f['microscope'].attrs

def ProcessCCDStack(CCDStack, MosaicBinning, OutputPath):

    BinnedImageShape = zoom(CCDStack[0,0,:,:], 1./MosaicBinning).shape

    StdImg = np.zeros(CCDStack.shape[0:2])
    AveImg = np.zeros(CCDStack.shape[0:2])
    MosaicImg = np.zeros((CCDStack.shape[0]*BinnedImageShape[0], CCDStack.shape[1]*BinnedImageShape[1]))

    for m in range(CCDStack.shape[0]):
        print "Processing row %d" % m
        for n in range(CCDStack.shape[1]):
            AveImg[m,n] = np.mean(CCDStack[m,n,:,:], axis=(0,1))
            StdImg[m,n] = np.std(CCDStack[m,n,:,:], axis=(0,1))
            MosaicImg[m*BinnedImageShape[0]:(m+1)*BinnedImageShape[0], n*BinnedImageShape[1]:(n+1)*BinnedImageShape[1]] = zoom(CCDStack[m,n,:,:], 1./MosaicBinning)

    plt.figure()
    plt.imshow(AveImg, interpolation='none')
    plt.title('Ave image')
    plt.colorbar()
    plt.savefig(os.path.join(OutputPath, 'Ave Image.png'))
    np.savetxt(os.path.join(OutputPath, 'Ave Image.txt'), AveImg)

    plt.figure()
    plt.imshow(StdImg, interpolation='none')
    plt.title('StDev image')
    plt.colorbar()
    plt.savefig(os.path.join(OutputPath, 'StDev Image.png'))
    np.savetxt(os.path.join(OutputPath, 'StDev Image.txt'), StdImg)

    plt.figure()
    plt.imshow(MosaicImg, interpolation='none', cmap='gray', extent=(0,MosaicImg.shape[0], 0,MosaicImg.shape[1]))
    plt.title('Mosaic image')
    plt.savefig(os.path.join(OutputPath, 'Mosaic Image.png'))
    np.savetxt(os.path.join(OutputPath, 'Mosaic Image.txt'), MosaicImg)

def PrintCCDCenterAndEdgeImages(CCDStack, OutputPath):
    # Get a short variable name for the CCD stack shape.
    s = CCDStack.shape

    # Allocate an image large enough for nine images.  The center, edges and corners.
    MosaicImg = np.zeros((3*s[2], 3*s[3]))
    for i, m in     [(0,0), (1, s[0]/2), (2,-1)]:
        for j, n in [(0,0), (1, s[1]/2), (2,-1)]:
            # Put that image in the mosaic.
            MosaicImg[i*s[2]:(i+1)*s[2], j*s[3]:(j+1)*s[3]] = CCDStack[m,n,:,:]

    plt.figure()
    plt.imshow(MosaicImg, interpolation='none', cmap='gray', extent=(0,MosaicImg.shape[0], 0,MosaicImg.shape[1]))
    plt.title('Center and Edge Mosaic')
    plt.savefig(os.path.join(OutputPath, 'Center and Edge Mosaic.png'))
    np.savetxt(os.path.join(OutputPath, 'Center and Edge Mosaic.txt'), MosaicImg)

    plt.figure()
    CenterImage = CCDStack[s[0]/2,s[1]/2,:,:]
    plt.imshow(CenterImage, interpolation='none', cmap='gray', extent=(0,s[2], 0,s[3]))
    plt.title('Center Image')
    plt.savefig(os.path.join(OutputPath, 'Center Image.png'))
    np.savetxt(os.path.join(OutputPath, 'Center Image.txt'), CenterImage)




# def PlotMeanAndStdOfStack(E, CCDStack, OutputPath):
#     # Reshape the stack so we can compute mean and average easily.
#     s = CCDStack.shape
#     TempStack = np.reshape(CCDStack, (s[0]*s[1], s[2]*s[3]))
#     # Mean
#     SMean = np.mean(TempStack, axis=0)
#     SStd = np.std(TempStack, axis=0)
#     # Make a plot showing the range of spectrum intensities, 2 sigma
#     (fig, ax) = plt.subplots()
#     ax.fill_between(E, SMean - SStd * 2, SMean + SStd * 2, facecolor='red', alpha=0.5)
#     # Plot the mean onto it and pretty it up.
#     QuickPlot.QuickPlot(E, SMean, xlim=[0, 10000], xlabel='eV', ylabel='Counts',
#                         title='Full Spectrum View, 2$\sigma$ variation', figax=(fig, ax))
#     plt.savefig(os.path.join(OutputPath, 'StackMeanAndStd.png'))


if __name__ == "__main__":

    # Read in the hd5 file.
    CCDStack, MicroscopeAttrs = LoadEMDFileCCD(ALCHEMIConfig.CCDFileName)

    # Make an output directory
    OutputPath, _ = os.path.split(ALCHEMIConfig.CCDFileName)
    OutputPath = os.path.join(OutputPath, 'CCDViews')
    if not os.path.exists(OutputPath):
        os.mkdir(OutputPath)

    #PlotMeanAndStdOfStack(E, EDSStack, OutputPath)

    ProcessCCDStack(CCDStack, ALCHEMIConfig.MosaicBinning, OutputPath)
    PrintCCDCenterAndEdgeImages(CCDStack, OutputPath)

    shutil.copy('ALCHEMIConfig.py', OutputPath)
    plt.show()