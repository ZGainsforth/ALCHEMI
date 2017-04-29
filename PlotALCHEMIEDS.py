## Created 2015, Zack Gainsforth
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
from scipy.optimize import curve_fit, fmin
from QuantSpectra import Section, GaussianComponent, LinearComponent
import h5py
import QuickPlot
import ALCHEMIConfig
import shutil, os

def LoadEMDFileEDS(StackName):
    f = h5py.File(StackName, 'r')
    # Access the stack as [m,n,energy]
    EDSStack = f['data']['ALCHEMI EDS']['data']

    # Load all the attributes too and print them out.
    for a in f['microscope'].attrs.keys():
        if a == 'EDX Delta':
            dEDS = float(f['microscope'].attrs[a])
        if a == 'EDX Offset':
            EDSOffset = float(f['microscope'].attrs[a])
        if a == 'EDX Spectrum Length':
            SLength = float(f['microscope'].attrs[a])
        print a + ": " + str(f['microscope'].attrs[a])

    # Generate an E axis for the spectra.
    dEDS = 10
    E = np.arange(0, SLength*dEDS, dEDS)

    return EDSStack, E, f['microscope'].attrs

def ProcessEDSStack(E, EDSStack, OutputPath):
    # First just grab the first spectrum so we can set up the sections for quant.
    DummySpectrum = np.vstack((E, EDSStack[0,0,:]))
    FullSpectrum = np.vstack((E, np.sum(np.sum(EDSStack, axis=0), axis=0)))

    # First we define the sections using the full spectrum.
    Sections = ALCHEMIConfig.DefineSections(FullSpectrum)

    # Show the user how the fits look.
    for n in Sections:
        print n.Name
        n.PlotSectionGraphically(OutputPath)

    # Now lock the FWHM of all the gaussians so they don't wander around.
    for s in Sections:
        for c in s.ComponentDict.values():
            if c.__class__.__name__ == 'GaussianComponent':
                c.FWHMLock = True

    #return

    # Figure out how many elements in total we are tracking.
    ICPNames = []
    for s in Sections:
        for c in s.ComponentDict.values():
            ICPNames.append(c.Name)
    ICPNames.append('Total counts')

    # Allocate a numpy array to hold the various quants.
    ICPs = np.zeros((len(ICPNames), EDSStack.shape[0], EDSStack.shape[1]))

    for m in range(EDSStack.shape[0]):
        print "Processing row %d" % m
        for n in range(EDSStack.shape[1]):
            Spectrum = np.vstack((E, EDSStack[n,m,:]))
            #Sections = ALCHEMIConfig.DefineSections(Spectrum)
            k = 0 # Index into the list of ICPs
            for s in Sections:
                s.UpdateSpectrum(ExperimentalSpectrum=Spectrum)
                s.Fit()
                # Printing the section graphically here can be a performance nightmare but really useful to see how it is fitting all the peaks.
                #s.PlotSectionGraphically(OutputPath=OutputPath+'/Ti/', FileName='%d,%d.png'%(m,n))
                for c in s.ComponentDict.values():
                    ICPs[k, m, n] = c.Area
                    if c.Name != ICPNames[k]:
                        print "Error in indexing the ICPs"
                    k+=1
            ICPs[k, m, n] = sum(EDSStack[n,m,:])
            if 'Total counts' != ICPNames[k]:
                print "Error in indexing total counts ICPs"

    for k in range(len(ICPNames)):
        plt.figure()
        plt.imshow(ICPs[k,:,:], interpolation='none')
        plt.title(ICPNames[k])
        plt.colorbar()
        plt.savefig(os.path.join(OutputPath, ICPNames[k] + '.png'))
        np.savetxt(os.path.join(OutputPath, ICPNames[k] + '.txt'), ICPs[k,:,:])


def PlotMeanAndStdOfStack(E, EDSStack, OutputPath):
    # Reshape the stack so we can compute mean and average easily.
    s = EDSStack.shape
    TempStack = np.reshape(EDSStack, (s[0] * s[1], s[2]))
    # Mean spectrum
    SMean = np.mean(TempStack, axis=0)
    SStd = np.std(TempStack, axis=0)
    # Make a plot showing the range of spectrum intensities, 2 sigma
    (fig, ax) = plt.subplots()
    ax.fill_between(E, SMean - SStd * 2, SMean + SStd * 2, facecolor='red', alpha=0.5)
    # Plot the mean onto it and pretty it up.
    QuickPlot.QuickPlot(E, SMean, xlim=[0, 10000], xlabel='eV', ylabel='Counts',
                        title='Full Spectrum View, 2$\sigma$ variation', figax=(fig, ax))
    plt.savefig(os.path.join(OutputPath, 'StackMeanAndStd.png'))


if __name__ == "__main__":

    # Read in the hd5 file.
    EDSStack, E, MicroscopeAttrs = LoadEMDFileEDS(ALCHEMIConfig.EDSFileName)

    # Make an output directory
    OutputPath, _ = os.path.split(ALCHEMIConfig.EDSFileName)
    OutputPath = os.path.join(OutputPath, 'ICPs')
    if not os.path.exists(OutputPath):
        os.mkdir(OutputPath)

    PlotMeanAndStdOfStack(E, EDSStack, OutputPath)

    ProcessEDSStack(E, EDSStack, OutputPath)

    shutil.copy('ALCHEMIConfig.py', OutputPath)
    plt.show()
