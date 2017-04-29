## Created 2015, Zack Gainsforth
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
from scipy.optimize import curve_fit, fmin
import QuickPlot
import os
#from numba import jit


verbose=False
def MyPrint(PrintStr):
    #from __future__ import print_function
    if verbose==True:
        print PrintStr

# First we need a list of different fit elements.  We need to fit bounded portions of the spectrum so we don't go
# down the rabbithole of trying to have perfect Brehmsstrahlung fit and such.  However, with overlapping peaks, or EELS
# we also need to be able to fit multiple peaks at one time.  So the highest tier is sections which contain components.

# This is the base component.  The only thing they have in common is unique name
class Component:
    Formula = 'x'

    def GetPlot(self, x):
        y = eval(self.Formula)
        # Constrain everything to be positive.
        y[y<0] = 0
        return y

    def __init__(self, Name=None):
        if Name==None:
            print 'Component needs a name, such as FeKa, FeEdge, or Bkg'
        self.Name = Name

    def __str__(self):
        # Always print out the basic class info and name first.
        PrintStr = "Class %s():\n" % self.__class__.__name__
        PrintStr += "\tName: %s\n" % self.Name

        # Then print out all the other local variables.
        LocalVars = dir(self)
        for v in LocalVars:
            # Ignore the python inherent variables.
            if v[0] == '_':
                continue
            # Ignore other variables that would clutter up the output.
            if v in ['Name', 'DefaultNames']:
                continue
            # Ignore methods.
            if callable(getattr(self, v)):
                continue
            # OK, it passes the filters, let's print it.
            PrintStr += "\t%s: %s\n" % (v, eval('self.'+v))

        return PrintStr

class ConstComponent(Component):
    Formula = 'x-x+self.Offset'
    sticky_i = 0

    def Guess(self, x, y):
        # We guess the line using a polyfit
        self.Offset = np.mean(y)
        # if True:#self.sticky_i == 0:
        #     print P, self.sticky_i
        #     self.sticky_i += 1

    def Fit(self, x, y):
        MyPrint( "Fitting " + self.Name)
        self.Guess(x,y) # Guess is already a polyfit ... that's about the speed we want!
        self.Area = np.sum(self.GetPlot(x))

    def __init__(self, Name=None, Offset=None):
        self.Name = 'Default'
        self.Offset = 0

        # Override any specific values only if the user passed them in.
        if Name is not None:
            self.Name = Name
        if Offset is not None:
            self.Offset = float(Offset)

class LinearComponent(Component):
    Formula = 'self.Slope*x + self.Offset'
    sticky_i = 0

    def Guess(self, x, y):
        # We guess the line using a polyfit
        P = np.polyfit(x, y, 1)
        self.Slope = P[0]
        self.Offset = P[1]
        # if True:#self.sticky_i == 0:
        #     print P, self.sticky_i
        #     self.sticky_i += 1

    def Fit(self, x, y):
        MyPrint( "Fitting " + self.Name)
        self.Guess(x,y) # Guess is already a polyfit ... that's about the speed we want!
        self.Area = np.sum(self.GetPlot(x))

    def __init__(self, Name=None, Offset=None, Slope=None):
        self.Name = 'Default'
        self.Offset = 0
        self.Slope = 0

        # Override any specific values only if the user passed them in.
        if Name is not None:
            self.Name = Name
        if Offset is not None:
            self.Offset = float(Offset)
        if Slope is not None:
            self.Slope = float(Slope)

class ExpComponent(Component):
    Formula = 'self.Amp*np.exp(x/self.Decay)'

    def Guess(self, x, y):
        print 'No Guess for exp yet!'
    def Fit(self, x, y):
        print 'No Fit for exp yet!'

    def __init__(self, Name=None, Amp=None, Decay=None):
        self.Name = 'Default'
        self.Amp = 0
        self.Decay = 1

        # Override any specific values only if the user passed them in.
        if Name is not None:
            self.Name = Name
        if Amp is not None:
            self.Amp = float(Amp)
        if Decay is not None:
            self.Decay = float(Decay)

class GaussianComponent(Component):
    Formula = 'self.Amp*np.exp(-((x-self.Center)**2/(self.FWHM**2*np.sqrt(2))))'
    FWHMLock = False
    CenterLock = False

    def Guess(self, x, y):
        # Given an x, and y for the experimental function, guess what the right amplitude is.
        if np.max(x) < self.Center or np.min(x) > self.Center:
            MyPrint( "Warning when guessing %s, given range doesn't include peak center." % self.Name)
            # We will continue, but the peak height could be off.
        # We guess the value of the function where x is closest to the center of the peak.
        self.Amp = y[np.argmin(np.abs(x-self.Center))]
        self.Area = 0

#    @jit(nopython=True)
    def FitFunc(self, x, Center, Amp, FWHM):
        if self.FWHMLock == True:
            FWHM = self.FWHM
        if self.CenterLock == True:
            Center = self.Center
        # Constrain the Amplitude to be positive.
        Penalization=1
        if Amp < 0:
            Penalization=Amp+1
        return Amp*np.exp(-((x-Center)**2/(FWHM**2*np.sqrt(2))))-Penalization**2

    def Fit(self, x, y):
        MyPrint( "Fitting " + self.Name)
        try:
            FitParams, FitCov = curve_fit(self.FitFunc, x, y, (self.Center, self.Amp, self.FWHM))
            self.Center = FitParams[0]
            self.Amp = FitParams[1]
            self.FWHM = FitParams[2]
            self.FitParams = FitParams
            self.FitCov = FitCov
            self.Area = np.sum(self.FitFunc(x, *FitParams))
        except:
            MyPrint( "Fitting " + self.Name + " failed.\n")

    # List of default values for known gaussian peaks.  Center, Amp, FWHM.
    DefaultNames = OrderedDict()
    DefaultNames['CK'] =   [280, 1, 30]
    DefaultNames['OK'] =   [530, 1, 40]
    DefaultNames['FeL'] =  [710, 1, 50]
    DefaultNames['CuL'] =  [935, 1, 50]
    DefaultNames['MgK'] =  [1270, 1, 50]
    DefaultNames['AlK'] =  [1490, 1, 50]
    DefaultNames['SiK'] =  [1740, 1, 35]
    DefaultNames['PK'] =   [2010, 1, 50]
    DefaultNames['SK'] =   [2310, 1, 50]
    DefaultNames['ClK'] =  [2621, 1, 60]
    DefaultNames['KK'] =   [3313, 1, 60]
    DefaultNames['CaK'] =  [3691, 1, 60]
    DefaultNames['TiKa'] = [4510, 1, 60]
    DefaultNames['TiKb'] = [4933, 1, 60]
    DefaultNames['VKa'] =  [4950, 1, 60]
    DefaultNames['VKb'] =  [5428, 1, 60]
    DefaultNames['CrKa'] = [5410, 1, 80]
    DefaultNames['CrKb'] = [5950, 1, 80]
    DefaultNames['FeKa'] = [6404, 1, 80]
    DefaultNames['FeKb'] = [7055, 1, 80]
    DefaultNames['NiKa'] = [7476, 1, 80]
    DefaultNames['NiKb'] = [8266, 1, 80]
    DefaultNames['CuKa'] = [8040, 1, 80]
    DefaultNames['CuKb'] = [8900, 1, 80]

    def __init__(self, Name=None, Center=None, Amp=None, FWHM=None):
        self.Name = 'Default'
        self.Center = 0 # eV
        self.Amp = 1 # counts, usually, but it is whatever units the spectrum has.
        self.FWHM = 100 #eV

        # If the user gives us a name, we will populate the values from the default list of known names.
        # Of course, the user can override any of the default values.
        if Name is not None:
            # Tell the user if his name matches a known name.
            if Name in self.DefaultNames:
                MyPrint( 'Loading default %s values for %s' % (self.__class__.__name__, Name))
                # Populate with the default values for that name.
                self.Name = Name
                self.Center = self.DefaultNames[Name][0]
                self.Amp = self.DefaultNames[Name][1]
                self.FWHM = self.DefaultNames[Name][2]

        # Override any specific values only if the user passed them in.
        if Name is not None:
            self.Name = Name
        if Center is not None:
            self.Center = float(Center)
        if Amp is not None:
            self.Amp = float(Amp)
        if FWHM is not None:
            self.FWHM = float(FWHM)

class EdgeComponent(GaussianComponent):
    Formula = '-self.Amp/np.pi*np.arctan((x-self.Center)/self.FWHM)'

    # List of default values for known gaussian peaks.  Center, Amp, FWHM.
    DefaultNames = OrderedDict()
    DefaultNames['CKEdge'] =   [284.2, 1, 30]
    DefaultNames['OKEdge'] =   [543.1, 1, 40]
    DefaultNames['FeLEdge'] =  [707, 1, 50]
    DefaultNames['CuLEdge'] =  [952, 1, 50]
    DefaultNames['MgKEdge'] =  [1303, 1, 50]
    DefaultNames['AlKEdge'] =  [1559, 1, 50]
    DefaultNames['SiKEdge'] =  [1839, 1, 30]
    DefaultNames['PKEdge'] =   [2146, 1, 50]
    DefaultNames['SKEdge'] =   [2472, 1, 50]
    DefaultNames['ClKEdge'] =  [2822, 1, 60]
    DefaultNames['KKEdge'] =   [3608, 1, 60]
    DefaultNames['CaKEdge'] =  [4038, 1, 60]
    DefaultNames['TiKEdge'] =  [4966, 1, 60]
    DefaultNames['CrKEdge'] = [5989, 1, 80]
    DefaultNames['FeKEdge'] = [7112, 1, 80]
    DefaultNames['CuKEdge'] = [8979, 1, 80]

class Section:

    def TrimSpectrum(self, S=None):
        if S is None:
            # User can trim a spectrum passed in, or he can trim the one stored in the class.
            S = self.ExperimentalSpectrum
        # S should be a 2xn array.
        StartPos = np.argmin(np.abs(S[0,:]-self.StarteV))
        EndPos = np.argmin(np.abs(S[0,:]-self.EndeV))
        return S[:,StartPos:EndPos]

    def GetPlot(self, x=None):
        if x is not None:
            y = np.zeros(len(x))
        else:
            S = self.TrimSpectrum(self.ExperimentalSpectrum)
            x = S[0,:]
            y = S[1,:]

        for c in self.ComponentDict.keys():
            y += self.ComponentDict[c].GetPlot(x)
        return y

    def GetEnergyAxis(self):
        S = self.TrimSpectrum(self.ExperimentalSpectrum)
        x = S[0,:]
        return x

    def PlotSectionGraphically(self, OutputPath=None, FileName=None):
        E = self.GetEnergyAxis()
        # Get the background
        ybkg = self.ComponentDict['Bkg'].GetPlot(E)
        # Get the experimental data
        yraw = self.TrimSpectrum()[1, :]
        # Get the fit data.
        yfit = self.GetPlot(E)
        # Plot these:
        (fig, ax) = QuickPlot.QuickPlot(E, ybkg, boldlevel=2)
        QuickPlot.QuickPlot(E, yraw, boldlevel=2, figax=(fig, ax))
        QuickPlot.QuickPlot(E, yfit, boldlevel=2, figax=(fig, ax),
                            title="Section: Name=%s, StarteV=%d, EndeV=%d" % (self.Name, self.StarteV, self.EndeV),
                            xlim=[self.StarteV, self.EndeV], legendstrs=['Background', 'Experimental', 'Fit'])

        AreaText = ''
        for n in self.ComponentDict.values():
            if getattr(n, 'Area', False):
                AreaText += n.Name + ": " + str(n.Area) + "\n"
        fig.text(0.15,0.7, AreaText)

        if OutputPath is not None:
            if FileName is None:
                FileName = self.Name
            plt.savefig(os.path.join(OutputPath, FileName))

        return (fig,ax)


    def Fit(self, x=None, y=None, Threshold=0.0001):
        MyPrint( 'Fitting section: ' + self.Name)

        # A user can pass in x and y to fit, but if not then we'll use the default portion of our stored spectrum.
        x = self.GetEnergyAxis()
        y = self.TrimSpectrum(self.ExperimentalSpectrum)[1, :]

        # Never do more than 100 fit iterations.  That's just be a hang.
        LastResidual=0
        for n in range(100):
            MyPrint( 'Iteration %d' % n)
            # Use a round-robin fit, fitting just one component at a time.
            for c in self.ComponentDict.keys():
                # First get our fit exempting the current component.
                CurPlot = np.zeros(len(x))
                for i in self.ComponentDict.keys():
                    if i == c:
                        continue
                    CurPlot += self.ComponentDict[i].GetPlot(x)

                # Remove that current plot from the experimental data and fit what's left
                self.ComponentDict[c].Fit(x,y-CurPlot)
            Residual = np.sum(np.abs(y-self.GetPlot(x)))
            MyPrint( 'Residual is %g' % Residual)
            if np.abs((Residual-LastResidual)/Residual) < Threshold:
                MyPrint( 'Change in residual is < %g.  Fit complete.' % Threshold)
                break
            LastResidual = Residual

    def UpdateSpectrum(self, ExperimentalSpectrum=None):
        if ExperimentalSpectrum is not None:
            if ExperimentalSpectrum.shape[0] == 2:
                self.ExperimentalSpectrum = ExperimentalSpectrum
                #print 'Spectrum updated'
            else:
                print 'ExperimentalSpectrum was not a two-row spectrum (should be a 2 x n array, top row is energy, bottom is intensity).'

    def __init__(self, Name=None, StarteV=None, EndeV=None, ComponentDict=None, ExperimentalSpectrum=None):
        self.Name = None
        self.StarteV = 0
        self.EndeV = 1000
        self.ComponentDict = OrderedDict()
        self.ExperimentalSpectrum = None

        if Name==None:
            print 'Section needs a name'
        self.Name = Name
        if StarteV is not None:
            self.StarteV = float(StarteV)
        if EndeV is not None:
            self.EndeV = float(EndeV)
        if ComponentDict is not None:
            self.ComponentDict = ComponentDict
        if ExperimentalSpectrum is not None:
            if ExperimentalSpectrum.shape[0] == 2:
                self.ExperimentalSpectrum = ExperimentalSpectrum
            else:
                print 'ExperimentalSpectrum was not a two-row spectrum (should be a 2 x n array, top row is energy, bottom is intensity).'

    def __str__(self):
        # Always print out the basic class info and name first.
        PrintStr = "Class %s():\n" % self.__class__.__name__
        PrintStr += "\tName: %s\n" % self.Name

        # Then print out all the other local variables.
        LocalVars = dir(self)
        for v in LocalVars:
            # Ignore the python inherent variables.
            if v[0] == '_':
                continue
            # Ignore other variables that would clutter up the output.
            if v in ['Name']:
                continue
            # Ignore methods.
            if callable(getattr(self, v)):
                continue
            # Special printing for the component dictionary.
            if v == 'ComponentDict':
                for c in self.ComponentDict.keys():
                    PrintStr += '\tContains component: %s\n' %c
                continue
            # OK, it passes the filters, let's print it.
            PrintStr += "\t%s: %s\n" % (v, eval('self.'+v))

        return PrintStr

def CreateEDSFitSection(SpectrumIn=None, SectionName=None, SectionRange=None, PeakNames=None, BkgName=None):
    """
    :param SpectrumIn: A spectrum which will be used to fit against.
    :param SectionName: The name of which section.  The exact text is unimportant.  E.g. Cr-Cu tells the user this fits the portion of the spectrum from Cr to Cu.
    :param SectionRange: e.g. [5000,9000] the range of energies to restrict the fit.  The spectrum is usually larger than the fitting range.
    :param PeakNames: A list of strings with GaussianPeaks to insert for known elements.  e.g. ['CrKa', 'CrKb']
    :param BkgNames: The name for the background, or None if no background.
    :return:  Returns the Section class after adding all the peaks and doing the fit.
    """

    # Make a section.
    a = Section(Name=SectionName, StarteV=SectionRange[0], EndeV=SectionRange[1], ExperimentalSpectrum=SpectrumIn)
    x = a.GetEnergyAxis()  # Get the x-axis for the portion of the spectrum bounded by StarteV to EndeV.
    yraw = a.TrimSpectrum(SpectrumIn)[1, :]
    # Let's define the lines we want.
    ComponentDict = OrderedDict()
    a.ComponentDict = ComponentDict

    # To guess the background, we have to first remove all the element lines.
    # if BkgName is not None:
    #a.ComponentDict[BkgName] = ConstComponent(Name=BkgName)
    a.ComponentDict[BkgName] = LinearComponent(Name=BkgName)
    a.ComponentDict[BkgName].Guess(x, yraw - a.GetPlot(x))
    # else:
    #     a.ComponentDict['Bkg'] = LinearComponent(Name='Bkg', Slope=0, Offset=0)

    # Add each line and then background we see into the component dictionary, and have it guess it's amplitude.
    for L in PeakNames:
        a.ComponentDict[L] = GaussianComponent(Name=L)
        a.ComponentDict[L].Guess(SpectrumIn[0, :], SpectrumIn[1, :])

    # And finally, fit it.
    a.Fit(x, yraw, Threshold=0.0001)

    # Return the now created class.
    return a


if __name__ == '__main__':
    # # Test the Gaussian Component class.
    # print "TESTING CLASS GaussianComponent()\n\n"
    # a = GaussianComponent()
    # print a
    # a = GaussianComponent(Name='FeKa')
    # print a
    # a = GaussianComponent(Name='FeKa', Center=10)
    # print a
    # a = GaussianComponent(Name='FeKa', Amp=1000)
    # print a
    # a = GaussianComponent(Name='FeKa', FWHM=10)
    # print a
    # a = GaussianComponent(Name='Oddball')
    # print a
    # a = GaussianComponent(Name='Oddball', Center=3, Amp=4, FWHM=5)
    # print a
    # x = np.linspace(-10,10,500)
    # y = a.GetPlot(x)
    # plt.figure()
    # plt.plot(x,y)
    # plt.title("GaussianComponent(Name='Oddball', Center=3, Amp=4, FWHM=5)")
    #
    # # Test the Edge Component class.
    # print "\nTESTING CLASS EdgeComponent()\n"
    # a = EdgeComponent()
    # print a
    # a = EdgeComponent(Name='FeKEdge')
    # print a
    # a = EdgeComponent(Name='FeKEdge', Center=10)
    # print a
    # a = EdgeComponent(Name='FeKEdge', Amp=1000)
    # print a
    # a = EdgeComponent(Name='FeKEdge', FWHM=10)
    # print a
    # a = EdgeComponent(Name='Oddball')
    # print a
    # a = EdgeComponent(Name='Oddball', Center=3, Amp=4, FWHM=5)
    # print a
    # y = a.GetPlot(x)
    # plt.figure()
    # plt.plot(x,y)
    # plt.title("EdgeComponent(Name='Oddball', Center=3, Amp=4, FWHM=5)")
    #
    # # Test the Linear Component class.
    # print "\nTESTING CLASS LinearComponent()\n"
    # a = LinearComponent()
    # print a
    # b = LinearComponent(Name='Bkg', Slope=5, Offset=10)
    # print b
    # y = b.GetPlot(x)
    # plt.figure()
    # plt.plot(x,y)
    # plt.title("LinearComponent(Name='Bkg', Slope=5, Offset=10)")
    #
    # # Test the Exponential Component class.
    # print "\nTESTING CLASS ExpComponent()\n"
    # a = ExpComponent()
    # print a
    # a = ExpComponent(Name='Bkg', Amp=30, Decay=20)
    # print a
    # y = a.GetPlot(x)
    # plt.figure()
    # plt.plot(x,y)
    # plt.title("ExpComponent(Name='Bkg', Amp=30, Decay=20)")

    # Test the Section class.
    print "\nTESTING CLASS Section()\n"
    S = np.genfromtxt('0.txt').T

    FeCrSection = CreateEDSFitSection(SpectrumIn=S, SectionName='Cr-Cu', SectionRange=[5000,9000], PeakNames=['CrKa', 'CrKb', 'FeKa', 'FeKb', 'CuKa', 'CuKb'], BkgName='Bkg')
    MgAlSection = CreateEDSFitSection(SpectrumIn=S, SectionName='Mg-Al', SectionRange=[1100,1600], PeakNames=['MgK', 'AlK'], BkgName='Bkg')
    OSection = CreateEDSFitSection(SpectrumIn=S, SectionName='O', SectionRange=[400,670], PeakNames=['OK'], BkgName='Bkg')

    for n in [FeCrSection, MgAlSection, OSection]:
        n.PlotSectionGraphically()

    plt.show()
