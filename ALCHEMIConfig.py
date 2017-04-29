from QuantSpectra import CreateEDSFitSection

'''
This is the only file you should have to edit for one ALCHEMI run to the next.  It contains the few lines of code that define what peaks we want to fit for this spectrum, and what the fitting ranges are.
All the non-hardcoded stuff should go there.
'''

#CCDFileName = '/Volumes/Zack/Desktop/20160628 - TitanX - Righter Ti-TiO2 spinel FIB 1/ALCHEMI 2/CCD/Stack.emd'
#MosaicBinning = 5

EDSFileName = '/Volumes/Desktop/20170412 - TitanX - GRA 95229 - Chondrule 1 - Altered chromite/ALCHEMI 5 - Altered region/EDS/Stack.emd'
def DefineSections(Spectrum):
    Sections = []
    #Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='Cr-Fe', SectionRange=[5000,7800], PeakNames=['CrKa', 'CrKb', 'FeKa', 'FeKb'], BkgName='Bkg'))
    #Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='Si', SectionRange=[1600,1950], PeakNames=['SiK'], BkgName='Bkg'))
    Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='Ti', SectionRange=[4290,4750], PeakNames=['TiKa'], BkgName='Bkg'))
    Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='V', SectionRange=[4800,5100], PeakNames=['VKa'], BkgName='Bkg'))
    Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='Mg', SectionRange=[1165,1367], PeakNames=['MgK'], BkgName='Bkg'))
    Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='Al', SectionRange=[1350,1600], PeakNames=['AlK'], BkgName='Bkg'))
    Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='O', SectionRange=[400,670], PeakNames=['OK'], BkgName='Bkg'))
    Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='Cr-Fe', SectionRange=[5000,7300], PeakNames=['CrKa', 'CrKb', 'FeKa', 'FeKb'], BkgName='Bkg'))
    return Sections

# EDSFileName = '/Users/Zack/Desktop/ALCHEMI 2 on chromite - 110 zone/StackEDS.emd'
#
# def DefineSections(Spectrum):
#     Sections = []
#     Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='Cr-Fe', SectionRange=[5000,7300], PeakNames=['CrKa', 'CrKb', 'FeKa', 'FeKb'], BkgName='Bkg'))
#     Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='Al', SectionRange=[1350,1600], PeakNames=['AlK'], BkgName='Bkg'))
#     Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='O', SectionRange=[400,670], PeakNames=['OK'], BkgName='Bkg'))
#
#     return Sections

# EDSFileName = '/Users/Zack/Desktop/Pedestal Pyrrhotite/ALCHEMI/Stack EDS.emd'
#
# def DefineSections(Spectrum):
#     Sections = []
#     Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='Cr-Cu', SectionRange=[5000,9000], PeakNames=['CrKa', 'CrKb', 'FeKa', 'FeKb', 'NiKa', 'NiKb', 'CuKa', 'CuKb',], BkgName='Bkg'))
#     Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='S', SectionRange=[2150,2400], PeakNames=['SK'], BkgName='Bkg'))
#     Sections.append(CreateEDSFitSection(SpectrumIn=Spectrum, SectionName='O', SectionRange=[400,670], PeakNames=['OK'], BkgName='Bkg'))
#
#     return Sections
