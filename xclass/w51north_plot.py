#!/usr/bin/env python

################################################################
#
# Python script to model the spectra of W51 (Adam)
#

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import codata
from astropy.io import fits
from astropy.wcs import WCS

plt.close()

################################################################
#
# Definition of some constants
#

C = codata.physical_constants
c = C['speed of light in vacuum'][0]
k = C['Boltzmann constant'][0]
h = C['Planck constant'][0]

################################################################
#
# User-defined variables
#

#
# which spectral chunks to plot

#chunks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
chunks = [0]

#
# setup for XCLASS

FreqStep = 0.1         # frequency resolution of the data (channel width) [in MHz]
Inter_Flag = True      # T/F for interferometric/single dish observations
TelescopeSize = 0.42   # size of the synthesized beam  [in arcsec] 
t_back_flag = True     # include the background temperature
tBack = 45.0           # background temperature [in K]
tslope = 0.0           # slope of the temperature
nH_flag = True         # include hydrogen column density
N_H = 1.0E+24          # hydrogen column density
beta_dust = 0.0        # dust spectral index
kappa_1300 = 0.0       # dust absorption coefficient
RestFreq = 0.0         # rest frequency in MHz (if the data needs to be shifted)
vLSR = 0.0             # line-of-sight velocity in km/s  (if the data needs to be shifted)
iso_flag = False       # include isotoplogues 

IsoTableFileName = "isonames.txt"

#MolfitsFileBaseNames = ["C2H3CN", "C2H5CN", "C2H5OH", "c-C3H2", "CH2OHCHO", "CH3COCH3", "aGg-CH2OH_2", "HCCCN", "C2H5OCHO", "n-C3H7CN", "c-C2H4O", "CCH", "CCS", "CH2NH", "CH3C3N", "CH3CCH", "CH3CHO", "CH3CN", "CH3COOH", "CH3NH2", "CH3OCH3", "CH3OCHO", "CH3OH", "CH3SH", "CS", "CH2CO", "H2CO", "H2CS", "H2S", "HCCNC", "HCN", "NH2CHO", "HCO+", "HCS+", "HDO", "HNC", "HNCO", "HNCS", "HOCN", "HOCO+", "HOC+", "HSCN", "N2H+", "NH2CN", "NH2D", "NS", "OCS", "SiO", "SO2", "SO", "HCOOH"]
#MolfitsFileBaseNames = ["C2H3CN", "C2H5CN", "C2H5OH", "c-C3H2", "CH2OHCHO", "CH3COCH3", "aGg-CH2OH_2", "HCCCN"]
#MolfitsFileBaseNames = ["C2H5OCHO", "n-C3H7CN", "c-C2H4O", "CCH"] #
#MolfitsFileBaseNames = ["CCS", "CH2NH", "CH3C3N", "CH3CCH", "CH3CHO", "CH3CN"] #
#MolfitsFileBaseNames = ["CH3COOH", "CH3NH2", "CH3OCH3", "CH3OCHO", "CH3OH"] #
#MolfitsFileBaseNames = ["CH3SH", "CS", "CH2CO", "H2CO", "H2CS", "H2S", "HCCNC"] #
#MolfitsFileBaseNames = ["HCN", "NH2CHO", "HCO+", "HCS+", "HDO", "HNC", "HNCO"] #
#MolfitsFileBaseNames = ["HNCS", "HOCN", "HOCO+", "HOC+", "HSCN", "N2H+", "NH2CN"] #
#MolfitsFileBaseNames = ["NH2D", "NS", "OCS", "SiO", "SO2", "SO", "HCOOH"] #
#MolfitsFileBaseNames = ["HCCCN"]
#MolfitsFileBaseNames = ["all_molecules"]
#MolfitsFileBaseNames = ["New_all__final__LM__call_1.out"]
MolfitsFileBaseNames = ["w51north"]

#colors = ["Red", "MediumBlue", "DarkRed", "DeepSkyBlue", "MediumVioletRed", "SpringGreen", "DarkViolet","DarkGreen", "Lime", "Chocolate", "SaddleBrown", "DodgerBlue", "DeepPink", "MediumVioletRed", "Teal", "MidnightBlue",  "LimeGreen", "OrangeRed", "ForestGreen", "Goldenrod", "Navy", "Purple", "Maroon", "DarkOrange", "Yellow", "Fuchsia", "Indigo", "GreenYellow", "Aqua", "SaddleBrown", "SlateBlue", "DarkGoldenrod", "HotPink"]
#colors = ["Grey", "LightGrey", "RosyBrown", "LightCoral", "IndianRed", "Brown", "FireBrick", "Maroon", "Red", "Tomato", "DarkSalmon", "Sienna", "Chocolate", "SandyBrown", "DarkOrange", "BurlyWood", "NavajoWhite", "BlanchedAlmond", "GoldenRod", "Gold", "Beige", "Yellow", "OliveDrab", "LawnGreen", "HoneyDew", "DarkSeaGreen", "Green", "Lime",  "Turquoise", "DarkSlateGrey", "Teal", "Cyan", "LightBlue", "SteelBlue", "RoyalBlue", "Lavender", "Navy", "Blue", "MediumPurple", "Indigo", "DarkOrchid", "Thistle", "Violet", "Purple", "Magenta", "Orchid", "DeepPink", "HotPink", "Crimson", "Pink"]
colors = ["Red"]

NumberProcessors = 20
NumHeaderLines = 0
NumLegendColumns = 4

################################

#FreqsMin = [159000, 161000, 163000, 165000, 167000, 169000, 171000, 173000, 175000, 177000, 179000, 181000, 183000, 185000, 187000, 189000, 191000, 193000, 195000, 197000, 199000, 201000, 203000, 205000, 207000, 209000]
FreqsMin = [218136]

#FreqsMax = [161000, 163000, 165000, 167000, 169000, 171000, 173000, 175000, 177000, 179000, 181000, 183000, 185000, 187000, 189000, 191000, 193000, 195000, 197000, 199000, 201000, 203000, 205000, 207000, 209000, 211000]
FreqsMax = [220000]

#
# intensity ranges in Kelvin	

IntensityMin = -1.5
				  
IntensityMax = +100.0

#
# setting up the matplotlib canvas
#
	
nx = 1
ny = len(chunks)

fig = plt.figure(figsize=(nx * 24, ny * 5), dpi=160)

#
# generate plot

for c, chunk in enumerate(chunks):
	
    pos = c + 1

    ax1 = fig.add_subplot(ny, nx, pos)

    #
    # get input	
    #
    FreqMin = FreqsMin[chunk]
    FreqMax = FreqsMax[chunk]

    #
    # XCLASS input parameters
    #
    RestFreq = 0.5 * (FreqMin + FreqMax)

    xLowerLimit = FreqMin / 1000.	
    xUpperLimit = FreqMax / 1000.

    yLowerLimit = IntensityMin
    yUpperLimit = IntensityMax

    FileName = "spectra/north_spw0_mean.fits"

    f1  = fits.open(FileName)
    w1 = WCS(f1[0].header)
    y1 = f1[0].data
    y1 = y1*145
    x1 = w1.wcs_pix2world(np.arange(f1[0].data.size), 0)[0] / 1e6
    f2 = fits.open('spectra/north_spw1_mean.fits')
    w2 = WCS(f2[0].header)
    y2 = f2[0].data
    y2 = y2*145
    x2 = w2.wcs_pix2world(np.arange(f2[0].data.size), 0)[0] / 1e6

    x = np.concatenate([x1,x2[x2>x1.max()]])
    y = np.concatenate([y1,y2[x2>x1.max()]])
    
    expdata = np.array([x,y]).T

    #
    # plot the data
    #
    dfreq = expdata[:, 0] / 1000.
    dvel = -(dfreq - RestFreq) / RestFreq / 1000 * c + vLSR	
    data = expdata[:,1]

    ax1.plot(dfreq, data, color='Gray', drawstyle = 'steps-mid', label="Data", lw=2) 

    #
    # run XCLASS if wanted
    #
	
    for m, MolfitsFileBaseName in enumerate(MolfitsFileBaseNames):
          
          MolfitsFileName = "%s.molfit" % MolfitsFileBaseName
          
          modeldata, log, TransEnergies, IntOptical, jobDir = myXCLASS()
          
          #
          # plot output
          #
          
          freq = modeldata[:,0] / 1000.
          mvel = modeldata[:,1]
          fit = modeldata[:,2]
          ax1.plot(freq, fit, color=colors[m], alpha=0.7, label=MolfitsFileBaseName, lw=2)

    ax1.set_xlim(xLowerLimit,xUpperLimit)
    ax1.set_ylim(yLowerLimit,yUpperLimit) 
    plt.xticks([i / 10. for i in range(FreqMin / 100, FreqMax / 100)])

    ax1.grid()

    #
    # setup the legend  --> have a legend for the first subplot (at the top) and the last subplot (at the bottom)

    if c == 0:
          legend = ax1.legend(bbox_to_anchor=(0., 1.05, 1., .102), loc=3, ncol=NumLegendColumns, mode="expand", borderaxespad=0.)

    elif c == (len(chunks) - 1):
          #
          # determine the number of lines in the legend
          
          if len(MolfitsFileBaseNames) % NumLegendColumns != 0:
                    numLines = (len(MolfitsFileBaseNames) / NumLegendColumns) + 1
          else:
                    numLines = len(MolfitsFileBaseNames) / NumLegendColumns
          #
          # adjust the location of the bottom legend
          print "numLines = ", numLines
          
          #legend2 = ax1.legend(bbox_to_anchor=(0., (-0.25 - (numLines / 10.)), 1., 1.), loc=3, ncol=NumLegendColumns, mode="expand", borderaxespad=0.)

    #
    # axis labels 

    ax1.set_ylabel (r"T$_{\rm MB}$ [K]")
    ax1.set_xlabel(r"Frequency [GHz]")

#fig.show()

plt.savefig("W51.png", )
plt.savefig("W51.pdf", bbox_extra_artists=(legend,legend,), bbox_inches='tight')
plt.close()

os.system("open W51.png")
