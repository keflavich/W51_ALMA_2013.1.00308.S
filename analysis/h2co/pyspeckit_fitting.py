"""
Copied from apex_cmz_h2co then cropped to just get the grids
"""
from astropy import log
from astropy.io import fits
from pyspeckit.spectrum import models
from pyspeckit.spectrum.models.model import SpectralModel
from paths import gpath

try:
    # create the Formaldehyde Radex fitter
    # This step cannot be easily generalized: the user needs to read in their own grids
    texgrid303 = fits.getdata(gpath('fjdu_pH2CO_303_tex_1kms.fits'))
    taugrid303 = fits.getdata(gpath('fjdu_pH2CO_303_tau_1kms.fits'))
    texgrid321 = fits.getdata(gpath('fjdu_pH2CO_321_tex_1kms.fits'))
    taugrid321 = fits.getdata(gpath('fjdu_pH2CO_321_tau_1kms.fits'))
    texgrid322 = fits.getdata(gpath('fjdu_pH2CO_322_tex_1kms.fits'))
    taugrid322 = fits.getdata(gpath('fjdu_pH2CO_322_tau_1kms.fits'))
    hdr = hdrb = fits.getheader(gpath('fjdu_pH2CO_303_tex_1kms.fits'))

    # # this deserves a lot of explanation:
    # # models.formaldehyde.formaldehyde_radex is the MODEL that we are going to fit
    # # models.model.SpectralModel is a wrapper to deal with parinfo, multiple peaks,
    # # and annotations
    # # all of the parameters after the first are passed to the model function 

    h2co_radex_fitter = SpectralModel(models.formaldehyde_mm.formaldehyde_mm_radex,
                                      5,
                                      parnames=['temperature','column','density','center','width'],
                                      parvalues=[50,12,4.5,0,1],
                                      parlimited=[(True,True), (True,True),
                                                  (True,True), (False,False),
                                                  (True,False)],
                                      parlimits=[(5,300), (11,17), (3,7), (0,0), (0,0)],
                                      parsteps=[0.01,0.01,0.1,0,0], fitunits='Hz',
                                      texgrid=((218.1,218.3,texgrid303),
                                               (218.35,218.55,texgrid322),
                                               (218.6,218.85,texgrid321)),
                                      # specify the frequency range over which the grid is valid (in GHz)
                                      taugrid=((218.1,218.3,taugrid303),
                                               (218.35,218.55,taugrid322),
                                               (218.6,218.85,taugrid321)),
                                      hdr=hdrb,
                                      shortvarnames=("T","N","n","v","\\sigma"),
                                      # specify the parameter names (TeX is OK)
                                      grid_vwidth=5.0,
                                      )
except IOError as ex:
    log.exception("Could not read files from disk: cannot load H2CO RADEX fitter")
    log.exception(ex)
