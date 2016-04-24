"""
Generate simulated images to synthetically observe to determine completeness in
various radial bins etc.
"""
import numpy as np
import os
from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
from astropy import wcs
from astropy import units as u
import radio_beam
from FITS_tools import hcongrid

def gridded_integrals(r_core, alpha, gridsize=100, plummer=False):

    zz,yy,xx = np.indices([gridsize]*3, dtype='float')
    center = gridsize/2.
    rr = np.sum([(ii-center)**2 for ii in (xx,yy,zz)], axis=0)**0.5

    if plummer:
        dens = (1+rr**2/r_core**2)**-2.5
    else:
        dens = (rr >= r_core)*(rr/r_core)**-alpha
        dens[(rr < r_core)] = 1.0

    img = dens.sum(axis=0)

    return img

def make_gaussian_sim_grid(shape, g_size, separation, amplitude_range,
                           random_offset=0):
    
    blank_im = np.zeros(shape)

    n_y = int(np.floor(shape[0] / separation))
    n_x = int(np.floor(shape[1] / separation))

    model = Gaussian2DKernel(g_size).array

    ampl_scale = amplitude_range[1]-amplitude_range[0]

    for xx in range(n_x):
        for yy in range(n_y):

            amplitude = np.random.rand()*ampl_scale + amplitude_range[0]

            xcen = (xx+0.5) * separation
            ycen = (yy+0.5) * separation
            halfsz = model.shape[0]/2

            if random_offset > 0:
                # to avoid gridding artifacts
                xcen += (np.random.rand()-0.5) * random_offset
                ycen += (np.random.rand()-0.5) * random_offset

            view = [slice(ycen-halfsz,ycen+halfsz),
                    slice(xcen-halfsz,xcen+halfsz)]

            if blank_im[view].shape != model.shape:
                #print("Skipping {0},{1}".format(xcen, ycen))
                continue

            blank_im[view] += model * amplitude

    return blank_im
    

def simulate_grid_for_fitsfile(fitsfilename, gridmaker=make_gaussian_sim_grid,
                               **kwargs):

    fh = fits.open(fitsfilename)
    mywcs = wcs.WCS(fh[0].header)
    pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5

    bm = radio_beam.Beam.from_fits_header(fh[0].header)

    base_size = bm.major.to(u.deg).value / pixscale

    simimage = gridmaker(shape=fh[0].data.shape, g_size=base_size, **kwargs)

    fh[0].data = simimage

    return fh

def synthetically_image_fitsfile(fitsfilename, base_name,
                                 msfile="w51_contvis_selfcal_3.ms"):
    """
    Assumes the image is already in the correct pixel space an units of..
    """

    assert 'fits' in fitsfilename, "Must have .fits as extension"
    ffile = fits.open(fitsfilename)
    hdr = ffile[0].header
    assert ffile[0].header['BUNIT'] == 'Jy/beam', "Must have Jy/beam units"

    sim_image = fitsfilename.replace(".fits",".image")
    # doesn't always work: unreliable = don't use.
    # rmresult = rmtables([sim_image])
    os.system('rm -rf {0}'.format(sim_image))
    importfits(fitsimage=fitsfilename,
               imagename=sim_image,
               overwrite=True,
               #defaultaxes=True,#['RA---TAN','DEC--TAN','FREQUENCY','STOKES'],
               #defaultaxesvalues=['2.909234541667E+02deg',
               #                   '1.451177222222E+01deg',
               #                   '233.9468GHz','I'],
               ## 18" = 1.22 lambda/D
               #beam=["{0}deg".format(18/3600.*distance_scaling),
               #      "{0}deg".format(18/3600.*distance_scaling),
               #      "0deg"],
              )
    print("FITS CDELT1={0}, CDELT2={1}".format(hdr['CDELT1'], hdr['CDELT2']))
    print("image CDELT1={0[value]}{0[unit]}, CDELT2={1[value]}{1[unit]}"
          .format(imhead(imagename=sim_image, mode='get', hdkey='CDELT1'),
                  imhead(imagename=sim_image, mode='get',
                         hdkey='CDELT2'),)
         )

    #imhead(sim_image, mode='put', hdkey='CDELT1', hdvalue={'value':hdr['CDELT1'], 'unit':'deg'})
    #imhead(sim_image, mode='put', hdkey='CDELT2', hdvalue={'value':hdr['CDELT2'], 'unit':'deg'})
    #imhead(sim_image, mode='put', hdkey='CRVAL1', hdvalue={'value':hdr['CRVAL1'], 'unit':'deg'})
    #imhead(sim_image, mode='put', hdkey='CRVAL2', hdvalue={'value':hdr['CRVAL2'], 'unit':'deg'})
    #imhead(sim_image, mode='put', hdkey='CRPIX1', hdvalue=hdr['CRPIX1'])
    #imhead(sim_image, mode='put', hdkey='CRPIX2', hdvalue=hdr['CRPIX2'])
    #imhead(sim_image, mode='put', hdkey='CTYPE3', hdvalue='FREQ')
    # can convert units and use this
    #imhead(sim_image, mode='put', hdkey='BUNIT', hdvalue='Jy/beam')
    #imhead(sim_image, mode='put', hdkey='BUNIT', hdvalue='Jy/pixel') #WRONG!!!
    exportfits(sim_image, sim_image+".fits",
               overwrite=True)
    hdr = fits.getheader(sim_image+".fits")
    print("CDELT1={0}, CDELT2={1}".format(hdr['CDELT1'], hdr['CDELT2']))

    os.system('rm -rf {0}'.format(sim_image))
    # try re-importing (this definitely should not fix any outstanding issues)
    importfits(fitsimage=fitsfilename,
               imagename=sim_image,
               overwrite=True,
               # The beam is OBVIOUSLY AND CORRECTLY in the header.
               # but if you leave this out, CASA complains loudly
               #beam=["{0}deg".format(18/3600.*distance_scaling),
               #      "{0}deg".format(18/3600.*distance_scaling),
               #      "0deg"],
              )
    ia.open(sim_image)
    print("BUNIT: ",ia.summary()['unit'])
    ia.close()

    #sm.openfromms("w51_contvis_selfcal_0.ms")
    sm.openfromms(msfile)
    #sm.openfromms("w51_test_onechan.ms")
    sm.setvp()
    success = sm.predict(sim_image)
    # TODO: get these from ASDM_CALWVR and WEATHER
    success2 = sm.setnoise(mode='tsys-atm', relhum=60.0, pwv='2mm', tatmos=265.0, )
    success3 = sm.corrupt()
    sm.done()
    sm.close()

    # problem:
    # plotms(vis='continuum_7m12m_noflag.ms', xaxis='uvdist', ydatacolumn='model')

    os.system('rm -rf {0}_model.ms'.format(base_name))
    assert split(vis=msfile,
                 outputvis="{0}_model.ms".format(base_name),
                 datacolumn='data')

    phasecenter = 'J2000 19h23m41.580 +14d30m41.37'

    os.system('rm -rf {0}_model_tclean_dirty*'.format(base_name))
    tclean(vis='{0}_model.ms'.format(base_name),
           imagename='{0}_model_tclean_dirty'.format(base_name),
           field='',
           spw='',
           specmode='mfs',
           deconvolver='clark',
           imsize = [1536,1536],
           cell= '0.1arcsec',
           weighting = 'uniform',
           phasecenter=phasecenter,
           #scales=[0,3,9,27,81],
           robust = -2.0,
           niter = 0,
           threshold = '1.0mJy',
           interactive = False,
           gridder = 'mosaic',
           savemodel='none',
           )
    exportfits(imagename='{0}_model_tclean_dirty.residual'.format(base_name),
               fitsimage='{0}_model_tclean_dirty.image.fits'.format(base_name),
               dropdeg=True, overwrite=True)

    os.system('rm -rf {0}_model_tclean_clean*'.format(base_name))
    tclean(vis='{0}_model.ms'.format(base_name),
           imagename='{0}_model_tclean_clean'.format(base_name),
           field='',
           spw='',
           specmode='mfs',
           deconvolver='clark',
           imsize = [1536,1536],
           cell= '0.1arcsec',
           weighting = 'uniform',
           phasecenter=phasecenter,
           #scales=[0,3,9,27,81],
           robust = -2.0,
           niter = 50000,
           threshold = '2.0mJy',
           interactive = False,
           gridder = 'mosaic',
           savemodel='none',
           )
    exportfits(imagename='{0}_model_tclean_clean.image'.format(base_name), fitsimage='{0}_model_tclean_clean.image.fits'.format(base_name),  dropdeg=True, overwrite=True)
    exportfits(imagename='{0}_model_tclean_clean.model'.format(base_name), fitsimage='{0}_model_tclean_clean.model.fits'.format(base_name),  dropdeg=True, overwrite=True)
    exportfits(imagename='{0}_model_tclean_clean.residual'.format(base_name), fitsimage='{0}_model_tclean_clean.residual.fits'.format(base_name),  dropdeg=True, overwrite=True)

    # msclean has been failing often, so I reduced niter from 50000 to 100
    # and increased threshold from 2.0 to 5.0 for speed
    os.system('rm -rf {0}_model_tclean_msclean*'.format(base_name))
    tclean(vis='{0}_model.ms'.format(base_name),
           imagename='{0}_model_tclean_msclean'.format(base_name),
           field='',
           spw='',
           specmode='mfs',
           deconvolver='multiscale',
           imsize = [1536,1536],
           cell= '0.1arcsec',
           weighting = 'uniform',
           phasecenter=phasecenter,
           scales=[0,3,9,27],
           robust = -2.0,
           niter = 100,
           threshold = '5.0mJy',
           interactive = False,
           gridder = 'mosaic',
           savemodel='none',
           )
    exportfits(imagename='{0}_model_tclean_msclean.image'.format(base_name), fitsimage='{0}_model_tclean_msclean.image.fits'.format(base_name),  dropdeg=True, overwrite=True)
    exportfits(imagename='{0}_model_tclean_msclean.model'.format(base_name), fitsimage='{0}_model_tclean_msclean.model.fits'.format(base_name),  dropdeg=True, overwrite=True)
    exportfits(imagename='{0}_model_tclean_clean.residual'.format(base_name), fitsimage='{0}_model_tclean_clean.residual.fits'.format(base_name),  dropdeg=True, overwrite=True)

    fits_in = fits.open(sim_image+".fits")
    dirty_header = fits.getheader('{0}_model_tclean_dirty.image.fits'.format(base_name))
    header_in = wcs.WCS(fits_in[0].header).celestial.to_header()
    header_in['NAXIS1'] = fits_in[0].data.shape[-1]
    header_in['NAXIS2'] = fits_in[0].data.shape[-2]
    reproj_fits_in = hcongrid.hastrom(fits_in[0].data.squeeze(),
                                      header_in,
                                      dirty_header)
    d_in = reproj_fits_in
    d_out_dirty = fits.getdata('{0}_model_tclean_dirty.image.fits'.format(base_name))
    d_out_cln = fits.getdata('{0}_model_tclean_clean.image.fits'.format(base_name))
    d_out_mscln = fits.getdata('{0}_model_tclean_msclean.image.fits'.format(base_name))

    print("Distance from dirty: {0}".format((np.nan_to_num(d_in-d_out_dirty)**2).sum()**0.5))
    print("Distance from cln: {0}".format((np.nan_to_num(d_in-d_out_cln)**2).sum()**0.5))
    print("Distance from mscln: {0}".format((np.nan_to_num(d_in-d_out_mscln)**2).sum()**0.5))
