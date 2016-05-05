import os
import numpy as np
import paths
from astropy.convolution import Tophat2DKernel
from simulated_cores import simulate_grid_for_fitsfile
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from astropy import log
import radio_beam

def setup_sims(maxscale=3):
    for scale in range(1,maxscale):
        fh = simulate_grid_for_fitsfile(paths.dpath("w51_te_continuum_best.fits"),
                                        separation=45*scale, amplitude_range=[0.001,0.1],
                                        scale=scale,
                                        random_offset=20)
        fh[0].header['BUNIT'] = 'Jy/beam'
        fh.writeto(paths.simpath("simimage_scale{0}_gaussian.fits".format(scale)),
                   clobber=True)

        fh = simulate_grid_for_fitsfile(paths.dpath("w51_te_continuum_best.fits"),
                                        separation=45*scale, amplitude_range=[0.001,0.1],
                                        kernel=Tophat2DKernel,
                                        scale=scale,
                                        random_offset=20)
        fh[0].header['BUNIT'] = 'Jy/beam'
        fh.writeto(paths.simpath("simimage_scale{0}_tophat.fits".format(scale)),
                   clobber=True)

def run_sims(maxscale=3):
    for scale in range(1,maxscale):
        synthetically_image_fitsfile(paths.simpath("simimage_scale{0}_gaussian.fits".format(scale)),
                                     base_name='casa_simimage_scale{0}_gaussian'.format(scale),
                                     cleanup=True)

        synthetically_image_fitsfile(paths.simpath("simimage_scale{0}_tophat.fits".format(scale)),
                                     base_name='casa_simimage_scale{0}_tophat'.format(scale),
                                     cleanup=True)

def ridiculous_tests():
    for unit in ('Jy/beam', 'Jy/pixel', 'Jy', 'K'):
        fh = simulate_grid_for_fitsfile(paths.dpath("w51_te_continuum_best.fits"),
                                        separation=45*3, amplitude_range=[0.001,0.1],
                                        scale=3,
                                        random_offset=20)
        assert fh[0].data.max() > 0
        fh[0].header['BUNIT'] = unit
        fh.writeto(paths.simpath("stupidtest_{0}.fits".format(unit.replace("/","_"))),
                   clobber=True)
        synthetically_image_fitsfile(paths.simpath("stupidtest_{0}.fits".format(unit.replace("/","_"))),
                                     base_name="casa_stupidtest_{0}".format(unit.replace("/","_")),
                                     cleanup=True
                                    )

        
def ratio_ims(inpfn, outfn):
    din = fits.getdata(inpfn).squeeze()
    dout = fits.getdata(outfn).squeeze()
    isfin = np.isfinite(din) & np.isfinite(dout)
    ispos = (din>0) & (dout>0)
    issignal = din > 0.0001
    mask = isfin & ispos & issignal
    ratio = dout[mask].sum() / din[mask].sum()
    return ratio

def analyze_sims(maxscale=3):
    gresults = {}
    tresults = {}
    for scale in range(1,maxscale):
        inpfn = "simimage_scale{0}_gaussian.fits".format(scale)
        outfn = 'casa_simimage_scale{0}_gaussian_model_tclean_clean.image.fits'.format(scale)
        if not os.path.exists(outfn):
            outfn = 'casa_simimage_scale{0}_gaussian_model_tclean_clean.residual.fits'.format(scale)
        gresults[scale] = ratio_ims(inpfn, outfn)
        inpfn = "simimage_scale{0}_tophat.fits".format(scale)
        outfn = 'casa_simimage_scale{0}_tophat_model_tclean_clean.image.fits'.format(scale)
        if not os.path.exists(outfn):
            outfn = 'casa_simimage_scale{0}_tophat_model_tclean_clean.residual.fits'.format(scale)
        tresults[scale] = ratio_ims(inpfn, outfn)
        
    return gresults, tresults

def synthetically_image_fitsfile(fitsfilename, base_name,
                                 msfile="simulation_continuum.ms",
                                 cleanup=True,
                                 imsize=[3072,3072],
                                 cell='0.05arcsec',
                                 niter=10000,
                                 phasecenter="J2000 19:23:41.629000 +14.30.42.38000",
                                 threshold='50mJy',
                                 addnoise=True,
                                 dividebyppbeam=True,
                                ):
    """
    Assumes the image is already in the correct pixel space an units of..
    """

    assert 'fits' in fitsfilename, "Must have .fits as extension"
    ffile = fits.open(fitsfilename)
    hdr = ffile[0].header
    #assert ffile[0].header['BUNIT'] in ("Jy/beam","Jy/pixel")
    print("input FITS file BUNIT: {0}".format(ffile[0].header['BUNIT']))
    if dividebyppbeam:
        mywcs = wcs.WCS(ffile[0].header)
        beam = radio_beam.Beam.from_fits_header(ffile[0].header)
        pixel_scale = np.abs(mywcs.celestial.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
        ppbeam = (beam.sr/(pixel_scale**2)).decompose().value
        log.info("Pixels per Beam to divide by: {0}".format(ppbeam))
        assert ppbeam > 1
                

    result = sm.openfromms(msfile)
    print("Result of openfromms: {0}".format(result))
    result = sm.setvp()
    print("Result of setvp: {0}".format(result))

    freqs = {}
    for spw in range(6):
        mymsmd = msmdtool()
        mymsmd.open(msfile)

        reffreq = mymsmd.reffreq(spw)['m0']['value']
        freqs[spw] = reffreq

        mymsmd.close()


    ffile[0].header['CTYPE1'] = 'RA---SIN'
    ffile[0].header['CTYPE2'] = 'DEC--SIN'
    ffile[0].header['CRPIX3'] = 1.
    ffile[0].header['CRVAL3'] = freqs[0]-8e9
    ffile[0].header['CUNIT3'] = 'Hz'
    ffile[0].header['CDELT3'] = 1e11
    ffile[0].header['CTYPE3'] = 'FREQ'
    ffile[0].header['CRPIX4'] = 1.
    ffile[0].header['CRVAL4'] = 1.
    ffile[0].header['CUNIT4'] = ''
    ffile[0].header['CDELT4'] = 1.
    ffile[0].header['CTYPE4'] = 'STOKES'
    ffile[0].header['RESTFRQ'] = freqs[0]
    for k in ffile[0].header:
        if ("PC" == k[:2] or "PV" == k[:2] or k[:3]=="OBS" or k[:3]=='TEL'
            or k == 'DATE-OBS'):
            del ffile[0].header[k]
    ffile[0].header['TELESCOP'] = 'NotReal'
    if dividebyppbeam:
        ffile[0].data /= ppbeam
    if ffile[0].data.ndim == 3:
        ffile[0].data = ffile[0].data[None,:,:,:]
    elif ffile[0].data.ndim == 2:
        ffile[0].data = ffile[0].data[None,None,:,:]
    ffile.writeto(fitsfilename.replace(".fits","_temp.fits"), clobber=True)
    hdr = ffile[0].header

    sim_image = fitsfilename.replace(".fits",".image")
    # doesn't always work: unreliable = don't use.
    # rmresult = rmtables([sim_image])
    os.system('rm -rf {0}'.format(sim_image))
    importfits(fitsimage=fitsfilename.replace(".fits","_temp.fits"),
               imagename=sim_image,
               overwrite=True,
               defaultaxes=True,#['RA---TAN','DEC--TAN','FREQUENCY','STOKES'],
               defaultaxesvalues=['{0}deg'.format(hdr['CRVAL1']),
                                  '{0}deg'.format(hdr['CRVAL2']),
                                  '{0}GHz'.format(hdr['CRVAL3']),
                                  'I'],
               beam=["{0}deg".format(hdr['BMAJ']),
                     "{0}deg".format(hdr['BMIN']),
                     "{0}deg".format(hdr['BPA'])],
              )
    print("FITS CDELT1={0}, CDELT2={1}".format(hdr['CDELT1'], hdr['CDELT2']))
    print("image CDELT1={0[value]}{0[unit]}, CDELT2={1[value]}{1[unit]}"
          .format(imhead(imagename=sim_image, mode='get', hdkey='CDELT1'),
                  imhead(imagename=sim_image, mode='get', hdkey='CDELT2'),)
         )
    print("image CDELT3={0[value]}{0[unit]}, CDELT4={1[value]}{1[unit]}"
          .format(imhead(imagename=sim_image, mode='get', hdkey='CDELT3'),
                  imhead(imagename=sim_image, mode='get', hdkey='CDELT4'),)
         )
    ia.open(sim_image)
    print("first read-in: BUNIT: {0}".format(ia.summary()['unit']))
    ia.close()

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
    print("CTYPE1={0}, CTYPE2={1}".format(hdr['CTYPE1'], hdr['CTYPE2']))
    print("CTYPE3={0}, CTYPE4={1}".format(hdr['CTYPE3'], hdr['CTYPE4']))

    #os.system('rm -rf {0}'.format(sim_image))
    ## try re-importing (this definitely should not fix any outstanding issues)
    #importfits(fitsimage=fitsfilename,
    #           imagename=sim_image,
    #           overwrite=True,
    #           # The beam is OBVIOUSLY AND CORRECTLY in the header.
    #           # but if you leave this out, CASA complains loudly
    #           beam=["{0}deg".format(hdr['BMAJ']),
    #                 "{0}deg".format(hdr['BMIN']),
    #                 "{0}deg".format(hdr['BPA'])],
    #          )
    #ia.open(sim_image)
    #print(".image.fits BUNIT: {0}".format(hdr['BUNIT']))
    #print("second read-in: BUNIT: {0}".format(ia.summary()['unit']))
    #ia.close()

    #sm.openfromms("w51_contvis_selfcal_0.ms")
    #sm.openfromms("w51_test_onechan.ms")
    success = sm.predict(sim_image)
    print("Result of predict: {0}".format(success))
    assert success

    # TODO: get these from ASDM_CALWVR and WEATHER
    if addnoise:
        success2 = sm.setnoise(mode='tsys-atm', relhum=60.0, pwv='2mm', tatmos=265.0, )
        print("Result of setnoise: {0}".format(success2))
        success3 = sm.corrupt()
        print("Result of corrupt: {0}".format(success3))
    sm.done()
    sm.close()

    # problem:
    # plotms(vis='continuum_7m12m_noflag.ms', xaxis='uvdist', ydatacolumn='model')

    os.system('rm -rf {0}_model.ms'.format(base_name))
    assert split(vis=msfile,
                 outputvis="{0}_model.ms".format(base_name),
                 datacolumn='data')


    #os.system('rm -rf {0}_model_tclean_dirty*'.format(base_name))
    #tclean(vis='{0}_model.ms'.format(base_name),
    #       imagename='{0}_model_tclean_dirty'.format(base_name),
    #       field='',
    #       spw='',
    #       specmode='mfs',
    #       deconvolver='clark',
    #       imsize = imsize,
    #       cell= cell,
    #       weighting = 'uniform',
    #       phasecenter=phasecenter,
    #       #scales=[0,3,9,27,81],
    #       robust = -2.0,
    #       niter = 0,
    #       threshold = '1.0mJy',
    #       interactive = False,
    #       gridder = 'mosaic',
    #       savemodel='none',
    #       )
    #exportfits(imagename='{0}_model_tclean_dirty.residual'.format(base_name),
    #           fitsimage='{0}_model_tclean_dirty.image.fits'.format(base_name),
    #           dropdeg=True, overwrite=True)

    print("Cleaning")
    os.system('rm -rf {0}_model_tclean_clean*'.format(base_name))
    tclean(vis='{0}_model.ms'.format(base_name),
           imagename='{0}_model_tclean_clean'.format(base_name),
           field='',
           spw='0',
           specmode='mfs',
           #deconvolver='clark',
           imsize = imsize,
           cell=cell,
           weighting = 'briggs',
           phasecenter=phasecenter,
           #scales=[0,3,9,27,81],
           robust = -2.0,
           niter = niter,
           threshold = threshold,
           interactive = False,
           pblimit=0.4,
           gridder='mosaic',
           interpolation='linear',
           savemodel='none',
           )

    exportfits(imagename='{0}_model_tclean_clean.image'.format(base_name), fitsimage='{0}_model_tclean_clean.image.fits'.format(base_name),  dropdeg=True, overwrite=True)
    exportfits(imagename='{0}_model_tclean_clean.model'.format(base_name), fitsimage='{0}_model_tclean_clean.model.fits'.format(base_name),  dropdeg=True, overwrite=True)
    exportfits(imagename='{0}_model_tclean_clean.residual'.format(base_name), fitsimage='{0}_model_tclean_clean.residual.fits'.format(base_name),  dropdeg=True, overwrite=True)
    #assert os.path.exists('{0}_model_tclean_clean.image'.format(base_name))

    # # msclean has been failing often, so I reduced niter from 50000 to 100
    # # and increased threshold from 2.0 to 5.0 for speed
    # os.system('rm -rf {0}_model_tclean_msclean*'.format(base_name))
    # tclean(vis='{0}_model.ms'.format(base_name),
    #        imagename='{0}_model_tclean_msclean'.format(base_name),
    #        field='',
    #        spw='',
    #        specmode='mfs',
    #        deconvolver='multiscale',
    #        imsize = [1536,1536],
    #        cell= '0.1arcsec',
    #        weighting = 'uniform',
    #        phasecenter=phasecenter,
    #        scales=[0,3,9,27],
    #        robust = -2.0,
    #        niter = 100,
    #        threshold = '5.0mJy',
    #        interactive = False,
    #        gridder = 'mosaic',
    #        savemodel='none',
    #        )
    # exportfits(imagename='{0}_model_tclean_msclean.image'.format(base_name), fitsimage='{0}_model_tclean_msclean.image.fits'.format(base_name),  dropdeg=True, overwrite=True)
    # exportfits(imagename='{0}_model_tclean_msclean.model'.format(base_name), fitsimage='{0}_model_tclean_msclean.model.fits'.format(base_name),  dropdeg=True, overwrite=True)
    # exportfits(imagename='{0}_model_tclean_msclean.residual'.format(base_name), fitsimage='{0}_model_tclean_clean.residual.fits'.format(base_name),  dropdeg=True, overwrite=True)

    if cleanup:
        for filetype in ('image', 'model', 'residual', 'pb', 'psf', 'mask',
                         'sumwt', 'weight', 'image'):
            os.system('rm -rf {0}_model_tclean_msclean.{1}'.format(base_name, filetype))
            os.system('rm -rf {0}_model_tclean_clean.{1}'.format(base_name, filetype))
            os.system('rm -rf {0}_model_tclean_dirty.{1}'.format(base_name, filetype))
            os.system('rm -rf '+'{0}_model.ms'.format(base_name))
        os.system('rm {0}'.format(fitsfilename.replace(".fits","_temp.fits")))

    #fits_in = fits.open(sim_image+".fits")
    #dirty_header = fits.getheader('{0}_model_tclean_dirty.image.fits'.format(base_name))
    #header_in = wcs.WCS(fits_in[0].header).celestial.to_header()
    #header_in['NAXIS1'] = fits_in[0].data.shape[-1]
    #header_in['NAXIS2'] = fits_in[0].data.shape[-2]
    #reproj_fits_in = hcongrid.hastrom(fits_in[0].data.squeeze(),
    #                                  header_in,
    #                                  dirty_header)
    #d_in = reproj_fits_in
    #d_out_dirty = fits.getdata('{0}_model_tclean_dirty.image.fits'.format(base_name))
    #d_out_cln = fits.getdata('{0}_model_tclean_clean.image.fits'.format(base_name))
    #d_out_mscln = fits.getdata('{0}_model_tclean_msclean.image.fits'.format(base_name))

    #print("Distance from dirty: {0}".format((np.nan_to_num(d_in-d_out_dirty)**2).sum()**0.5))
    #print("Distance from cln: {0}".format((np.nan_to_num(d_in-d_out_cln)**2).sum()**0.5))
    #print("Distance from mscln: {0}".format((np.nan_to_num(d_in-d_out_mscln)**2).sum()**0.5))

def make_blank(msfile="simulation_continuum.ms", addnoise=True,
               imsize=[3072,3072], cell='0.05arcsec', niter=10000,
               phasecenter="J2000 19:23:41.629000 +14.30.42.38000",
               threshold='50mJy',):
    """
    Make a blank (noise-only) field for reference
    """

    tb.open(msfile)
    d = tb.getcol('DATA')
    d[:] = 0
    tb.putcol('DATA', d)
    tb.close()

    result = sm.openfromms(msfile)
    print("Result of openfromms: {0}".format(result))
    result = sm.setvp()
    print("Result of setvp: {0}".format(result))
    #sm.openfromms("w51_contvis_selfcal_0.ms")
    #sm.openfromms("w51_test_onechan.ms")
    #success = sm.predict(complist="")
    #print("Result of predict: {0}".format(success))
    #assert success

    # TODO: get these from ASDM_CALWVR and WEATHER
    if addnoise:
        success2 = sm.setnoise(mode='tsys-atm', relhum=60.0, pwv='2mm', tatmos=265.0, )
        print("Result of setnoise: {0}".format(success2))
        success3 = sm.corrupt()
        print("Result of corrupt: {0}".format(success3))
    sm.done()
    sm.close()

    os.system('rm -rf blank_model.ms')
    assert split(vis=msfile,
                 outputvis="blank_model.ms",
                 datacolumn='data')

    print("Cleaning")
    os.system('rm -rf blank_model_tclean_clean*')
    tclean(vis='blank_model.ms',
           imagename='blank_model_tclean_clean',
           field='',
           spw='0',
           specmode='mfs',
           #deconvolver='clark',
           imsize = imsize,
           cell=cell,
           weighting = 'briggs',
           phasecenter=phasecenter,
           #scales=[0,3,9,27,81],
           robust = -2.0,
           niter = niter,
           threshold = threshold,
           interactive = False,
           pblimit=0.4,
           gridder='mosaic',
           interpolation='linear',
           savemodel='none',
           )

    base_name = 'blank'
    exportfits(imagename='{0}_model_tclean_clean.image'.format(base_name), fitsimage='{0}_model_tclean_clean.image.fits'.format(base_name),  dropdeg=True, overwrite=True)
    exportfits(imagename='{0}_model_tclean_clean.model'.format(base_name), fitsimage='{0}_model_tclean_clean.model.fits'.format(base_name),  dropdeg=True, overwrite=True)
    exportfits(imagename='{0}_model_tclean_clean.residual'.format(base_name), fitsimage='{0}_model_tclean_clean.residual.fits'.format(base_name),  dropdeg=True, overwrite=True)
