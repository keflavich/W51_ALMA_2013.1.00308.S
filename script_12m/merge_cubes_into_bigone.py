"""
http://docs.astropy.org/en/stable/io/fits/appendix/faq.html#how-can-i-create-a-very-large-fits-file-from-scratch
"""
from astropy import log
from astropy import wcs
from astropy.io import fits
import numpy as np
import glob
import re
import os
from astropy.utils.console import ProgressBar
from spectral_cube import SpectralCube

frange = {0: [218136., 218575.],
          1: [218422., 220268.],
          2: [230436., 232310.],
          3: [233041., 234915.],
         }
fstep = {0:130., # kHz
         1:500., # kHz
         2:500., # kHz
         3:500., # kHz
        }
nchans_total = {ii: int(np.abs(np.diff(frange[ii])/fstep[ii]*1000)[0])
                for ii in frange}


# Extract the appropriate pixel indices from the file name.
# A more sophisticated approach is probably better, in which the individual
# cubes are inspected for their start/end frequencies.
# But, on the other hand, for this process to make any sense at all, you
# have to have done the original cube imaging right
def getinds(fn):
    inds = re.search('channels([0-9]*)to([0-9]*)', fn).groups()
    return [int(ii) for ii in inds]

def make_spw_cube(spw='spw{0}', spwnum=0, fntemplate='w51pointing32',
                  overwrite_existing=False, bmaj_limits=None,
                  fnsuffix="", filesuffix='image.fits',
                  slices=None,
                  cropends=False,
                  minimize=True,
                  skip_failures=False,
                  add_beam_info=True):
    """
    Parameters
    ----------
    spw : str
        String template for the input/output name
    spwnum : int
        The spectral window number
    fntemplate : str
        Filename template (goes into the glob)
    overwrite_existing : bool
        Overwrite data in the output cube?
    cropends: bool or int
        Number of pixels to crop off the ends of an image
    minimize: bool
        Compute the spatial minimal subcube before building the cube?  Slices
        for all subsequent cubes will be computed from the first cube.
    """
    spw = spw.format(spwnum)

    big_filename = '{1}_{0}{2}_lines.fits'.format(spw, fntemplate, fnsuffix)

    header_fn = glob.glob('piece_of_{1}_cube{2}.{0}.channels0to*.{3}'.format(spw, fntemplate, fnsuffix, filesuffix))
    if len(header_fn) != 1:
        raise ValueError("Found too many or too few matches: {0}".format(header_fn))
    else:
        header_fn = header_fn[0]

    # First set up an empty file
    if not os.path.exists(big_filename):

        log.info("Creating new large cube from {0}".format(header_fn))

        if minimize:
            cube0 = SpectralCube.read(header_fn)
            if slices is None:
                slices = cube0.subcube_slices_from_mask(cube0.mask,
                                                        spatial_only=True)
            # use the calculated 3rd dimension, plus the difference of the
            # x and y slices
            #header['NAXIS2'] = slices[1].stop-slices[1].start
            #header['NAXIS1'] = slices[2].stop-slices[2].start
            header = cube0[slices].header
        else:
            header = fits.getheader(header_fn)

        # Make an arbitrary, small data before prepping the header
        data = np.zeros((100, 100), dtype=np.float32)
        hdu = fits.PrimaryHDU(data=data, header=header)
        cdelt_sign = np.sign(hdu.header['CDELT3'])
        # Set the appropriate output size (this can be extracted from the LISTOBS)
        header['NAXIS3'] = nchans_total[spwnum]
        if cdelt_sign == -1:
            ind0, ind1 = getinds(header_fn)
            # if ind1 == cubeshape, need to add 1 more pixel
            # but, if cropping is intended & cube is 1 pixel bigger...
            # (does this make sense?!)
            extra = ind1 == cube0.shape[0]
            header['CRPIX3'] = nchans_total[spwnum] - ind1 + extra

        shape = (header['NAXIS3'], header['NAXIS2'], header['NAXIS1'])

        # Write to disk
        header.tofile(big_filename)
        # Using the 'append' io method, update the *header*
        with open(big_filename, 'rb+') as fobj:
            # Seek past the length of the header, plus the length of the
            # data we want to write.
            # The -1 is to account for the final byte that we are about to
            # write:
            # 'seek' works on bytes, so divide #bits / (bytes/bit)
            fobj.seek(len(header.tostring()) + (shape[0] *
                                                shape[1] *
                                                shape[2] *
                                                int(np.abs(header['BITPIX'])/8)) -
                      1)
            fobj.write(b'\0')

        big_cube = SpectralCube.read(big_filename)
        header_cube = SpectralCube.read(header_fn)
        # in both cases, SpectralCube sorts the extrema
        assert big_cube.spectral_extrema[0] == header_cube.spectral_extrema[0]
        assert np.all(big_cube.wcs.wcs.cdelt == header_cube.wcs.wcs.cdelt)


    # Find the appropriate files (this is NOT a good way to do this!  Better to
    # provide a list.  But wildcards are quick & easy...
    files = glob.glob("piece_of_{1}_cube{2}.{0}.chan*{3}"
                      .format(spw,fntemplate,fnsuffix,filesuffix))
    log.info("Files to be merged: ")
    log.info(str(files))

    # open the file in update mode (it should have the right dims now)
    hdul = fits.open(big_filename, mode='update')

    if add_beam_info:
        shape = hdul[0].data.shape[0]
        if len(hdul) > 1 and isinstance(hdul[1], fits.BinTableHDU):
            pass
        else:
            hdul.append(fits.BinTableHDU(np.recarray(shape,
                                                     names=['BMAJ','BMIN','BPA','CHAN','POL'],
                                                     formats=['f4','f4','f4','i4','i4'])))

    bad = {}

    for fn in ProgressBar(files):
        hdr = fits.getheader(fn)
        mwcs = wcs.WCS(hdr)
        (f0,f1), = mwcs.sub([wcs.WCSSUB_SPECTRAL]).wcs_pix2world((1,hdr['NAXIS3']), 1)
        ind0,ind1 = getinds(fn)
        log.info("{0}->{1} {2}".format([ind0,ind1], [f0,f1], fn))

        if ind0 > 0:
            if skip_failures and ind1-ind0 != fits.getdata(fn).shape[0]:
                print("{0} has wrong size".format(fn))
                bad[fn] = 'ind1-ind0={0}, shape0={1}'.format(ind1-ind0,
                                                             fits.getdata(fn).shape[0])
                continue
            else:
                # this might not be exactly right... but I think it is.
                assert ind1-ind0 == fits.getdata(fn).shape[0], "Data cube has wrong size."
        else:
            if skip_failures and ind1 > fits.getdata(fn).shape[0]:
                print("{0} has wrong spectral shape".format(fn))
                bad[fn] = "ind1 {0} > shape {1}".format(ind1,
                                                        fits.getdata(fn).shape[0])
                continue
            else:
                # only worry about the 2nd index being in range
                assert ind1 <= fits.getdata(fn).shape[0]

        if slices is None and minimize:
            # this is wrong, but I'm not ready to delete it yet
            #cube0 = SpectralCube.read(fn)
            #slices = cube0.subcube_slices_from_mask(cube0.mask,
            #                                        spatial_only=True)

            cube0 = SpectralCube.read(header_fn)
            slices = cube0.subcube_slices_from_mask(cube0.mask,
                                                    spatial_only=True)
        elif slices is None:
            slices = (slice(None),)*3

        if cropends:
            # don't crop 1st or last pixel in full cube
            if ind0 > 0:
                ind0 = ind0 + cropends
                dataind0 = cropends
                extra = 0
            else:
                # because I forgot to reduce nchan, there is an "extra" pixel
                # when we start at zero (there should not be a corresponding one
                # when we end too late)
                dataind0 = 0
                extra = 1

            if ind1 < nchans_total[spwnum] - 1:
                ind1 = ind1 - cropends
                dataind1 = - cropends - extra
            else:
                dataind1 = None

        if 'cdelt_sign' not in locals():
            cdelt_sign = np.sign(fits.getheader(fn)['CDELT3'])
            log.warn("cdelt_sign was not defined: overwriting a"
                     " previously-existing file.  "
                     "This may not be what you want; the data could be going "
                     "opposite the parent cube.  Check that the original "
                     "header is OK. sign(CDELT) is now {0}, "
                     "while for the big header it is {1}"
                     .format(cdelt_sign,
                             np.sign(fits.getheader(big_filename)['CDELT3'])))
        if cdelt_sign == -1:
            ind1, ind0 = (nchans_total[spwnum] - ind0 - 1,
                          nchans_total[spwnum] - ind1 - 1)
        plane = hdul[0].data[ind0]
        if np.all(plane == 0) or overwrite_existing:
            log.info("Replacing indices {0}->{2} {1}"
                     .format(getinds(fn), fn, (ind0,ind1)))

            data = fits.getdata(fn)

            if bmaj_limits is not None:
                beamtable = fits.open(fn)[1]
                ok_beam = ((beamtable.data['BMAJ'] > bmaj_limits[0]) &
                           (beamtable.data['BMAJ'] < bmaj_limits[1]))
                data[~ok_beam] = np.nan
            if add_beam_info:
                beamtable = fits.open(fn)[1]
                hdul[1].data[ind0:ind1] = beamtable.data[dataind0:dataind1]

            if data[dataind0:dataind1, slices[1], slices[2]].shape != hdul[0].data[ind0:ind1,:,:].shape:
                bad[fn] = "bigshape {0}, inpshape {1}".format(data[dataind0:dataind1, slices[1], slices[2]].shape,
                                                              hdul[0].data[ind0:ind1,:,:].shape)
                if skip_failures:
                    print("{0} has wrong shape".format(fn))
                    continue
                else:
                    raise ValueError("{0} has wrong shape".format(fn))


            hdul[0].data[ind0:ind1,:,:] = data[dataind0:dataind1, slices[1], slices[2]]
            hdul.flush()

    if skip_failures:
        return bad
