"""
Generate simulated images to synthetically observe to determine completeness in
various radial bins etc.
"""
import numpy as np
from astropy import log
from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
from astropy import wcs
from astropy import units as u
import radio_beam

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

def make_sim_grid(shape, g_size, separation, amplitude_range,
                  kernel=Gaussian2DKernel, random_offset=0):
    
    blank_im = np.zeros(shape)

    n_y = int(np.floor(shape[0] / separation))
    n_x = int(np.floor(shape[1] / separation))

    model = kernel(g_size).array
    model /= model.max()

    ampl_scale = amplitude_range[1]-amplitude_range[0]

    total_n = n_x*n_y
    skipct = 0

    for xx in range(n_x):
        for yy in range(n_y):

            amplitude = np.random.rand()*ampl_scale + amplitude_range[0]

            xcen = (xx+0.5) * separation
            ycen = (yy+0.5) * separation
            halfsz = model.shape[0]/2.

            if random_offset > 0:
                # to avoid gridding artifacts
                xcen += (np.random.rand()-0.5) * random_offset
                ycen += (np.random.rand()-0.5) * random_offset

            view = [slice(ycen-halfsz,ycen+halfsz),
                    slice(xcen-halfsz,xcen+halfsz)]

            if blank_im[view].shape != model.shape:
                #print("Skipping {0},{1}".format(xcen, ycen))
                skipct += 1
                continue

            blank_im[view] += model * amplitude

    log.info("Skipped {0} models".format(skipct))
    if (skipct/float(total_n)) > 0.2:
        raise ValueError("Skipped more than 20% of models")

    return blank_im


def simulate_grid_for_fitsfile(fitsfilename, gridmaker=make_sim_grid,
                               scale=1, kernel=Gaussian2DKernel,
                               **kwargs):

    fh = fits.open(fitsfilename)
    mywcs = wcs.WCS(fh[0].header)
    pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5

    bm = radio_beam.Beam.from_fits_header(fh[0].header)

    base_size = bm.major.to(u.deg).value / pixscale

    simimage = gridmaker(shape=fh[0].data.shape, g_size=base_size*scale,
                         kernel=kernel,
                         **kwargs)

    fh[0].data = simimage

    return fh
