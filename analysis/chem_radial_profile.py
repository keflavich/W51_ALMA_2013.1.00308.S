"""
Given the chemical maps made from chem_images, perform radial profile calculations
"""
import paths
import pylab as pl
import os
import itertools
from astropy.io import fits
import radio_beam
import numpy as np
from astropy import units as u
import pyregion
import image_tools
from astropy import wcs
from astropy.nddata import Cutout2D

import re
import glob

region_names = {'e2': 'e2_exclude_e2w.reg'}

for regfn,region,fignum,imtype,suffix in (
    ('e2_exclude_e2w.reg','e2',1,'m0',''),
    ('e2_exclude_e2w.reg','e2',2,'m0','_merge'),
    ('e2_exclude_e2w.reg','e2',3,'m0','_merge_natural'),
    ('e2_exclude_e2w.reg','e2',1,'max',''),
    ('e2_exclude_e2w.reg','e2',2,'max','_merge'),
    ('e2_exclude_e2w.reg','e2',3,'max','_merge_natural'),
   ):

    reg = pyregion.open(paths.rpath(regfn))

    fig = pl.figure(fignum)
    fig.clf()
    ax = pl.gca()
    ax.set_title(region+suffix+" "+imtype)
    path_template = paths.dpath("chemslices/chemical_{2}_slabs*_{0}*{1}.fits"
                                .format(region, suffix, imtype))

    files = glob.glob(path_template)

    linestyles = {name: itertools.cycle(['-'] + ['--'] + [':'] + ['-.'])
                  for name in region_names}

    linere = re.compile("chemical_{0}_slabs_[^_]*_(.*?)(_merge.fits|.fits)"
                        .format(imtype))


    for fn in files:
        if 'merge' in fn and 'merge' not in suffix:
            # this is a way to skip _merge when looking for chemname.fits
            continue
        fh = fits.open(fn)

        linestyle = next(linestyles[region])

        if 'BMAJ' not in fh[0].header:
            print("File {0} does not have BMAJ".format(fn))
            continue
        try:
            beam = radio_beam.Beam.from_fits_header(fh[0].header)
        except KeyError:
            print("File {0} doesn't have beam info in the header".format(fn))
            continue

        mywcs = wcs.WCS(fh[0].header)
        pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5
        ppbeam = (beam.sr/(pixscale**2*u.deg**2)).decompose().value / u.beam
        print("fn  {0} ppbeam={1:0.2f}".format(fn, ppbeam))

        data = fh[0].data
        mask = reg.get_mask(fh[0]) & np.isfinite(data)

        center = mywcs.wcs_world2pix(reg[0].coord_list[0], reg[0].coord_list[1], 0)

        nr, bins, rprof = image_tools.radialprofile.azimuthalAverage(data,
                                                                     mask=mask,
                                                                     center=center,
                                                                     binsize=1.0,
                                                                     return_nr=True)

        species = linere.search(fn).groups()[0]

        ax.plot(bins*pixscale*3600., rprof/ppbeam,
                label=species, linestyle=linestyle)
        ax.set_ylabel("Azimuthally Averaged Flux (Jy)")
        ax.set_xlabel("Radius (arcsec)")

        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.savefig(paths.fpath("chemslices/radialprofile_{2}_{0}{1}.png"
                            .format(region, suffix, imtype)),
                bbox_inches='tight')
