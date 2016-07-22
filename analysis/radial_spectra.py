import paths
from astropy import coordinates
import pylab as pl
import spectral_cube
from spectral_cube import SpectralCube
import os
from scipy import ndimage
import itertools
from astropy.io import fits
import radio_beam
import numpy as np
from astropy import units as u
import pyregion
import image_tools
from astropy import wcs

import re
import glob


def extract_radial_spectrum(cube, coordinate, excludemask, radial_bins=[(0,1),(1,2)]):


    yy,xx = np.indices(cube.shape[1:])
    rr = ((yy-coordinate[1])**2 + (xx-coordinate[0])**2)**0.5

    spectra = {}
    for inner_bin_radius,outer_bin_radius in radial_bins:
        radial_mask = (rr > inner_bin_radius) & (rr < outer_bin_radius)
        radial_mask &= ~excludemask

        avspec = cube.with_mask(radial_mask).mean(axis=(1,2))

        spectra[(inner_bin_radius,outer_bin_radius)] = avspec

    return spectra

def spectra_from_cubefn(cubefn, reg, bins_arcsec, coordinate):
    cube = SpectralCube.read(cubefn)

    pixcoordinate = cube.wcs.celestial.wcs_world2pix(coordinate.ra.deg,
                                                     coordinate.dec.deg,
                                                     0)

    pixscale = (cube.wcs.celestial.pixel_scale_matrix.diagonal()**2).sum()**0.5

    includemask = reg.get_mask(header=cube.wcs.celestial.to_header(),
                               shape=cube.shape[1:])


    spectra = extract_radial_spectrum(cube, pixcoordinate, ~includemask,
                                      radial_bins=bins_arcsec/(pixscale*3600))

    return spectra

if __name__ == "__main__":


    coordinate = coordinates.SkyCoord("19:23:43.961",
                                      "+14:30:34.56",
                                      frame='fk5',
                                      unit=(u.hour, u.deg))
    bins_ends_arcsec = np.linspace(0,2.25,7)
    bins_arcsec = np.array(list(zip(bins_ends_arcsec[:-1], bins_ends_arcsec[1:])))
    reg = pyregion.open(paths.rpath('e2_exclude_e2w.reg'))


    for spw in (0,1,2,3):
        for suffix in ('','_hires'):

            cubefn = paths.dpath('merge/fullcube_cutouts/e2cutout_full_W51_7m12m_spw{0}{1}_lines.fits'
                                .format(spw,suffix))
            print(cubefn)

            spectra = spectra_from_cubefn(cubefn, reg, bins_arcsec, coordinate)

            for bins, (key, spectrum) in zip(bins_arcsec, spectra.items()):
                include = np.isfinite(spectrum) & np.array([(bm.major < 1*u.arcsec) &
                                                            (bm.minor < 1*u.arcsec)
                                                            for bm in spectrum.beams])
                avg_beam = spectral_cube.cube_utils.average_beams(spectrum.beams,
                                                                  includemask=include)
                spectrum.meta['beam'] = avg_beam
                spectrum.write(paths.merge_spath('e2e_radial_bin_{0:0.2f}to{1:0.2f}_7m12m_spw{2}{3}.fits'
                                                 .format(bins[0], bins[1], spw,
                                                         suffix)),
                               overwrite=True
                              )

    for spw in (0,1,2,3):

        cubefn = paths.dpath('12m/fullcube_cutouts/e2cutout_full_W51_spw{0}_lines.fits'
                             .format(spw))
        print(cubefn)

        spectra = spectra_from_cubefn(cubefn, reg, bins_arcsec, coordinate)

        pl.figure(1).clf()

        for bins, (key, spectrum) in zip(bins_arcsec, spectra.items()):
            include = np.isfinite(spectrum) & np.array([(bm.major < 1*u.arcsec) &
                                                        (bm.minor < 1*u.arcsec)
                                                        for bm in spectrum.beams])
            avg_beam = spectral_cube.cube_utils.average_beams(spectrum.beams,
                                                              includemask=include)
            spectrum.meta['beam'] = avg_beam
            spectrum.write(paths.spath('e2e_radial_bin_{0:0.2f}to{1:0.2f}_spw{2}.fits'
                                       .format(bins[0], bins[1], spw)),
                           overwrite=True
                          )


            pl.plot(spectrum.spectral_axis.to(u.GHz).value, spectrum.value,
                    label='{0:0.1f}-{1:0.1f}'.format(*bins))

        pl.xlabel("Frequency (GHz)")
        pl.ylabel("Intensity (Jy)")
        pl.gca().ticklabel_format(useOffset=False)
        pl.gca().get_xaxis().get_major_formatter().set_scientific(False)
        pl.gca().get_yaxis().get_major_formatter().set_scientific(False)
        
        pl.legend(loc='best')
        pl.savefig(paths.fpath('radial_spectra/e2e_radial_spectra_spw{0}.png'
                               .format(spw)))
