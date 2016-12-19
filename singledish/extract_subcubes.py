import os
from astropy import units as u
from spectral_cube import SpectralCube
from lines_to_extract import line_to_image_list

inpath = '/Volumes/seagate_passport/w51-apex/processed/merge'
outpath = '/Volumes/seagate_passport/w51-apex/processed/merge/cutouts'

fullcubenames = [
                 "W51_291GHz_merge.fits", "W51_293GHz_merge.fits",
                 "W51_218GHz_merge.fits", "W51_232GHz_merge.fits",
                 "W51_217GHz_merge.fits", "W51_12CO_merge.fits",
                ]

for fcn in fullcubenames:
    cube = SpectralCube.read(os.path.join(inpath, fcn))

    for linename, freqstr in line_to_image_list:
        freq = u.Quantity(freqstr)

        if freq > cube.spectral_axis.min() and freq < cube.spectral_axis.max():
            slab = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio',
                                           rest_value=freq).spectral_slab(30*u.km/u.s,
                                                                          90*u.km/u.s)
            outfn = os.path.join(outpath,
                                 "W51_{0}_{1}GHz.fits".format(linename,
                                                              freq.to(u.GHz).value)
                                )
            if not os.path.exists(outfn):
                slab.write(outfn)
            else:
                print("{0} already exists".format(outfn))
