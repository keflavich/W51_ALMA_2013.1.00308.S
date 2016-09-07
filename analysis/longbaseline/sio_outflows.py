from spectral_cube import SpectralCube
from astropy import units as u

northcube = SpectralCube.read('/Volumes/passport/alma/w51/longbaseline/W51northcax.SPW0_ALL_medsub_cutout.fits')
northvcube = northcube.with_spectral_unit(u.km/u.s, rest_value=217.10498*u.GHz,
                                          velocity_convention='radio')

northslab = northvcube.spectral_slab(-20*u.km/u.s, 180*u.km/u.s)
northmed = northslab.median(axis=0)
northmslab = northslab-northmed

northsioblue = northmslab.spectral_slab(-32*u.km/u.s, 55*u.km/u.s).moment0()
northsioblue.write('/Users/adam/work/w51/alma/FITS/longbaseline/SiO_m32to55kms_north.fits')

northsiored = northmslab.spectral_slab(74*u.km/u.s, 118*u.km/u.s).moment0()
northsiored.write('/Users/adam/work/w51/alma/FITS/longbaseline/SiO_74to118kms_north.fits')




e2cube = SpectralCube.read('/Volumes/passport/alma/w51/longbaseline/W51e2cax.SPW0_ALL_medsub_cutout.fits')
e2vcube = e2cube.with_spectral_unit(u.km/u.s, rest_value=217.10498*u.GHz,
                                    velocity_convention='radio')

e2slab = e2vcube.spectral_slab(-20*u.km/u.s, 180*u.km/u.s)
e2slab.allow_huge_operations = True
e2med = e2slab.median(axis=0)
e2mslab = e2slab-e2med

e2sioblue = e2mslab.spectral_slab(-32*u.km/u.s, 55*u.km/u.s).moment0()
e2sioblue.write('/Users/adam/work/w51/alma/FITS/longbaseline/SiO_m32to55kms_e2.fits')

e2siored = e2mslab.spectral_slab(74*u.km/u.s, 118*u.km/u.s).moment0()
e2siored.write('/Users/adam/work/w51/alma/FITS/longbaseline/SiO_74to118kms_e2.fits')
