from spectral_cube import SpectralCube
from astropy import units as u


for line, freq in (('SiO',217.10498*u.GHz),
                   ('HC3N', 218.32472*u.GHz),
                   ('H2CO303', 218.22219*u.GHz),
                  ):

    e8cube = SpectralCube.read('/Volumes/passport/alma/w51/longbaseline/W51e8cax.SPW0_ALL_medsub_cutout.fits')
    e8vcube = e8cube.with_spectral_unit(u.km/u.s, rest_value=freq,
                                        velocity_convention='radio')

    e8slab = e8vcube.spectral_slab(-100*u.km/u.s, 210*u.km/u.s)

    e8slab.write('/Volumes/passport/alma/w51/longbaseline/W51e8cax.{0}cutout.fits'
                 .format(line))
