"""
Make moment maps of cutouts

(this was the first version using early data; moments_rd2 (round 2) is more up
to date and probably better)
"""
import re
from astropy import units as u
from astropy import coordinates
from spectral_cube import SpectralCube
import pyregion
import paths
import pylab as pl

repl = re.compile("W51[en]2?cax")

for ii in pl.get_fignums():
    pl.close(ii)

for source,cubefn in [('e2', "W51e2cax.CH3CN_K3_nat.image.fits"),
                      ('e2', "W51e2cax.CH3CN_K3_nat_all.image.fits"),
                      ('e2', "W51e2cax.CH3CN_K8.image.pbcor.fits"),
                      ('e2', "W51e2cax.H30alpha.image.pbcor.fits"),
                      ('e8', "W51e2cax.CH3CN_K3_nat.image.fits"),
                      ('e8', "W51e2cax.CH3CN_K3_nat_all.image.fits"),
                      ('e8', "W51e2cax.CH3CN_K8.image.pbcor.fits"),
                      ('e8', "W51e2cax.H30alpha.image.pbcor.fits"),
                      ('north', "W51ncax.H30alpha.image.pbcor.fits"),
                      ('north', "W51ncax.CH3CN_K8.image.pbcor.fits"),
                      ]:
    suffix = ".image.fits" if ".image.fits" in cubefn else ".image.pbcor.fits"

    outfn = repl.sub("W51{0}cax".format(source), cubefn).replace(suffix,"_medsub_cutout.fits")

    cube = SpectralCube.read(outfn)
    print(outfn, cube)

    if 'K3_nat_all' in cubefn:
        linefn = 'CH3CN_K3_nat_all'
        rest_values = [('CH3CN_K7', 220.53932*u.GHz),
                       ('CH3CN_K6', 220.59442*u.GHz),
                       ('CH3CN_K5', 220.64108*u.GHz),
                       ('CH3CN_K4', 220.67929*u.GHz),]
    elif 'K3_nat' in cubefn:
        linefn = 'CH3CN_K3_nat'
        rest_values = [('CH3CN_K3',220.70902*u.GHz)]
    elif 'K8' in cubefn:
        linefn = 'CH3CN_K8'
        rest_values = [('CH3CN_K8',220.47581*u.GHz)]
    elif 'H30alpha' in cubefn:
        linefn = 'H30alpha'
        rest_values = [('H30alpha',231.90093*u.GHz)]



    for line,rest_value in rest_values:
        ch3cn = cube.with_spectral_unit(u.GHz).with_spectral_unit(u.km/u.s, rest_value=rest_value, velocity_convention='radio')
        ch3cn_slab = ch3cn.spectral_slab(54*u.km/u.s, 64*u.km/u.s)
        print(ch3cn_slab)

        vmap_ch3cn = ch3cn_slab.with_mask((ch3cn_slab > 5*u.mJy)|(ch3cn_slab < -5*u.mJy)).moment1()

        vmap_ch3cn.quicklook()
        vmap_ch3cn.FITSFigure.show_colorscale(cmap='viridis', vmin=54, vmax=64)

        vmap_ch3cn.FITSFigure.save(outfn.replace(".fits","_mom1.png").replace(linefn, line))

        vdispmap_ch3cn = ch3cn_slab.with_mask((ch3cn_slab > 5*u.mJy)|(ch3cn_slab < -5*u.mJy)).moment2()

        vdispmap_ch3cn.quicklook()
        vdispmap_ch3cn.FITSFigure.show_colorscale(cmap='viridis', vmin=0, vmax=5)

        vdispmap_ch3cn.FITSFigure.save(outfn.replace(".fits","_mom2.png").replace(linefn, line))
