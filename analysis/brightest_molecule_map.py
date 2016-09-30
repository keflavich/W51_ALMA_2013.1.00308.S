import numpy as np
import paths
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import constants
from line_to_image_list import line_to_image_list
from radio_beam import Beam

template = 'merge/fullcube_cutouts/{sourcename}cutout_full_W51_7m12m_spw{0}_{res}lines.fits'

for res in ('','hires_'):
    for sourcename in ('e2','e8','north'):
        for spw in (0,1,2,3):
            cube = SpectralCube.read(paths.dpath(template.format(spw,
                                                                 sourcename=sourcename,
                                                                 res=res)))
            argmax = cube.argmax(axis=0)

            peakfreqmap = cube.spectral_axis[argmax]

            linefreqs = np.array([float(x[1].strip('GHz')) for x in line_to_image_list]) * (1-60*u.km/u.s/constants.c) * u.GHz
            linenames = np.array([x[0] for x in line_to_image_list])
            mapping = {x:y for x,y in zip(linenames, linefreqs)}
            closest = np.argmin(np.abs(peakfreqmap[None,:,:] - linefreqs[:,None,None]), axis=0)
            vdiff = np.min(np.abs(peakfreqmap[None,:,:] - linefreqs[:,None,None]), axis=0) / (220*u.GHz) * constants.c
            nomatch = vdiff > 10*u.km/u.s

            names = linenames[closest]
            freqs = u.Quantity([u.Quantity([mapping[y] for y in n]) for n in names])

            indices = np.unique(closest)
            unames = linenames[indices]
            uinds = np.arange(len(unames))

            closest_flat = np.empty_like(closest, dtype='float')
            for ii,jj in enumerate(indices):
                closest_flat[closest == jj] = ii
            closest_flat[nomatch] = np.nan

            import pylab as pl
            cm = pl.get_cmap('jet', len(unames))
            cm.set_bad('grey')
            pl.clf()
            pl.imshow(closest_flat, cmap=cm, vmin=-0.5, vmax=uinds.max()+0.5)
            cb = pl.colorbar()
            cb.set_ticks(uinds)
            cb.set_ticklabels(unames)

            pl.savefig(paths.fpath('peak_chem/{sourcename}_spw{spw}_{res}peakspecies.png'.format(sourcename=sourcename, spw=spw, res=res)))

            bad_beams = np.array([bm.major > 0.5*u.arcsec for bm in cube.beams], dtype='bool')
            cube = cube.with_mask((~bad_beams)[:,None,None])

            cube.beams = [Beam(0.35*u.arcsec, 0.35*u.arcsec, 0.0*u.deg)
                          if bm.major > 0.5*u.arcsec
                          else bm
                          for bm in cube.beams]

            jtok_equiv = cube._average_beams(1000).jtok_equiv(cube.with_spectral_unit(u.GHz).spectral_axis)
            kcube = cube.to(u.K, jtok_equiv)
            kcube.beam_threshold = 100000

            print("cube max")
            kmax = kcube.max(axis=0)
            kmax.write(paths.dpath('merge/moments/{sourcename}_spw{spw}_{res}max.fits'.format(sourcename=sourcename,
                                                                                              spw=spw,
                                                                                              res=res)),
                       overwrite=True
                      )

            print("cube median")
            kmed = kcube.median(axis=0)
            kmed.write(paths.dpath('merge/moments/{sourcename}_spw{spw}_{res}median.fits'.format(sourcename=sourcename,
                                                                                                 spw=spw,
                                                                                                 res=res)),
                       overwrite=True)

            print("cube 25th percentile")
            kpct25 = kcube.percentile(25, axis=0)
            kpct25.write(paths.dpath('merge/moments/{sourcename}_spw{spw}_{res}25thpct.fits'.format(sourcename=sourcename,
                                                                                                    spw=spw,
                                                                                                    res=res)),
                         overwrite=True)

            kmaxm25 = kmax - kpct25
            kmaxm25.write(paths.dpath('merge/moments/{sourcename}_spw{spw}_{res}max-25thpct.fits'.format(sourcename=sourcename,
                                                                                                         spw=spw,
                                                                                                         res=res)),
                          overwrite=True)
