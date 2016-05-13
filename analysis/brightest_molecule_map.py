import numpy as np
import paths
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import constants
from line_to_image_list import line_to_image_list

sourcename = 'e2'
for spw in (0,1,2,3):
    cube = SpectralCube.read(paths.dpath('merge/fullcube_cutouts/{sourcename}cutout_full_W51_7m12m_spw{0}_hires_lines.fits'.format(spw, sourcename=sourcename)))
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

    pl.savefig(paths.fpath('peak_chem/{sourcename}_spw{spw}_peakspecies.png'.format(sourcename=sourcename, spw=spw)))
