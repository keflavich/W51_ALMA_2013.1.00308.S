import numpy as np
import paths
import copy
from astropy.table import Table
from astropy import units as u
from astropy import coordinates
import pylab as pl
from astropy.io import fits
from astropy import wcs
from astropy.wcs import WCS as WCSaxes
from astropy.convolution import convolve, Gaussian2DKernel
import matplotlib

if int(matplotlib.__version__[0]) >= 2:
    pl.rcParams['figure.dpi'] = 75.
    pl.rcParams['savefig.dpi'] = 300.
    pl.rcParams['axes.labelsize'] = 16
    pl.rcParams['axes.titlesize'] = 16
    pl.rcParams['xtick.labelsize'] = 12
    pl.rcParams['ytick.labelsize'] = 12
    pl.rcParams['axes.linewidth'] = 0.15
    tick_fontsize = 16
    markersize = 0.5
else:
    tick_fontsize = 16
    markersize = 2

core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')
cores = coordinates.SkyCoord(core_phot_tbl['RA'], core_phot_tbl['Dec'],
                             frame='fk5')

w51moxc = Table.read(paths.tpath('w51_moxc.ipac'), format='ascii.ipac')
cmcontsrc = Table.read(paths.vpath('tables/EVLA_VLA_PointSourcePhotometry.ipac'),
                       format='ascii.ipac')
cmok = (cmcontsrc['Frequency'] == 5.9) & (cmcontsrc['Epoch'] == '3')
cmcoords = coordinates.SkyCoord(cmcontsrc['gracen'][cmok],
                                cmcontsrc['gdeccen'][cmok], frame='fk5')



hdu = fits.open(paths.vpath("data/W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits"))[0]

mywcs = wcs.WCS(hdu.header).sub([wcs.WCSSUB_CELESTIAL])
wcsaxes = WCSaxes(mywcs.to_header())

fig = pl.figure(1)
fig.clf()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)
im = ax.imshow(hdu.data.squeeze()*1e3, cmap=pl.cm.gray_r, origin='lower')
ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
ra.set_axislabel("RA (J2000)")
dec.set_axislabel("Dec (J2000)")
lims = ax.axis()
mylims = ((290.94465, 14.492588), (290.90729,14.530606))
(x1,y1),(x2,y2) = mywcs.wcs_world2pix(mylims, 0)

clims = mywcs.wcs_pix2world([[lims[0],lims[2]], [lims[1],lims[3]]], 0)
bins = [np.linspace(clims[1,0], clims[0,0], 200),
        np.linspace(clims[0,1], clims[1,1], 200)]

tr_fk5 = ax.get_transform("fk5")

coredots, = ax.plot(cores.ra, cores.dec, 'b.', transform=tr_fk5, markersize=2,
                    zorder=50, alpha=0.7)
cmdots, = ax.plot(cmcoords.ra, cmcoords.dec, 'g.', transform=tr_fk5, alpha=0.7,
                  zorder=1, markersize=2)

dots, = ax.plot(w51moxc['RAJ2000'], w51moxc['DEJ2000'], 'r.', transform=tr_fk5,
                markersize=1.5, zorder=60, alpha=0.7)
#ax.axis(lims)
ax.axis([x1,x2,y1,y2])

scalebar_right = coordinates.SkyCoord("19h23m40s", "14d29m45s", frame='fk5')
length = (0.5*u.pc / (5400*u.pc)).to(u.deg, u.dimensionless_angles())
ax.plot([scalebar_right.ra.deg, (scalebar_right.ra-length).deg]*u.deg,
        [(scalebar_right.dec).deg]*2*u.deg,
        'k',
        transform=tr_fk5,
        zorder=100, linewidth=2)
ax.text((scalebar_right.ra-length/2).deg, (scalebar_right.dec+0.001*u.deg).deg,
        "0.5 pc", color='k', transform=tr_fk5, ha='center')

fig.savefig(paths.fpath("moxc_points_on_cband_withcores.png"), bbox_inches='tight')
fig.savefig(paths.fpath("moxc_points_on_cband_withcores.pdf"), bbox_inches='tight')

H,bx,by = np.histogram2d(w51moxc['RAJ2000'], w51moxc['DEJ2000'], bins=bins)
H2 = convolve(H, Gaussian2DKernel(2))
cx = (bx[1:]+bx[:-1])/2.
cy = (by[1:]+by[:-1])/2.
con = ax.contour(cx,
                 cy,
                 H2.T, transform=tr_fk5,
                 levels=[0.12, 0.18, 0.24, 0.30, 0.36, 0.42],
                 colors=['r']*10,
                 linewidths=[0.5]*10,
                 alpha=0.5,
                 #interpolation='bicubic',
                )


H,bx,by = np.histogram2d(cores.ra.deg, cores.dec.deg, bins=bins)
H2 = convolve(H, Gaussian2DKernel(2))
cx = (bx[1:]+bx[:-1])/2.
cy = (by[1:]+by[:-1])/2.
con = ax.contour(cx,
                 cy,
                 H2.T, transform=tr_fk5,
                 levels=[0.12, 0.18, 0.24, 0.30, 0.36, 0.42],
                 colors=['b']*10,
                 linewidths=[0.5]*10,
                 alpha=0.5,
                 #interpolation='bicubic',
                )



#ax.axis(lims)
ax.axis([x1,x2,y1,y2])
fig.savefig(paths.fpath("moxc_points_contours_on_cband_withcores.png"), bbox_inches='tight')
fig.savefig(paths.fpath("moxc_points_contours_on_cband_withcores.pdf"), bbox_inches='tight')
dots.set_visible(False)
fig.savefig(paths.fpath("moxc_contours_on_cband_withcores.png"), bbox_inches='tight')
fig.savefig(paths.fpath("moxc_contours_on_cband_withcores.pdf"), bbox_inches='tight')
pl.draw()
pl.show()
