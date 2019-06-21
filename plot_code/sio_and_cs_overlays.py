import numpy as np
from astropy.io import fits
from astropy import wcs
#from astropy.visualization import AsinhStretch, imshow_norm, LogStretch, LinearStretch, MinMaxInterval, PercentileInterval
import astropy.visualization as vis
from astropy import convolution
from astropy import coordinates
from spectral_cube import SpectralCube
from astropy import units as u
import reproject
import paths
import pylab as pl

line = 'SiO'
e2siored = fits.open('/Users/adam/work/w51/alma/FITS//moments/w51_LB_SiO_red74to118_masked_cutoute2e.fits')
e2sioblue = fits.open('/Users/adam/work/w51/alma/FITS//moments/w51_LB_SiO_bluem32to55_masked_cutoute2e.fits')
line = 'CS21'
e2CSj21blue = fits.open('/Users/adam/work/w51/alma/FITS/longbaseline/{line}_m32to55kms_e2.fits'.format(line=line))
e2CSj21red = fits.open('/Users/adam/work/w51/alma/FITS/longbaseline/{line}_74to118kms_e2.fits'.format(line=line))
line = 'SiOJ21'
e2sioj21blue = fits.open('/Users/adam/work/w51/alma/FITS/longbaseline/{line}_m32to55kms_e2.fits'.format(line=line))
e2sioj21red = fits.open('/Users/adam/work/w51/alma/FITS/longbaseline/{line}_74to118kms_e2.fits'.format(line=line))

cont3mm = fits.open('/Users/adam/work/w51/alma/FITS/longbaseline/w51e2_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits.gz')

e2CSj21cube = SpectralCube.read('/Users/adam/work/w51/alma/FITS/longbaseline/velo_cutouts/w51e2e_csv0_j2-1_r0.5_medsub.fits').to(u.K)
cs21max = e2CSj21cube.spectral_slab(-32*u.km/u.s, 118*u.km/u.s).max(axis=0)

cs10cube = SpectralCube.read('/Users/adam/work/w51/vla_q/FITS/cutouts/W51e2e_CS_cutout.fits').to(u.K)
cs10max = cs10cube.spectral_slab(-32*u.km/u.s, 118*u.km/u.s).max(axis=0)
#cs10blue = cs10cube.spectral_slab(-32*u.km/u.s, 55*u.km/u.s).moment0(axis=0)
#cs10red = cs10cube.spectral_slab(74*u.km/u.s, 118*u.km/u.s).moment0(axis=0)
cs10mid = cs10cube.spectral_slab(55*u.km/u.s, 74*u.km/u.s).moment0(axis=0)

wcs_sio54 = wcs.WCS(e2siored[0].header)

cont3mm_proj,_ = reproject.reproject_interp((cont3mm[0].data.squeeze(),
                                             wcs.WCS(cont3mm[0].header).celestial),
                                            e2siored[0].header)


fig = pl.figure(1)
fig.set_size_inches(8,8)
fig.clf()
ax = pl.subplot(projection=wcs_sio54)

contnorm = vis.ManualInterval(-0.0005, 0.01)
siorednorm = vis.ManualInterval(-0.05, 0.15)
siobluenorm = vis.ManualInterval(-0.15, 0.60)

rgbim = np.array([siorednorm(e2siored[0].data),
                  contnorm(cont3mm_proj),
                  siobluenorm(e2sioblue[0].data)])

ax.imshow(rgbim.T.swapaxes(0,1), origin='lower', interpolation='none')

#csblue_reproj,_ = reproject.reproject_interp(e2sioj21blue, e2siored[0].header)
ax.contour(#e2CSj21blue[0].data,
           cs21max.value,
           #transform=ax.get_transform(wcs.WCS(e2CSj21blue[0].header)),
           transform=ax.get_transform(cs21max.wcs),
           levels=[2000, 4000, 6000], colors=[(1,1,1,0.9)]*10,
           linewidths=[0.5]*10)

# ADD CORRECTION from measured_vla_alma_offset
ww_cs10_corr = cs10max.wcs
ww_cs10_corr.wcs.radesys = 'ICRS'
#ww_cs10_corr.wcs.crval[0] -= -2.14775e-02/3600
#ww_cs10_corr.wcs.crval[1] += 1.98125e-02/3600

ax.contour(cs10max.value,
           transform=ax.get_transform(ww_cs10_corr),
           levels=[2000,4000,6000],
           colors=[(0,0,0,0.9)]*10,
           linewidths=[0.5]*10)

#csblue_smooth = convolution.convolve_fft(e2CSj21blue[0].data,
#                                         convolution.Gaussian2DKernel(2))
#ax.contour(csblue_smooth,
#           transform=ax.get_transform(wcs.WCS(e2CSj21blue[0].header)),
#           levels=[0.075], colors=[(1,1,1,0.9)]*10,
#           linewidths=[0.5]*10)
#
#
#csred_smooth = convolution.convolve_fft(e2CSj21red[0].data,
#                                        convolution.Gaussian2DKernel(2))
#ax.contour(csred_smooth,
#           transform=ax.get_transform(wcs.WCS(e2CSj21red[0].header)),
#           levels=[0.05, 0.1, 0.2], colors=[(1,0.5,0.2,0.9)]*4,
#           linewidths=[0.5]*10,
#          )

ax.contour(e2sioj21blue[0].data,
           transform=ax.get_transform(wcs.WCS(e2sioj21blue[0].header)),
           levels=[0.1, 0.2, 0.5, 0.8], colors=[(0.5,0.5,1,0.9)]*10,
           linewidths=[0.5]*10)
ax.contour(e2sioj21red[0].data,
           transform=ax.get_transform(wcs.WCS(e2sioj21red[0].header)),
           levels=[0.05, 0.1,], colors=[(1,0.5,0.5,0.9)]*10,
           linewidths=[0.5]*10)

def make_scalebar(ax, left_side, length, color='w', linestyle='-', label='',
                  fontsize=12, text_offset=0.1*u.arcsec, coordsys='icrs'):
    axlims = ax.axis()
    lines = ax.plot(u.Quantity([left_side.ra, left_side.ra-length]),
                    u.Quantity([left_side.dec]*2),
                    color=color, linestyle=linestyle, marker=None,
                    transform=ax.get_transform(coordsys),
                   )
    txt = ax.text((left_side.ra-length/2).to(u.deg).value,
                  (left_side.dec+text_offset).to(u.deg).value,
                  label,
                  verticalalignment='bottom',
                  horizontalalignment='center',
                  transform=ax.get_transform(coordsys),
                  color=color,
                  fontsize=fontsize,
                 )
    ax.axis(axlims)
    return lines,txt

make_scalebar(ax, coordinates.SkyCoord('19:23:43.95 +14:30:34.0', unit=(u.hour, u.deg), frame='icrs'),
              length=0.2*u.arcsec, label='0.2" = 1000 AU', fontsize=14, text_offset=0.02*u.arcsec)


ax.axis((234.0200833580487, 476.95387462618044, 268.73084003615594, 511.66463130428764))
ax.coords[0].set_ticklabel(size=16)
ax.coords[1].set_ticklabel(size=16)
ax.set_xlabel("RA (ICRS)", fontsize=18)
ax.set_ylabel("Dec (ICRS)", fontsize=18)
fig.savefig(paths.fpath('W51e2e_sio_outflow_with_CS_contours.pdf'), bbox_inches='tight')
fig.savefig(paths.fpath('W51e2e_sio_outflow_with_CS_contours.png'), bbox_inches='tight')
