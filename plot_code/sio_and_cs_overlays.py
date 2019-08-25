import numpy as np
import radio_beam
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

#cont3mm = fits.open('/Users/adam/work/w51/alma/FITS/longbaseline/w51e2_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits.gz')
cont1mm = fits.open('/Users/adam/work/w51/alma/FITS/longbaseline/W51e2_cont_uniform.image.tt0.pbcor.fits')

e2CSj21cube = SpectralCube.read('/Users/adam/work/w51/alma/FITS/longbaseline/velo_cutouts/w51e2e_csv0_j2-1_r0.5_medsub.fits').to(u.K)
cs21max = e2CSj21cube.spectral_slab(-32*u.km/u.s, 118*u.km/u.s).max(axis=0)

cs10cube = SpectralCube.read('/Users/adam/work/w51/vla_q/FITS/cutouts/W51e2e_CS_cutout.fits').to(u.K)
cs10max = cs10cube.spectral_slab(-32*u.km/u.s, 118*u.km/u.s).max(axis=0)
#cs10blue = cs10cube.spectral_slab(-32*u.km/u.s, 55*u.km/u.s).moment0(axis=0)
#cs10red = cs10cube.spectral_slab(74*u.km/u.s, 118*u.km/u.s).moment0(axis=0)
cs10mid = cs10cube.spectral_slab(55*u.km/u.s, 74*u.km/u.s).moment0(axis=0)

wcs_sio54 = wcs.WCS(e2siored[0].header)

#cont3mm_proj,_ = reproject.reproject_interp((cont3mm[0].data.squeeze(),
#                                             wcs.WCS(cont3mm[0].header).celestial),
#                                            e2siored[0].header)
cont1mm_proj,_ = reproject.reproject_interp((cont1mm[0].data.squeeze(),
                                             wcs.WCS(cont1mm[0].header).celestial),
                                            e2siored[0].header)

fig = pl.figure(1)
fig.set_size_inches(8,8)
fig.clf()
ax = pl.subplot(projection=wcs_sio54)
pixscale = np.mean(wcs.utils.proj_plane_pixel_scales(wcs_sio54))*u.deg

contnorm = vis.ManualInterval(-0.0005, 0.007)
siorednorm = vis.ManualInterval(-0.05, 0.15)
siobluenorm = vis.ManualInterval(-0.15, 0.60)

rgbim = np.array([siorednorm(e2siored[0].data),
                  contnorm(cont1mm_proj),
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


beam_sio21 = radio_beam.Beam.from_fits_header(e2sioj21red[0].header)
beam_sio54 = radio_beam.Beam.from_fits_header(e2siored[0].header)
beam_cont = radio_beam.Beam.from_fits_header(cont1mm[0].header)

crd_21 = coordinates.SkyCoord('19:23:44.0025 +14:30:34.02', unit=(u.hour, u.deg), frame='icrs')
crd_54 = coordinates.SkyCoord('19:23:43.9975 +14:30:34.02', unit=(u.hour, u.deg), frame='icrs')
crd_cont = coordinates.SkyCoord('19:23:43.994 +14:30:34.02', unit=(u.hour, u.deg), frame='icrs')
pix_21 = wcs_sio54.wcs_world2pix(crd_21.ra.deg, crd_21.dec.deg, 0)
pix_54 = wcs_sio54.wcs_world2pix(crd_54.ra.deg, crd_54.dec.deg, 0)
pix_cont = wcs_sio54.wcs_world2pix(crd_cont.ra.deg, crd_cont.dec.deg, 0)
ell21 = beam_sio21.ellipse_to_plot(pix_21[0], pix_21[1], pixscale)
ell21.set_facecolor('none')
ell21.set_edgecolor('w')
ell54 = beam_sio54.ellipse_to_plot(pix_54[0], pix_54[1], pixscale)
ell54.set_facecolor('blue')
ell54.set_edgecolor('blue')
ellcont = beam_cont.ellipse_to_plot(pix_cont[0], pix_cont[1], pixscale)
ellcont.set_facecolor('green')
ellcont.set_edgecolor('green')

ax.add_artist(ell21)
ax.add_artist(ell54)
ax.add_artist(ellcont)


ax.axis((234.0200833580487, 476.95387462618044, 268.73084003615594, 511.66463130428764))
ax.coords[0].set_ticklabel(size=16)
ax.coords[1].set_ticklabel(size=16)
ax.set_xlabel("RA (ICRS)", fontsize=18)
ax.set_ylabel("Dec (ICRS)", fontsize=18)




## Add inset axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import astropy.visualization


xslice = slice(365, 385)
yslice = slice(340, 365)
xslice = slice(None)
yslice = slice(None)

axins = zoomed_inset_axes(ax, zoom=4, loc=2,
                          bbox_to_anchor=[0.15,0.85],
                          bbox_transform=fig.transFigure,
                          axes_class=astropy.visualization.wcsaxes.core.WCSAxes,
                          axes_kwargs=dict(wcs=wcs_sio54[xslice, yslice]))

rgbim[0,:,:] = 0
rgbim[2,:,:] = 0
axins.imshow(rgbim.T.swapaxes(0,1)[xslice, yslice, :],
             origin='lower', interpolation='none',
            )

mark_inset(parent_axes=ax, inset_axes=axins,
           loc1=1, loc2=3, fc="none", ec="0.5",
           lw=0.5)
axins.set_xticklabels([])
axins.set_yticklabels([])

lon = axins.coords['ra']
lat = axins.coords['dec']
lon.set_ticklabel_visible(False)
lat.set_ticklabel_visible(False)
lon.set_ticks_visible(False)
lat.set_ticks_visible(False)

axins.contour(cs10max.value,
              transform=axins.get_transform(ww_cs10_corr),
              levels=[2000,4000,6000],
              colors=[(0,0,0,0.9)]*10,
              linewidths=[0.5]*10)
axins.contour(cs21max.value,
              transform=axins.get_transform(cs21max.wcs),
              levels=[2000, 4000, 6000],
              colors=[(1,1,1,0.9)]*10,
              linewidths=[0.5]*10)

axins.axis((338.0, 360.0, 364.0, 386.0,))
#axins.axis((0,19,0,19))

ellcont = beam_cont.ellipse_to_plot(342, 367, pixscale)
ellcont.set_facecolor('green')
ellcont.set_edgecolor('green')

axins.add_artist(ellcont)


fig.savefig(paths.fpath('W51e2e_sio_outflow_with_CS_contours.pdf'), bbox_inches='tight')
fig.savefig(paths.fpath('W51e2e_sio_outflow_with_CS_contours.png'), bbox_inches='tight')
