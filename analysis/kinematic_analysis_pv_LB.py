""" LB = long baseline """
import numpy as np

import pvextractor
import os
import glob
import paths
from astropy import units as u
from astropy import coordinates
from astropy import wcs
from spectral_cube import SpectralCube
import pyregion
from line_to_image_list import line_to_image_list

# use outflow_meta b/c higher precision than ds9 reg
from outflow_meta import e2e, e8, north, lacy
from line_point_offset import offset_to_point

import pylab as pl

# not used, but could be to mark location of central "source" with a vertical line
# e2ereg = pyregion.open(paths.rpath('w51e2e.reg'))[0]
# w51e2e = coordinates.SkyCoord(e2ereg.coord_list[0]*u.deg,
#                               e2ereg.coord_list[1]*u.deg, frame='fk5')

for ii,direction in enumerate(('perpco', 'perpsio')):
    diskycoorddict = {}
    for source in ('e2e','e8','north','lacy'):
        diskycoord_list = pyregion.open(paths.rpath("{0}_disk_pvextract.reg"
                                                    .format(source)))[ii].coord_list
        diskycoords = coordinates.SkyCoord(["{0} {1}".format(diskycoord_list[jj],
                                                             diskycoord_list[jj+1])
                                            for jj in range(0,
                                                            len(diskycoord_list),
                                                            2)], unit=(u.deg,
                                                                       u.deg),
                                           frame='fk5')
        diskycoorddict[source] = diskycoords

    diskycoorddict['e2'] = diskycoorddict['e2e']


    for name, cutoutname, source, vrange in (
        ('lacy', 'northwest', lacy, (50,75)),
        ('north', 'north', north, (45,75)),
        ('e8', 'e8', e8, (45,75)),
        ('e2', 'e2', e2e, (45,70)),
       ):

        diskycoords = diskycoorddict[name]

        for fn in glob.glob(paths.lbpath("W51{0}cax*ALL_medsub_cutout.fits".format(cutoutname))):

            print(fn)
            try:
                cube = SpectralCube.read(fn).with_spectral_unit(u.GHz)
            except TypeError as ex:
                print(ex)
                continue

            extraction_path = pvextractor.Path(diskycoords, 0.05*u.arcsec)

            for line, restfreq, velocity_res, spw in line_to_image_list:

                basename = line

                frq = float(restfreq.strip('GHz')) * u.GHz
                vcube = cube.with_spectral_unit(u.km/u.s,
                                                velocity_convention='radio',
                                                rest_value=frq)
                svcube = vcube.spectral_slab((vrange[0]-1)*u.km/u.s,
                                             (vrange[1]+1)*u.km/u.s)
                if svcube.shape[0] <= 5:
                    print("SKIPPING {3} {0} {2}: {1}".format(line, restfreq, direction, name))
                    continue

                print("Extracting {3} {0} {2}: {1}".format(line, restfreq, direction, name))

                extracted = pvextractor.extract_pv_slice(svcube, extraction_path)
                #extracted.writeto(outfn, clobber=True)

                ww = wcs.WCS(extracted.header)
                ww.wcs.cdelt[1] /= 1000.0
                ww.wcs.crval[1] /= 1000.0
                ww.wcs.cunit[1] = u.km/u.s
                ww.wcs.cdelt[0] *= 3600
                ww.wcs.cunit[0] = u.arcsec

                fig = pl.figure(1)
                fig.clf()
                ax = fig.add_axes([0.15, 0.1, 0.8, 0.8],projection=ww)
                ax.imshow(extracted.data, cmap='viridis', vmin=-0.005)
                ax.set_xlabel("Offset [\"]")
                ax.set_ylabel("$V_{LSR}$ [km/s]")
                #ax.set_aspect(25)
                ax.set_aspect(extracted.data.shape[1]/extracted.data.shape[0])

                good_limits = (np.argmax(np.isfinite(extracted.data.max(axis=0))),
                               extracted.data.shape[1] -
                               np.argmax(np.isfinite(extracted.data.max(axis=0)[::-1])))
                leftmost_position = ww.wcs_pix2world(good_limits[0],
                                                     vrange[0],
                                                     0)[0]*u.arcsec

                trans = ax.get_transform('world')
                length = (1000*u.au / (5400*u.pc)).to(u.deg, u.dimensionless_angles())
                endpoints_x = u.Quantity([0.5*u.arcsec, 0.5*u.arcsec+length]) + leftmost_position
                ax.plot(endpoints_x.to(u.arcsec),
                        [vrange[0]+5]*2*u.km/u.s,
                        'w',
                        transform=trans,
                        zorder=100, linewidth=2)
                ax.text(endpoints_x.mean().value,
                        (vrange[0]+6),
                        "1000 au", color='w', transform=trans, ha='center')

                origin = offset_to_point(source.ra.deg,
                                         source.dec.deg,
                                         extraction_path)*u.deg

                ax.vlines(origin.to(u.arcsec).value, vrange[0]-5, vrange[1]+5,
                          color='w', linestyle='--', linewidth=2.0,
                          alpha=0.6, transform=trans)

                ax.set_ylim(ww.wcs_world2pix(0,vrange[0],0)[1],
                            ww.wcs_world2pix(0,vrange[1],0)[1])
                ax.set_xlim(good_limits)

                fig.savefig(paths.fpath('lbpv/{0}/{0}_{1}_{2}_lbpv.png'
                                        .format(name, basename, direction)),
                            bbox_inches='tight')


            #outflow_coords = coordinates.SkyCoord(["19:23:44.127 +14:30:32.30", "19:23:43.822 +14:30:36.64"], unit=(u.hour, u.deg), frame='fk5')
            #outflowpath = pvextractor.Path(outflow_coords, 0.2*u.arcsec)
            #extracted = pvextractor.extract_pv_slice(medsub, outflowpath)
            #extracted.writeto('W51e2_PV_outflowaxis_spw{0}.fits'.format(ii), clobber=True)

    # outflow_coords = coordinates.SkyCoord(["19:23:44.127 +14:30:32.30", "19:23:43.822 +14:30:36.64"], unit=(u.hour, u.deg), frame='fk5')
    # outflowpath = pvextractor.Path(outflow_coords, 0.2*u.arcsec)
    # 
    # if __name__ == '__main__':
    #     files = glob.glob("/Volumes/passport/w51-alma/*.fits")
    # 
    #     for fn in files:
    #         try:
    #             cube = SpectralCube.read(fn)
    #         except:
    #             continue
    #         cube.allow_huge_operations=True
    #         print(cube)
    # 
    #         out = os.path.splitext(os.path.split(fn)[1])[0]
    #         print(out)
    # 
    #         diskaxis_out = 'W51e2_PV_diskaxis_{0}.fits'.format(out)
    #         if not os.path.exists(diskaxis_out):
    #             extracted = pvextractor.extract_pv_slice(cube, diskypath)
    #             extracted.writeto(diskaxis_out, clobber=True)
    # 
    #         outflowaxis_out = 'W51e2_PV_outflowaxis_{0}.fits'.format(out)
    #         if not os.path.exists(outflowaxis_out):
    #             extracted = pvextractor.extract_pv_slice(cube, outflowpath)
    #             extracted.writeto(outflowaxis_out, clobber=True)
