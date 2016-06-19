"""
Create some simple grids for the 1mm para-H2CO lines

Used in the W51 ALMA H2CO project


Grid shape is [Temperature, Density, Column Density]
"""
from __future__ import print_function
import pyradex
import pyradex.fjdu
import numpy as np
from astropy.utils.console import ProgressBar
from astropy.io import fits
from astropy import log
import warnings
import os
# Make sure warnings are only shown once so the progressbar doesn't get flooded
warnings.filterwarnings('once')

ntemp,ndens,ncol = 50,20,30

temperatures = np.linspace(10,250,ntemp)
densities = np.linspace(2.5,7,ndens)
columns = np.linspace(11, 16.1, ncol)
abundance = 3e-9 # approx 10^-8.5
opr = 0.01 # assume primarily para
opr = 3
fortho = opr/(1+opr)

if not os.path.exists('ph2co-h2.dat'):
    import urllib.request
    urllib.request.urlretrieve('http://home.strw.leidenuniv.nl/~moldata/datafiles/ph2co-h2.dat')

def compute_grid(densities=densities, temperatures=temperatures,
                 columns=columns, fortho=fortho, deltav=1.0,
                 escapeProbGeom='lvg', Radex=pyradex.Radex,
                 run_kwargs={'reuse_last': True, 'reload_molfile': False}):

    # Initialize the RADEX fitter with some reasonable parameters
    R = Radex(species='ph2co-h2',
              column=1e14,
              temperature=50,
              escapeProbGeom=escapeProbGeom,
              collider_densities={'oH2':2e4*fortho,'pH2':2e4*(1-fortho)})

    R.run_radex()
    R.maxiter = 200

    # get the table so we can look at the frequency grid
    #table = R.get_table()

    # Get warnings about temperature early
    R.temperature = temperatures.min()
    R.temperature = temperatures.max()

    # Target frequencies:
    #table[np.array([6,1,11])].pprint()

    key_303 = 2
    key_321 = 9
    key_322 = 12
    key_404 = 3
    key_423 = 16
    key_422 = 19
    #key_303 = np.where((table['upperlevel'] == '3_0_3') &
    #                   (table['frequency'] > 218) &
    #                   (table['frequency'] < 220))[0]
    #key_321 = np.where((table['upperlevel'] == '3_2_1') &
    #                   (table['frequency'] > 218) &
    #                   (table['frequency'] < 220))[0]
    #key_322 = np.where((table['upperlevel'] == '3_2_2') &
    #                   (table['frequency'] > 218) &
    #                   (table['frequency'] < 220))[0]

    # used to assess where the grid failed
    bad_pars = []

    ntemp = len(temperatures)
    ndens = len(densities)
    ncols = len(columns)

    shape = [ntemp,ndens,ncols,]

    pars = dict(
        taugrid_303=np.full(shape, np.nan),
        texgrid_303=np.full(shape, np.nan),
        fluxgrid_303=np.full(shape, np.nan),
        taugrid_321=np.full(shape, np.nan),
        texgrid_321=np.full(shape, np.nan),
        fluxgrid_321=np.full(shape, np.nan),
        taugrid_322=np.full(shape, np.nan),
        texgrid_322=np.full(shape, np.nan),
        fluxgrid_322=np.full(shape, np.nan),
        taugrid_404=np.full(shape, np.nan),
        texgrid_404=np.full(shape, np.nan),
        fluxgrid_404=np.full(shape, np.nan),
        taugrid_423=np.full(shape, np.nan),
        texgrid_423=np.full(shape, np.nan),
        fluxgrid_423=np.full(shape, np.nan),
        taugrid_422=np.full(shape, np.nan),
        texgrid_422=np.full(shape, np.nan),
        fluxgrid_422=np.full(shape, np.nan),
    )

    for iTem,tt in enumerate(ProgressBar(temperatures)):
        R.temperature = tt
        for iDens,dd in enumerate(densities):
            R.density = {'oH2':10**dd*fortho,'pH2':10**dd*(1-fortho)}
            for iCol,cc in enumerate(columns):
                #R.abundance = abundance # reset column to the appropriate value
                R.column_per_bin = 10**cc
                R.deltav = deltav
                #niter = R.run_radex(reuse_last=False, reload_molfile=True)
                niter = R.run_radex(**run_kwargs)

                if niter == R.maxiter:
                    bad_pars.append([tt,dd,cc])

                TI = R.source_line_surfbrightness
                pars['taugrid_303'][iTem,iDens,iCol] = R.tau[key_303]
                pars['texgrid_303'][iTem,iDens,iCol] = R.tex[key_303].value
                pars['fluxgrid_303'][iTem,iDens,iCol] = TI[key_303].value
                pars['taugrid_321'][iTem,iDens,iCol] = R.tau[key_321]
                pars['texgrid_321'][iTem,iDens,iCol] = R.tex[key_321].value
                pars['fluxgrid_321'][iTem,iDens,iCol] = TI[key_321].value
                pars['taugrid_322'][iTem,iDens,iCol] = R.tau[key_322]
                pars['texgrid_322'][iTem,iDens,iCol] = R.tex[key_322].value
                pars['fluxgrid_322'][iTem,iDens,iCol] = TI[key_322].value
                pars['taugrid_404'][iTem,iDens,iCol] = R.tau[key_404]
                pars['texgrid_404'][iTem,iDens,iCol] = R.tex[key_404].value
                pars['fluxgrid_404'][iTem,iDens,iCol] = TI[key_404].value
                pars['taugrid_422'][iTem,iDens,iCol] = R.tau[key_422]
                pars['texgrid_422'][iTem,iDens,iCol] = R.tex[key_422].value
                pars['fluxgrid_422'][iTem,iDens,iCol] = TI[key_422].value
                pars['taugrid_423'][iTem,iDens,iCol] = R.tau[key_423]
                pars['texgrid_423'][iTem,iDens,iCol] = R.tex[key_423].value
                pars['fluxgrid_423'][iTem,iDens,iCol] = TI[key_423].value

    return (TI, pars, bad_pars)

def makefits(data, btype, densities=densities, temperatures=temperatures,
             columns=columns, ):

    newfile = fits.PrimaryHDU(data=data)
    newfile.header.update('BTYPE',  btype)


    newfile.header.update('CRVAL1',  min(columns))
    newfile.header.update('CRPIX1',  1)
    newfile.header.update('CDELT1', columns[1]-columns[0])
    newfile.header.update('CTYPE1',  'LOG-COLU')

    newfile.header.update('CRVAL2',  min(densities))
    newfile.header.update('CRPIX2',  1)
    newfile.header.update('CDELT2', densities[1]-densities[0])
    newfile.header.update('CTYPE2',  'LOG-DENS')

    newfile.header.update('CRVAL3',  (min(temperatures)))
    newfile.header.update('CRPIX3',  1)
    if len(np.unique(temperatures)) == 1:
        newfile.header.update('CTYPE3',  'ONE-TEMP')
        newfile.header.update('CDELT3', temperatures[0])
    else:
        newfile.header.update('CTYPE3',  'LIN-TEMP')
        newfile.header.update('CDELT3', (np.unique(temperatures)[1]) - (np.unique(temperatures)[0]))
    return newfile

if __name__ == "__main__":
    import re
    from paths import gpath
    bt = re.compile("tex|tau|flux")

    (fTI, fpars, fbad_pars) = compute_grid(Radex=pyradex.fjdu.Fjdu,
                                           run_kwargs={})
    
    for pn in fpars:
        btype = bt.search(pn).group()
        ff = makefits(fpars[pn], btype, densities=densities,
                      temperatures=temperatures, columns=columns)
        outfile = 'fjdu_pH2CO_{line}_{type}_{dv}.fits'.format(line=pn[-3:],
                                                              type=btype,
                                                              dv='1kms')
        ff.writeto(gpath(outfile),
                   clobber=True)
        print(outfile)

    ff = makefits(fpars['fluxgrid_321']/fpars['fluxgrid_303'], 'ratio',
                  densities=densities, temperatures=temperatures,
                  columns=columns)
    outfile = 'fjdu_pH2CO_{line}_{type}_{dv}.fits'.format(line='321to303',
                                                          type='ratio',
                                                          dv='1kms')
    ff.writeto(gpath(outfile), clobber=True)

    (TI, pars, bad_pars) = compute_grid()

    for pn in pars:
        btype = bt.search(pn).group()
        ff = makefits(pars[pn], btype, densities=densities,
                      temperatures=temperatures, columns=columns)
        outfile = 'pH2CO_{line}_{type}_{dv}.fits'.format(line=pn[-3:],
                                                         type=btype, dv='1kms')
        ff.writeto(gpath(outfile),
                   clobber=True)
        print(outfile)

    ff = makefits(pars['fluxgrid_321']/pars['fluxgrid_303'], 'ratio',
                  densities=densities, temperatures=temperatures,
                  columns=columns)
    outfile = 'pH2CO_{line}_{type}_{dv}.fits'.format(line='321to303',
                                                     type='ratio', dv='1kms')
    ff.writeto(gpath(outfile), clobber=True)

    log.info("FJDU had {0} bad pars".format(len(fbad_pars)))
    log.info("RADEX had {0} bad pars".format(len(bad_pars)))
    

    # look at differences
    for pn in pars:
        btype = bt.search(pn).group()
        outfile = 'pH2CO_{line}_{type}_{dv}.fits'.format(line=pn[-3:],
                                                         type=btype, dv='1kms')
        header = fits.getheader(gpath(outfile))
        im1 = fits.getdata(gpath('fjdu_'+outfile))
        im2 = fits.getdata(gpath(outfile))
        hdu = fits.PrimaryHDU(data=im1-im2, header=header)
        hdu.writeto(gpath('diff_fjdu-radex_'+outfile), clobber=True)
