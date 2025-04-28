from astropy.io import fits
import numpy as np
import os
import glob
import datetime
import sys
sys.path.append('.')
sys.path.append('/orange/adamginsburg/w51/W51_ALMA_2013.1.00308.S/script_longbaseline_b3')
from source_ids import sources_fmtd, source_field_mapping

if 'workdir' not in locals():
    workdir = os.getenv('WORK_DIR', '/red/adamginsburg/w51/')
os.chdir(workdir)

def makefits(myimagebase, cleanup=True):
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    # exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image

    if cleanup:
        for suffix in ('psf', 'weight', 'sumwt', 'pb', 'model', 'residual',
                       'mask', 'image', 'workdirectory'):
            os.system('rm -rf {0}.{1}'.format(myimagebase, suffix))

# These variables can be set from the command line when calling the script
# e.g., casa -c "sourcename='w51n'; spw=0; robust=0.5; suffix='clarkclean1e5'; niter=100000; execfile('scriptForImaging_line_cutouts_single.py')"
# or through environment variables: SOURCENAME, SPW, ROBUST, SUFFIX, NITER

# Check for environment variables first, then local variables
sourcename = os.getenv('SOURCENAME', None)
if sourcename is None and 'sourcename' in locals():
    sourcename = locals()['sourcename']
if sourcename is None:
    print("sourcename not defined, defaulting to 'w51n'")
    sourcename = 'w51n'

spw = os.getenv('SPW', None)
if spw is None and 'spw' in locals():
    spw = locals()['spw']
if spw is None:
    print("spw not defined, defaulting to 0")
    spw = 0
else:
    spw = int(spw)

robust = os.getenv('ROBUST', None)
if robust is None and 'robust' in locals():
    robust = locals()['robust']
if robust is None:
    print("robust not defined, defaulting to 0.5")
    robust = 0.5
else:
    robust = float(robust)

suffix = os.getenv('SUFFIX', None)
if suffix is None and 'suffix' in locals():
    suffix = locals()['suffix']
if suffix is None:
    print("suffix not defined, defaulting to 'clarkclean1e5'")
    suffix = 'clarkclean1e5'

niter = os.getenv('NITER', None)
if niter is None and 'niter' in locals():
    niter = locals()['niter']
if niter is None:
    print("niter not defined, defaulting to 100000")
    niter = int(1e5)
else:
    niter = int(niter)

print(f"Running with parameters: sourcename={sourcename}, spw={spw}, robust={robust}, suffix={suffix}, niter={niter}")

mslist = ['calibrated_final.ms']

# Define SPW mapping based on the SPW index
spw_mapping = [
    (0, 4, 8, 12, 16),   # index 0
    (1, 5, 9, 13, 17),   # index 1
    (2, 6, 10, 14, 18),  # index 2
    (3, 7, 11, 15, 19)   # index 3
]

spws = spw_mapping[spw]

# Create listobs if it doesn't exist
for ms in mslist:
    listobs_file = ms + '.listobs'
    if not os.path.exists(listobs_file):
        listobs(ms, listfile=listobs_file, overwrite=True)

imagename = f'W51{sourcename}_only.B3.robust{robust}.spw{spw}.{suffix}.1024'

# Check if file already exists and has the correct coordinate system
if os.path.exists("{0}.image.pbcor.fits".format(imagename)):
    if fits.getheader("{0}.image.pbcor.fits".format(imagename))['RADESYS'] == 'ICRS':
        print("Skipping completed file {0}".format(imagename))
        sys.exit(0)
    else:
        print("Redoing {0} because it's in fk5.".format(imagename))

if any([os.path.exists(f"{imagename}.{suffix}") for suffix in ('image', 'model', 'residual', 'psf', 'pb', 'sumwt', 'weight')]):
    print("Skipping {0} because imaging appears to be ongoing".format(imagename))
    sys.exit(0)

print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
tclean(vis=mslist,
       imagename=imagename,
       datacolumn='data',
       spw=",".join(['{0}'.format(ss) for ss in spws]),
       field=source_field_mapping[sourcename],
       specmode='cube',
       outframe='LSRK',
       threshold='2.5mJy',
       imsize=[1024, 1024],
       cell=['0.01arcsec'],
       niter=niter,
       cycleniter=-1, # -1 is default
       #cyclefactor=0.0001, # set very small: try to prevent major cycles
       phasecenter=sources_fmtd[sourcename],
       deconvolver='clark',
       gridder='standard',
       weighting='briggs',
       robust=robust,
       pbcor=True,
       pblimit=0.2,
       savemodel='none',
       #chanchunks=1,
       #parallel=True,
       interactive=False)
if os.path.exists(f"{imagename}.image"):
    makefits(imagename)
else:
    print(f"tclean failed for {imagename}")
    sys.exit(1)