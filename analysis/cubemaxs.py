from spectral_cube import SpectralCube
import glob

for fn in glob.glob("full*fits"):
    cube = SpectralCube.read(fn)
    print(cube)
    cube.beam_threshold = 1e20
    max = cube.max(axis=0, how='slice')
    med = cube.median(axis=0, iterate_rays=True)
    max.write('moments/{0}_max.fits'.format(fn.replace(".fits","")))
    med.write('moments/{0}_med.fits'.format(fn.replace(".fits","")))
