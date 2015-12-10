import pyregion
from spectral_cube import SpectralCube

tmplt = "full_W51_spw{0}_lines.fits"

region_list = pyregion.open("cores.reg")

for spw in (0,1,2,3):
    cube = SpectralCube.read(tmplt.format(spw))
    print(cube)

    for reg in region_list:
        if 'text' not in reg.attr[1]:
            continue
        name = reg.attr[1]['text']
        if name:
            print("Extracting {0} from {1}".format(name, spw))
            SL = pyregion.ShapeList([reg])
            sc = cube.subcube_from_ds9region(SL)
            spec = sc.mean(axis=(1,2))
            spec.hdu.writeto("spectra/{0}_spw{1}_mean.fits".format(name, spw),
                             clobber=True)
