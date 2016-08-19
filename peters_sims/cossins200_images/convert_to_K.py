from astropy.io import fits
from astropy import units as u
from astropy import constants
import radio_beam
from astropy.convolution import convolve_fft

def convert_to_K(fitsfilename):
    hdr = fits.getheader(fitsfilename)
    #pix_area = (hdr['CDELT1']*u.cm)**2
    factor = (constants.c**2 / (hdr['CRVAL3']*u.Hz)**2 / (2*constants.k_B)).to(u.K/(u.erg/u.s/u.cm**2/u.Hz))
    data = fits.getdata(fitsfilename)
    dataK = data*factor.value
    hdr['BUNIT'] = 'K'
    return fits.PrimaryHDU(data=dataK, header=hdr)

def convolve_to_beam(fitsfilename, beam=radio_beam.Beam(0.04*u.arcsec), distance=5400*u.pc):
    hdr = fits.getheader(fitsfilename)
    pix_area = (hdr['CDELT1']*u.cm)**2
    pix_area_arcsec = (pix_area/distance**2).to(u.arcsec**2, u.dimensionless_angles())
    kernel = beam.as_kernel(pix_area_arcsec**0.5)

    data = fits.getdata(fitsfilename)
    smoothed = convolve_fft(data, kernel)
    return fits.PrimaryHDU(data=smoothed, header=hdr)


if __name__ == "__main__":
    import glob
    for fn in glob.glob("*.fits"):
        if "_K.fits" in fn or "_K_smoothed.fits" in fn:
            continue
        fh = convert_to_K(fn)
        fh.writeto(fn.replace(".fits","_K.fits"), clobber=True)
        fh = convolve_to_beam(fn.replace(".fits","_K.fits"))
        fh.writeto(fn.replace(".fits","_K_smoothed.fits"))
