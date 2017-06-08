import regions
from astropy import units as u
import paths
from gaussfit_catalog import gaussfit_catalog, gaussfit_image
from astropy.table import Table,Column

def tryint(x):
    try:
        return int(x)
    except:
        return x

def data_to_table(fit_data):
    names = fit_data.keys()
    numnames = [ii for ii,nm in enumerate(names)]
    stripnames = [nm for nm in names]
    stripnames = [fullname for nnm,fullname in sorted(zip(numnames,stripnames))]
    names = [fullname for nnm,fullname in sorted(zip(numnames,names))]
    namecol = Column(name='Name', data=stripnames)
    colnames = ['amplitude', 'center_x', 'center_y', 'fwhm_x', 'fwhm_y', 'pa',
                'chi2', 'chi2/n', 'e_amplitude', 'e_center_x', 'e_center_y',
                'e_fwhm_x', 'e_fwhm_y', 'e_pa', 'success',]
    columns = [Column(name=k, data=[fit_data[entry][k].value
                                    if hasattr(fit_data[entry][k],'value')
                                    else fit_data[entry][k]
                                    for entry in names],
                      unit=(fit_data[names[0]][k].unit
                            if hasattr(fit_data[names[0]][k], 'unit')
                            else None))
               for k in colnames]

    return Table([namecol]+columns)



if __name__ == "__main__":

    for regfn, contfn, name in (
        ('w51north_protostars.reg', 'W51n_cont_uniform.image.tt0.pbcor.fits', 'NorthUniform'),
        ('w51north_protostars.reg', 'W51n.cont.image.allEB.fits', 'NorthRobust'),
        ('w51e_protostars.reg', 'W51e2_cont_uniform.image.tt0.pbcor.fits', 'W51eUniform'),
        ('w51e_protostars.reg', 'W51e2_cont_briggs.image.fits', 'W51eRobust'),
    ):

        regs = regions.read_ds9(paths.rpath(regfn))

        contfn = paths.dpath('longbaseline/'+contfn)

        fit_data = gaussfit_catalog(contfn, regs, radius=0.1*u.arcsec,
                                    prefix=name+"_",
                                    max_radius_in_beams=5,
                                    max_offset_in_beams=2,
                                    savepath=paths.fpath('longbaseline/gaussfits'))

        tbl = data_to_table(fit_data)

        tbl.rename_column("chi2/n", "chi2_n")
        tbl.write(paths.tpath("gaussian_fit_table_{0}.ipac".format(name)),
                  format='ascii.ipac',
                  overwrite=True)
