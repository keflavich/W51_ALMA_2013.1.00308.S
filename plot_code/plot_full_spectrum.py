"""
Create a one-page, printable version of a whole spectrum, ideally with line IDs
"""
import numpy as np
import pyspeckit
import line_to_image_list
from astropy import units as u
import paths
import pylab as pl
import matplotlib as mpl

line_ids = {line_to_image_list.labeldict[linename]:
            float(freqstr.strip('GHz'))*u.GHz
            for linename,freqstr,_,_ in line_to_image_list.line_to_image_list
            if linename in line_to_image_list.labeldict}

def plot_whole_spectrum(spectra, line_id=line_ids, velocity=55*u.km/u.s,
                        fignum=1, figsize=(8.3,11.7), figname=None,
                        title=None, do_line_ids=True,
                        fontsize=4, linewidth=0.25):
    pl.close(fignum)
    fig = pl.figure(fignum, figsize=figsize)
    
    with mpl.rc_context(rc={'font.size': 12,
                            'axes.labelsize':'medium',
                            'axes.titlesize':'medium'}):

        ax = fig.add_subplot(7,1,1)
        spectra[0].xarr.convert_to_unit(u.GHz)
        spectra[0].xarr.convert_to_unit(u.GHz)
        spectra[0].plotter(axis=ax)
        assert spectra[0].plotter.xlabel == 'Frequency (GHz)'
        if do_line_ids:
            spectra[0].plotter.line_ids(list(line_id.keys()),
                                        list(line_id.values()),
                                        velocity_offset=velocity,
                                        label1_size=fontsize,
                                        plot_kwargs={'linewidth':linewidth},
                                        max_iter=10,
                                       )

        for ii in (1,2,3):
            ax = fig.add_subplot(7,1,2*ii)
            spectra[ii].xarr.convert_to_unit(u.GHz)
            cropspec = spectra[ii][:spectra[ii].shape[0]/2]
            cropspec.plotter(axis=ax)
            if do_line_ids:
                cropspec.plotter.line_ids(list(line_id.keys()),
                                          list(line_id.values()),
                                          velocity_offset=velocity,
                                          label1_size=fontsize,
                                          plot_kwargs={'linewidth':linewidth},
                                          #annotate_kwargs={'fontsize':fontsize},
                                          max_iter=10,
                                         )

            ax = fig.add_subplot(7,1,2*ii+1)
            cropspec = spectra[ii][spectra[ii].shape[0]/2:]
            cropspec.plotter(axis=ax)
            if do_line_ids:
                cropspec.plotter.line_ids(list(line_id.keys()),
                                          list(line_id.values()),
                                          velocity_offset=velocity,
                                          label1_size=fontsize,
                                          plot_kwargs={'linewidth':linewidth},
                                          max_iter=10,
                                          #annotate_kwargs={'fontsize':fontsize},
                                         )

        for ii in range(1,8):
            ax = pl.subplot(7,1,ii)
            ax.set_ylabel("$T_B$ [K]")

        if title is not None:
            pl.subplot(7,1,1).set_title(title)

        if figname is not None:
            pl.savefig(paths.fpath(figname), dpi=300, bbox_inches='tight',
                       bbox_extra_artists=[])

if __name__ == "__main__":

    import glob
    import radio_beam
    from astropy.table import Table

    myvtbl = Table.read(paths.tpath('core_velocities.txt'),
                        # write as ascii.fixed_width
                        format='ascii.fixed_width', delimiter='|')

    line_table = Table.read(paths.apath('full_line_table.csv'))
    all_line_ids = {"{0}_{1}".format(row['Species'], row['Resolved QNs']):
                    (row['Freq-GHz'] if row['Freq-GHz']
                     else row['Meas Freq-GHz'])*u.GHz
                    for row in line_table}
    all_line_ids.update(line_ids)

    for row in myvtbl:

        speclist = [pyspeckit.Spectrum(fn) for fn in
                    glob.glob(paths.spath("{0}_spw*_peak.fits".format(row['source'])))]
        if len(speclist) == 0:
            continue

        for sp in speclist:
            beam = radio_beam.Beam.from_fits_header(sp.header)
            sp.data *= beam.jtok(sp.xarr)
            sp.unit='K'

        plot_whole_spectrum(speclist,
                            title=row['source'],
                            line_id=all_line_ids,
                            figname='fullspectra/ALLLINES_{0}.png'.format(row['source']),
                            velocity=row['velocity']*u.km/u.s,
                           )



    # default: show only line-to-image-list lines
    for row in myvtbl:

        speclist = [pyspeckit.Spectrum(fn) for fn in
                    glob.glob(paths.spath("{0}_spw*_peak.fits".format(row['source'])))]
        if len(speclist) == 0:
            continue

        for sp in speclist:
            beam = radio_beam.Beam.from_fits_header(sp.header)
            sp.data *= beam.jtok(sp.xarr)
            sp.unit='K'

        plot_whole_spectrum(speclist,
                            title=row['source'],
                            figname='fullspectra/{0}.png'.format(row['source']),
                            velocity=row['velocity']*u.km/u.s,
                           )


    # Methanol lines (for identification purposes - similar to ch3oh spectral
    # fit overlays...)
    from astroquery.splatalogue import Splatalogue
    Splatalogue.LINES_LIMIT=5000
    methanol_lines = Splatalogue.query_lines(218*u.GHz, 235*u.GHz,
                                            energy_max=3000,
                                            energy_type='eu_k',
                                            chemical_name='Methanol')

    from generic_lte_molecule_model import LTEModel
    methanolmodel = LTEModel(chemical_name='Methanol')
    # REQUIRES USING INTENSITY, not Aij
    methanolmodel_osu = LTEModel(chemical_name='Methanol', line_lists=['OSU'],
                                 freq_type='Meas Freq-GHz')

    for row in myvtbl:

        speclist = [pyspeckit.Spectrum(fn) for fn in
                    glob.glob(paths.spath("{0}_spw*_peak.fits".format(row['source'])))]
        if len(speclist) == 0:
            continue

        for sp in speclist:
            beam = radio_beam.Beam.from_fits_header(sp.header)
            sp.data *= beam.jtok(sp.xarr)
            sp.unit='K'

        methanol_line_ids = {"{0}_{1}".format(row['Species'], row['Resolved QNs']):
                            (row['Freq-GHz'] if row['Freq-GHz']
                             else row['Meas Freq-GHz'])*u.GHz
                            for row in methanol_lines}

        figname='fullspectra/methanol_{0}.png'.format(row['source'])
        plot_whole_spectrum(speclist,
                            title=row['source'],
                            line_id=methanol_line_ids,
                            figname=figname,
                            velocity=row['velocity']*u.km/u.s,
                           )

        spectra = pyspeckit.Spectra(speclist)
        spectra.xarr.convert_to_unit(u.GHz)
        mod = methanolmodel.lte_model(spectra.xarr,  row['velocity']*u.km/u.s,
                                      5*u.km/u.s, 500*u.K, 1e16*u.cm**-2)

        median = np.nanmedian(spectra.data)

        for ii in range(1,8):
            pl.subplot(7,1,ii)
            pl.plot(spectra.xarr, mod+median, color='r', linewidth=1, alpha=0.5)

        pl.savefig(paths.fpath(figname[:-4]+"_model.png"), dpi=300,
                   bbox_inches='tight', bbox_extra_artists=[])

    # Acetone lines (for identification purposes)
    from astroquery.splatalogue import Splatalogue
    Splatalogue.LINES_LIMIT=5000
    acetone_lines = Splatalogue.query_lines(218*u.GHz, 235*u.GHz,
                                            energy_max=800,
                                            energy_type='eu_k',
                                            chemical_name='Acetone')

    from generic_lte_molecule_model import LTEModel
    acetonemodel = LTEModel(chemical_name='Acetone')

    for row in myvtbl:

        speclist = [pyspeckit.Spectrum(fn) for fn in
                    glob.glob(paths.spath("{0}_spw*_peak.fits".format(row['source'])))]
        if len(speclist) == 0:
            continue

        for sp in speclist:
            beam = radio_beam.Beam.from_fits_header(sp.header)
            sp.data *= beam.jtok(sp.xarr)
            sp.unit='K'

        acetone_line_ids = {"{0}_{1}".format(row['Species'], row['Resolved QNs']):
                            (row['Freq-GHz'] if row['Freq-GHz']
                             else row['Meas Freq-GHz'])*u.GHz
                            for row in acetone_lines}

        figname='fullspectra/Acetone_{0}.png'.format(row['source'])
        plot_whole_spectrum(speclist,
                            title=row['source'],
                            line_id=acetone_line_ids,
                            figname=figname,
                            velocity=row['velocity']*u.km/u.s,
                           )

        spectra = pyspeckit.Spectra(speclist)
        spectra.xarr.convert_to_unit(u.GHz)
        mod = acetonemodel.lte_model(spectra.xarr,  row['velocity']*u.km/u.s, 5*u.km/u.s,
                                     500*u.K, 1e16*u.cm**-2)

        median = np.nanmedian(spectra.data)

        plot_whole_spectrum(speclist,
                            title=row['source'],
                            line_id=acetone_line_ids,
                            do_line_ids=False,
                            figname=None,
                            velocity=row['velocity']*u.km/u.s,
                           )

        for ii in range(1,8):
            pl.subplot(7,1,ii)
            pl.plot(spectra.xarr, mod+median, color='r', linewidth=1, alpha=0.5)

        pl.savefig(paths.fpath(figname[:-4]+"_model.png"), dpi=300,
                   bbox_inches='tight', bbox_extra_artists=[])


    # Formamide lines (for identification purposes)
    from astroquery.splatalogue import Splatalogue
    Splatalogue.LINES_LIMIT=5000
    formamide_lines = Splatalogue.query_lines(218*u.GHz, 235*u.GHz,
                                              energy_max=800,
                                              energy_type='eu_k',
                                              chemical_name='Formamide')

    from generic_lte_molecule_model import LTEModel
    formamidemodel = LTEModel(chemical_name='Formamide')

    for row in myvtbl:

        speclist = [pyspeckit.Spectrum(fn) for fn in
                    glob.glob(paths.spath("{0}_spw*_peak.fits".format(row['source'])))]
        if len(speclist) == 0:
            continue

        for sp in speclist:
            beam = radio_beam.Beam.from_fits_header(sp.header)
            sp.data *= beam.jtok(sp.xarr)
            sp.unit='K'

        formamide_line_ids = {"{0}_{1}".format(row['Species'], row['Resolved QNs']):
                              (row['Freq-GHz'] if row['Freq-GHz']
                               else row['Meas Freq-GHz'])*u.GHz
                              for row in formamide_lines}

        figname='fullspectra/Formamide_{0}.png'.format(row['source'])
        plot_whole_spectrum(speclist,
                            title=row['source'],
                            line_id=formamide_line_ids,
                            figname=figname,
                            velocity=row['velocity']*u.km/u.s,
                           )

        spectra = pyspeckit.Spectra(speclist)
        spectra.xarr.convert_to_unit(u.GHz)
        mod = formamidemodel.lte_model(spectra.xarr,  row['velocity']*u.km/u.s,
                                       5*u.km/u.s, 500*u.K, 1e16*u.cm**-2)

        median = np.nanmedian(spectra.data)

        plot_whole_spectrum(speclist,
                            title=row['source'],
                            line_id=formamide_line_ids,
                            do_line_ids=False,
                            figname=None,
                            velocity=row['velocity']*u.km/u.s,
                           )

        for ii in range(1,8):
            pl.subplot(7,1,ii)
            pl.plot(spectra.xarr, mod+median, color='r', linewidth=1, alpha=0.5)

        pl.savefig(paths.fpath(figname[:-4]+"_model.png"), dpi=300,
                   bbox_inches='tight', bbox_extra_artists=[])
