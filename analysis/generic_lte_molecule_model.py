import numpy as np
from pyspeckit.spectrum.models import lte_molecule
from astropy import units as u
from astropy import constants
from astroquery.splatalogue import Splatalogue

try:
    from .ch3oh_rotational_diagram_maps import nupper_of_kkms
except (SystemError,ImportError):
    from ch3oh_rotational_diagram_maps import nupper_of_kkms


class LTEModel(object):
    def __init__(self, chemical_name, energy_max=2500,
                 nu_min=216*u.GHz, nu_max=235*u.GHz):
        slaim = Splatalogue.query_lines(nu_min, nu_max, chemical_name=chemical_name,
                                        energy_max=energy_max, energy_type='eu_k',
                                        line_lists=['SLAIM'],
                                        show_upper_degeneracy=True)
        self.freqs = np.array(slaim['Freq-GHz'])*u.GHz
        self.aij = slaim['Log<sub>10</sub> (A<sub>ij</sub>)']
        self.deg = slaim['Upper State Degeneracy']
        self.EU = (np.array(slaim['E_U (K)'])*u.K*constants.k_B).to(u.erg).value

    def lte_model(self, xarr, vcen, width, tex, column, background=None, tbg=2.73):

        if hasattr(tex,'unit'):
            tex = tex.value
        if hasattr(tbg,'unit'):
            tbg = tbg.value
        if hasattr(column, 'unit'):
            column = column.value
        if column < 25:
            column = 10**column
        if hasattr(vcen, 'unit'):
            vcen = vcen.value
        if hasattr(width, 'unit'):
            width = width.value

        ckms = constants.c.to(u.km/u.s).value

        freq = xarr.to(u.Hz).value # same unit as nu below
        model = np.zeros_like(xarr).value

        freqs_ = self.freqs.to(u.Hz).value

        #Q = specmodel.calculate_partitionfunction(result.data['States'],
        #                                          temperature=tex)[ch3ocho.Id]

        # use a very approximate Q_rot instead of a well-determined one
        Q = (self.deg * np.exp(-self.EU*u.erg / (constants.k_B * tex*u.K))).sum()

        for A, g, nu, eu in zip(self.aij, self.deg, freqs_, self.EU):
            taudnu = lte_molecule.line_tau_cgs(tex,
                                               column,
                                               Q,
                                               g,
                                               nu,
                                               eu,
                                               10**A)
            width_dnu = width / ckms * nu
            effective_linewidth_dnu = (2 * np.pi)**0.5 * width_dnu
            fcen = (1 - vcen/ckms) * nu
            tauspec = (np.exp(-(freq - fcen)**2 / (2 * width_dnu**2)) *
                       taudnu/effective_linewidth_dnu)
            jnu = (lte_molecule.Jnu_cgs(nu, tex)-lte_molecule.Jnu_cgs(nu, tbg))

            model = model + jnu*(1-np.exp(-tauspec))

        if background is not None:
            return background-model
        return model
