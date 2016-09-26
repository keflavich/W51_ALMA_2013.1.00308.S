import numpy as np
from astropy import log
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
                 nu_min=216*u.GHz, nu_max=235*u.GHz, line_lists=['SLAIM'],
                 freq_type='Freq-GHz',
                 linelimit=100000,
                ):
        Splatalogue.LINES_LIMIT=linelimit
        line_ids = Splatalogue.get_species_ids().find(chemical_name)
        if len(line_ids) > 1:
            mwts = [int(lid[:3]) for lid in line_ids]
            if len(set(mwts)) != 1:
                raise ValueError("Chemical name {0} matches too many "
                                 "different lines: {1}".format(chemical_name,
                                                               line_ids))
            self.molwt = mwts[0]
        elif len(line_ids) == 1:
            self.molwt = int(list(line_ids.keys())[0][:3])
        else:
            raise ValueError("No matching molecules")
        linecat = Splatalogue.query_lines(nu_min, nu_max,
                                          chemical_name=chemical_name,
                                          energy_max=energy_max,
                                          energy_type='eu_k',
                                          line_strengths=['ls1','ls2','ls3','ls4','ls5'],
                                          line_lists=line_lists,
                                          show_upper_degeneracy=True)
        self.linecat = linecat
        self.freqs = np.array(linecat[freq_type])*u.GHz
        self.aij = linecat['Log<sub>10</sub> (A<sub>ij</sub>)']
        self.deg = [x if x > 0 else 1 for x in linecat['Upper State Degeneracy']]
        self.SijMu2 = linecat['S<sub>ij</sub>&#956;<sup>2</sup> (D<sup>2</sup>)']
        self.EU = (np.array(linecat['E_U (K)'])*u.K*constants.k_B).to(u.erg).value

        ## crappy debug hack
        #if 'OSU' in line_lists:
        #    self.EU = (np.array(linecat['E_L (K)'])*u.K*constants.k_B).to(u.erg).value

        slaim_query = Splatalogue.query_lines(1*u.Hz, 10000*u.GHz,
                                              chemical_name=chemical_name,
                                              energy_max=energy_max,
                                              energy_type='eu_k',
                                              line_lists=['SLAIM'],
                                              show_upper_degeneracy=True)
        self.slaim = slaim_query
        self.all_EU = (slaim_query['E_U (K)']*u.K*constants.k_B).to(u.erg).value
        self.all_freq = slaim_query['Freq-GHz']
        self.all_deg = slaim_query['Upper State Degeneracy']

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

        freq = xarr.to(u.Hz) # same unit as nu below
        model = np.zeros_like(xarr).value

        freqs_ = self.freqs.to(u.Hz)

        #Q = specmodel.calculate_partitionfunction(result.data['States'],
        #                                          temperature=tex)[ch3ocho.Id]

        # use a very approximate Q_rot instead of a well-determined one
        self.Q = Q = (self.all_deg * np.exp(-self.all_EU*u.erg /
                                            (constants.k_B * tex*u.K))).sum()

        for A, g, nu, eu, sijmu2 in zip(self.aij, self.deg, freqs_, self.EU, self.SijMu2):

            # skip lines that don't have an entry in the appropriate frequency
            # column (THIS IS BAD - it means we're not necessarily
            # self-consistently treating freqs above...)
            if nu == 0:
                continue

            width_dnu = width / ckms * nu
            effective_linewidth_dnu = (2 * np.pi)**0.5 * width_dnu
            fcen = (1 - vcen/ckms) * nu
            if np.isfinite(A) and A!=0 and g>0:
                taudnu = lte_molecule.line_tau(tex=tex*u.K,
                                               total_column=column*u.cm**-2,
                                               partition_function=Q, degeneracy=g,
                                               frequency=u.Quantity(nu, u.Hz),
                                               energy_upper=u.Quantity(eu,u.erg),
                                               einstein_A=(10**A)*u.s**-1,
                                              )
                #log.info("Line: {0} aij: {1}".format(nu, A))
                tauspec = (np.exp(-(freq - fcen)**2 / (2 * (width_dnu**2))) *
                           taudnu/effective_linewidth_dnu)
            else:
                #log.info("Line: {0} SijMu2: {1} tex: {2} degen: {3}".format(nu,
                #                                                            sijmu2,
                #                                                            tex,
                #                                                            g))
                tau = lte_molecule.line_tau_nonquantum(tex=tex*u.K,
                                                       total_column=column*u.cm**-2,
                                                       partition_function=Q,
                                                       degeneracy=g,
                                                       frequency=u.Quantity(nu, u.Hz),
                                                       energy_upper=u.Quantity(eu, u.erg),
                                                       SijMu2=sijmu2*u.debye**2,
                                                       molwt=self.molwt*u.Da,)
                tauspec = (np.exp(-(freq - fcen)**2 / (2 * (width_dnu**2))) *
                           tau)

            if np.any(np.isnan(tauspec)):
                raise ValueError("NaN encountered")

            jnu = (lte_molecule.Jnu_cgs(nu.to(u.Hz).value, tex) -
                   lte_molecule.Jnu_cgs(nu.to(u.Hz).value, tbg))

            model = model + jnu*(1-np.exp(-tauspec))

        if background is not None:
            return background-model
        return model
