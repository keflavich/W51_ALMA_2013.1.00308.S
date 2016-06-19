import inspect
import time
import warnings

import numpy as np
from scipy.ndimage.interpolation import map_coordinates
from astropy import units as u
from astropy import log
from scipy import stats

from h2co_modeling import grid_fitter

class generic_paraH2COmodel(object):
    def grid_getmatch_321to303(self, ratio, eratio, chi2_thresh=1):
            match,indbest,chi2r = grid_fitter.grid_getmatch(ratio, eratio,
                                                            self.modelratio1,
                                                            chi2_thresh=chi2_thresh)
            return chi2r

    def grid_getmatch_322to321(self, ratio, eratio, chi2_thresh=1):
            match,indbest,chi2r = grid_fitter.grid_getmatch(ratio, eratio,
                                                            self.modelratio2,
                                                            chi2_thresh=chi2_thresh)
            return chi2r

    def grid_getmatch_303(self, taline303, etaline303, chi2_thresh=1):
            match,indbest,chi2r = grid_fitter.grid_getmatch(taline303, etaline303,
                                                            self.tline303,
                                                            chi2_thresh=chi2_thresh)
            return chi2r

    def grid_getmatch_321(self, taline321, etaline321, chi2_thresh=1):
            match,indbest,chi2r = grid_fitter.grid_getmatch(taline321, etaline321,
                                                            self.tline321,
                                                            chi2_thresh=chi2_thresh)
            return chi2r

    def grid_getmatch_322(self, taline322, etaline322, chi2_thresh=1):
            match,indbest,chi2r = grid_fitter.grid_getmatch(taline322, etaline322,
                                                            self.tline322,
                                                            chi2_thresh=chi2_thresh)
            return chi2r

    @property
    def chi2(self):
        return self._chi2

    @chi2.setter
    def chi2(self, value):
        self._chi2 = value
        self._likelihood = np.exp(-value/2)
        self._likelihood /= self._likelihood.sum()

    @property
    def likelihood(self):
        return self._likelihood

    def chi2_fillingfactor(self, tline, etline, lineid):
        """
        Return a chi^2 value for each model parameter treating the specified
        line brightness as a lower limit

        Parameters
        ----------
        tline : float
            The line brightness temperature
        lineid : int
            The line id, one of 303,321,322
        """
        chi2 = ((self.tline[lineid] - tline)/etline)**2 * (self.tline[lineid] < tline)
        return chi2

    def chi2_column(self, logh2column, elogh2column, h2coabundance, linewidth):
        """
        Linewidth discussion:
            the grid is in cm^-2 / km/s /pc
            the grid linewidth is a gradient, km/s/pc
            X = N(H2CO) * gradient / density
            if the line is wider than the gradient, the column
            is spread among multiple 'cells'
            ncells = linewidth / gradient
            N(H2) = n(H2) * ncells * cellsize
            N(H2) = N_real_tot / X
                  = N_real / X * ncells
                  = N_lvg * G * L / X * ncells
                  = N_lvg * linewidth * L / X
                  L = 1 pc
            What is N_lvg? G=gradient L=length
            N_lvg = N_real / G / L
            X = N_real / (n(H2) * L)
              = N_lvg * G / (n(H2) * L->cm)

        ncells = linewidth(FWHM) / grid_gradient (NOT the integral)
        """

        h2fromh2co = (self.columnarr
                      + np.log10(linewidth)
                      - h2coabundance)
        chi2_h2 = ((h2fromh2co-logh2column)/elogh2column)**2

        return chi2_h2

    def chi2_abundance(self, logabundance, elogabundance):
        """
        See chi2_column for linewidth discussion
        """
        model_logabundance = (self.columnarr
                              + np.log10(self.grid_linewidth)
                              - np.log10(u.pc.to(u.cm))
                              - self.densityarr)
        chi2X = ((model_logabundance-logabundance)/elogabundance)**2
        return chi2X


    def get_parconstraints(self,
                           nsigma=1):
        """
        If parameter constraints have been set with set_constraints or
        set_constraints_fromrow

        Parameters
        ----------
        nsigma : float
            The number of sigmas to go out to when determining errors
        """
        if not hasattr(self, 'chi2'):
            raise AttributeError("Run set_constraints first")

        row = {}

        inds = np.argsort(self.likelihood.flat)
        cdf = np.cumsum(self.likelihood.flat[inds])
        frac_above = (stats.norm.cdf(nsigma)-stats.norm.cdf(-nsigma))
        cdfmin = np.argmin(np.abs(cdf - (1-frac_above)))
        sigma_like = self.likelihood.flat[inds][cdfmin]

        indbest = np.argmax(self.likelihood)
        # Compute the *marginal* 1-sigma regions
        for parname,pararr,ax in zip(('temperature','column','density'),
                                     (self.temparr,self.columnarr,self.densityarr),
                                     (0,2,1)):
            row['{0}_chi2'.format(parname)] = pararr.flat[indbest]
            row['expected_{0}'.format(parname)] = ((pararr*self.likelihood).sum() / self.likelihood.sum())

            axes = tuple(x for x in (0,1,2) if x != ax) 
            like = self.likelihood.sum(axis=axes)
            cdf_inds = np.argsort(like)
            ppf = 1-like[cdf_inds].cumsum()
            cutoff_like = like[cdf_inds[np.argmin(np.abs(ppf-frac_above))]]
            selection = like > cutoff_like

            slc = [slice(None) if x==ax else 0 for x in (0,1,2)]
            pararr = pararr[slc]

            if np.abs(like[selection].sum() - frac_above) > 0.05:
                # we want the sum of the likelihood to be right!
                #import ipdb; ipdb.set_trace()
                warnings.warn("Likelihood is not self-consistent.")

            if np.count_nonzero(selection) > 0:
                row['{0:1.1s}min1sig_chi2'.format(parname)] = pararr[selection].min()
                row['{0:1.1s}max1sig_chi2'.format(parname)] = pararr[selection].max()
            else:
                row['{0:1.1s}min1sig_chi2'.format(parname)] = np.nan
                row['{0:1.1s}max1sig_chi2'.format(parname)] = np.nan

        for parname in ('logh2column', 'elogh2column', 'logabundance',
                        'elogabundance'):
            if hasattr(self, parname):
                row[parname] = getattr(self, parname)

        self._parconstraints = row

        return row

    @property
    def parconstraints(self):
        if not hasattr(self,'_parconstraints') or self._parconstraints is None:
            return self.get_parconstraints()
        else:
            return self._parconstraints


    def denstemplot(self):
        self.parplot('dens','tem')

    def denscolplot(self):
        self.parplot('col','dens')

    def coltemplot(self):
        self.parplot('col','tem')
