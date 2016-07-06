from h2co_modeling import lte_model
from constrain_parameters import paraH2COmodel
import paths
import numpy as np
import pylab as pl

modelratio321303 = lte_model.T_321/lte_model.T_303
modelratio322303 = lte_model.T_322/lte_model.T_303

def ratio_to_temperature(data, modelratio=modelratio321303):
    vals = data[np.isfinite(data)]
    tems = np.interp(vals, modelratio[np.isfinite(modelratio)],
                     np.array(lte_model.tem)[np.isfinite(modelratio)])
    newr = data.copy()
    newr[np.isfinite(data)] = tems
    return newr

data = np.linspace(0.01, 0.7)

fig = pl.figure(1)
pl.clf()

ax = fig.gca()

ax.plot(data, ratio_to_temperature(data), alpha=0.5, linewidth=2,
        label='LTE',
       )
ax.set_ylim(0,150)
ax.set_xlabel("Ratio H$_2$CO $3_{2,1}-2_{2,0} / 3_{0,3}-2_{0,2}$")
ax.set_ylabel("Temperature (K)")

fig.savefig(paths.fpath('h2co_temperature_vs_ratio_lte.png'))

pm = paraH2COmodel()

densind = np.argmin(np.abs(pm.densityarr[0,:,0]-4.5))
colind = np.argmin(np.abs(pm.columnarr[0,:,0]-13.5))

ax.plot(pm.modelratio1[:,densind,colind],
        pm.temparr[:,densind,colind],
        'r',
        alpha=0.5,
        linewidth=2,
        label='$n=10^{4.5}$ cm$^{-3}$, $N=10^{13.5}$ cm$^{-2}$',
       )

pl.legend(loc='best')
fig.savefig(paths.fpath('h2co_temperature_vs_ratio_lte_and_radex.png'))
