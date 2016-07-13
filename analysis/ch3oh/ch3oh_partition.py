import numpy as np
import paths
import pylab as pl
from astroquery.vamdc import Vamdc
from vamdclib import specmodel

ch3oh = Vamdc.query_molecule('CH3OH')

temperatures = np.linspace(10,300)

partition_func = [specmodel.calculate_partitionfunction(ch3oh.data['States'],
                                                        temperature=tex)['XCDMS-149']
                  for tex in temperatures]


pl.matplotlib.rc_file('pubfiguresrc')
pl.figure(1).clf()
pl.plot(temperatures, partition_func)
pl.xlabel("T (K)")
pl.ylabel("$Q_{rot}$")
pl.savefig(paths.fpath("chemistry/ch3oh_partition_function.png"))
