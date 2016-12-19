import numpy as np
import requests

def get_dust_opacity():
    rslt = requests.get('https://hera.ph1.uni-koeln.de/~ossk/Jena/tables/mrn5')
    with open('dustkappa_mrn5.inp','w') as fh:
        lines = rslt.content.rstrip().split("\n")
        wav,opac = np.array([list(map(float, line.split())) for line in lines]).T
        beta = np.log(opac[-1]/opac[-2])/np.log(wav[-2]/wav[-1])
        beta = 1.5
        lastwav = 4000.
        const = opac[-1] * wav[-1]**beta
        lastopac = const * lastwav**-beta
        fh.write("{0:10d}\n".format(1))
        fh.write("{0:10d}\n".format(len(lines)+1))
        for line in lines:
            fh.write("{0}\n".format(line))
        fh.write(" {0:9.3e} {1:9.3e}\n".format(lastwav,lastopac))

    dust_type = 'mrn5'

    with open('dustopac.inp', 'w') as fh:
        fh.write("2               Format number of this file\n")
        fh.write("1               Nr of dust species\n")
        fh.write("============================================================================\n")
        fh.write("1               Way in which this dust species is read\n")
        fh.write("0               0=Thermal grain, 1=Quantum heated\n")
        fh.write("{0}      Extension of name of dustkappa_***.inp file\n".format(dust_type))
        fh.write("----------------------------------------------------------------------------\n")
