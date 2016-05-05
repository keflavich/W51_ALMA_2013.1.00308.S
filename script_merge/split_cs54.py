import numpy as np


restfrq = 231.221 # GHz
# vcen = 55 #km/s


ckms = 2.997924580e5
vstart = 30 #km/s
vend = 80 # km/s
vwidth = 1.265 #km/s
nchan = int(np.ceil((vend-vstart)/vwidth))
frqstart = restfrq * (1-vend/ckms)
frqwidth = vwidth / ckms * restfrq

cvel(vis='w51_concat_7m12m.spw2.merge',
     outputvis='w51_concat_7m12m_13cs54.ms',
     field='w51',
     passall=False,
     spw='',
     mode='frequency',
     start='{0}GHz'.format(frqstart),
     width='{0}GHz'.format(frqwidth),
     nchan=nchan, restfreq='{0}GHz'.format(restfrq), outframe='LSRK',
     veltype='radio',
     interpolation='linear',
     hanning=False,
     )
