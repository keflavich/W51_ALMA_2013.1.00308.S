fixvis(vis='w51_spw3_continuum_flagged.split',
       outputvis='w51_spw3_continuum_fixvis.ms',
       field=','.join([str(x-4) for x in (24,25,31,32,33,39,40)]),
       datacolumn='data',
       phasecenter='J2000 19h23m43.896 +14d30m28.26',
     )

plotms(vis='w51_spw3_continuum_fixvis.ms', gridrows=4, gridcols=4, xaxis='time',
       yaxis='phase', timerange='2015/04/20~2015/04/28', avgchannel='1000',
       iteraxis='baseline', coloraxis='field', uvrange='75~10000', field='20,21,27,28,29,35,36')
