
spws_12m = {0: '0,4',
            1: '1,5',
            2: '2,6',
            3: '3,7',
           }

for spw in range(4):
    with open('linechannels12m_spw{0}'.format(spw),'r') as f:
        linechannels = f.read()


    finalvis='w51_concat.ms.split.cal'
    contvis='w51_spw{0}_continuum_flagged.split'.format(spw)
    flagmanager(vis=finalvis,mode='save',
                versionname='before_cont_flags')

    summary = flagdata(vis=finalvis, mode='summary')
    print("# of flagged points before continuum flagging: {0}".format(summary['total']))
    flagdata(vis=finalvis,mode='manual',
             spw=linechannels,flagbackup=False)
    summary = flagdata(vis=finalvis, mode='summary')
    print("# of flagged points after continuum flagging: {0}".format(summary['total']))

    split(vis=finalvis,
          spw=spws_12m[spw],
          field='w51',
          outputvis=contvis,
          width=[192,192],
          datacolumn='data')


    flagmanager(vis=finalvis,mode='restore',
                versionname='before_cont_flags')

    summary = flagdata(vis=finalvis, mode='summary')
    print("# of flagged points after restoration: {0}".format(summary['total']))

    split(vis=finalvis,
          spw=spws_12m[spw],
          field='w51',
          outputvis='w51_spw{0}_continuum_noflag.split'.format(spw),
          width=[192,192],
          datacolumn='data')
