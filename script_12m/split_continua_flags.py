"""
This appears broken:
2015-12-10 16:59:17     WARN    flagmanager:::: Version name 'before_cont_flags' already exist. Will rename it to before_cont_flags.old.1449766757
# of flagged points before continuum flagging: 20091340800.0
# of flagged points after continuum flagging: 20091340800.0
....10....20....30....40....50....60....70....80....90....100%
# of flagged points after restoration: 20091340800.0
....10....20....30....40....50....60....70....80....90....100%
2015-12-10 17:29:53     WARN    flagmanager:::: Version name 'before_cont_flags' already exist. Will rename it to before_cont_flags.old.1449768593
"""

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
    print("# of flagged points before continuum flagging: {flagged}/{total}".format(**summary))
    flagdata(vis=finalvis,mode='manual',
             spw=linechannels,flagbackup=False)
    summary = flagdata(vis=finalvis, mode='summary')
    print("# of flagged points before continuum flagging: {flagged}/{total}".format(**summary))

    split(vis=finalvis,
          spw=spws_12m[spw],
          field='w51',
          outputvis=contvis,
          width=[192,192],
          datacolumn='data')


    flagmanager(vis=finalvis,mode='restore',
                versionname='before_cont_flags')

    summary = flagdata(vis=finalvis, mode='summary')
    print("# of flagged points before continuum flagging: {flagged}/{total}".format(**summary))

    split(vis=finalvis,
          spw=spws_12m[spw],
          field='w51',
          outputvis='w51_spw{0}_continuum_noflag.split'.format(spw),
          width=[192,192],
          datacolumn='data')
