split(vis='w51_concat.ms.split.cal',
      outputvis='w51_test_small.ms',
      field='31,32,33,39,40,24,25',
      spw='3:1000~1200,7:1000~1200',
      datacolumn='data',
     )

cvel(vis='w51_concat.ms.split.cal',
     outputvis='w51_test_small.ms',
     field='31,32,33,39,40,24,25',
     spw='3,7',
     mode='frequency',
     nchan=100,
     start='233540.032MHz',
     width='0.488281MHz',
     )
