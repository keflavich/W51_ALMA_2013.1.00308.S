import numpy as np
import aplpy
import paths
import matplotlib
matplotlib.rc_file('pubfiguresrc')


matplotlib.pyplot.figure(1).clf()
F = aplpy.FITSFigure(paths.dpath('chemslices/chemical_m0_slabs_e2_CH3OH1029-936_merge.fits'),
                     figure=matplotlib.pyplot.figure(1))
F.show_grayscale(invert=True, vmax=1100, vmin=-20)
F.show_contour(paths.dpath('W51_te_continuum_best.fits'),
               colors=['r']*11,
               levels=[0.015, 0.0256944, 0.0577778, 0.11125, 0.186111,
                       0.282361, 0.4, ])
F.save(paths.fpath("continuum_contours_on_ch3oh1029.png"))

matplotlib.pyplot.figure(2).clf()
F = aplpy.FITSFigure(paths.dpath('W51_te_continuum_best.fits'),
                     figure=matplotlib.pyplot.figure(2))
F.recenter(290.93315, 14.5097, 2.8/3600.)
F.show_grayscale(invert=True, vmax=0.43, vmin=-0.01, stretch='arcsinh')
F.show_contour(paths.dpath('chemslices/chemical_m0_slabs_e2_CH3OH1029-936_merge.fits'),
               colors=['b']*11,
               levels=np.linspace(200,1100,6), layer='CH3OH')
#F.ticks.set_xspacing(0.0005)
F.save(paths.fpath("ch3oh1029_contours_on_continuum.png"))

matplotlib.pyplot.figure(3).clf()
F = aplpy.FITSFigure(paths.dpath('12m/moments/CH3OH_e2e8_cutout_temperaturemap.fits'),
                     figure=matplotlib.pyplot.figure(3))
F.show_colorscale(vmax=600, vmin=0, cmap='hot')
F.show_contour(paths.dpath('W51_te_continuum_best.fits'),
               colors=['r']*11,
               levels=[0.015, 0.0256944, 0.0577778, 0.11125, 0.186111,
                       0.282361, 0.4, ])
F.save(paths.fpath("chemistry/continuum_contours_on_ch3oh_temperature.png"))
