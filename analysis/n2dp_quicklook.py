import glob,pyspeckit
for fn in glob.glob("*spw2*fits"):
    sp = pyspeckit.Spectrum(fn)
    sp.xarr.convert_to_unit(u.km/u.s, refX=231.32222428571433*u.GHz)
    x = sp.xarr
    x0 = x[0].value
    sp.plotter(xmin=20, xmax=90, figure=pl.figure(1))
    sp.xarr.convert_to_unit(u.GHz)
    sp.xarr.convert_to_unit(u.km/u.s, refX=231.22069*u.GHz)
    assert x0 != sp.xarr[0].value
    sp.plotter(xmin=20, xmax=90, figure=pl.figure(1), clear=False, color='r')
    sp.plotter.axis.set_xlabel("Velocity (km/s)")
    sp.plotter.savefig("n2dp/{0}".format(fn.replace(".fits","_n2dp.png")))
