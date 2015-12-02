import numpy as np
from calibrated_configuration import field, phasecenter, cell, imsize, weighting, robust, threshold, spws_12m, spws_7m, nchans_total, frange, fstep


ncubes_per_window = 20


for spwnum in '1320':
    spwnum = int(spwnum)

    concatvis = 'w51_concat_7m12m.spw{0}.merge'.format(spwnum)
    if not os.path.exists(concatvis):
        print("{0} does not exist - run cvel_and_merge first".format(concatvis))
        break

    print "# running clean on all lines in spw{0}".format(spwnum)
    nchans_total_thiscube = nchans_total[spwnum]
    nchans_per_cube = int(nchans_total_thiscube/ncubes_per_window)
    for ii in range(ncubes_per_window):
        start = nchans_per_cube*ii
        end = nchans_per_cube*(ii+1)
        if end > nchans_total_thiscube:
            end = nchans_total_thiscube
        output = 'piece_of_full_W51_7m12m_cube.spw{0}.channels{1}to{2}'.format(spwnum, start, end)

        # Channel-based gridding has major bugs when dealing with CVEL'd data
        # It is therefore necessary to compute the frequency gridding by hand
        startfreq = "{0}GHz".format(frange[spwnum][0]/1e3 + start * fstep[spwnum]/1e6)
        width = "{0}kHz".format(fstep[spwnum])


        # LINE IMAGING (MOSAIC MODE)
        if not os.path.exists(output+".image"):
            print "Imaging {0}".format(output)
            os.system('rm -rf ' + output + '*')
            clean(vis = concatvis,
                  imagename = output,
                  field = '',
                  spw = '', # there should be only one
                  imagermode = 'mosaic',
                  mode = 'frequency',
                  width = width,
                  start = startfreq,
                  nchan = nchans_per_cube,
                  # it is not clear whether chaniter=True causes problems,
                  # but it seems that chaniter=False is faster and works pretty
                  # well.
                  # Followup: in cubes with bright lines, chaniter=False may cause
                  # overall dips (non-convergent clean, maybe?)
                  chaniter = True,
                  veltype = 'radio',
                  outframe = 'LSRK',
                  interactive = F,
                  niter = 5000,
                  imsize = imsize,
                  cell = cell,
                  psfmode='clark',
                  weighting = weighting,
                  phasecenter = phasecenter,
                  robust = robust,
                  threshold = threshold,
                  pbcor = F,
                  # Do not store the models on disk.  First, it takes time,
                  # second, we don't care about them!
                  usescratch=F)

          
        myimagebase = output
        # I've given up on primary beam correction, at least for now
        exportfits(imagename=myimagebase+'.image',
                   fitsimage=myimagebase+'.image.fits',
                   overwrite=True,
                   dropdeg=True)

