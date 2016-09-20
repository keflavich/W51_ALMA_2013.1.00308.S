import os
import numpy as np
import paths
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import constants
from astropy.io import fits
from astropy import log
from line_to_image_list import line_to_image_list
import pyregion

log.warning("Seems that the fraction is too low; some lines should be 100% threecore"
            " but are not")

threecore_reg = pyregion.open(paths.rpath("three1ascores.reg"))
twelvem_ptgs = pyregion.open(paths.rpath('12m_pointings.reg'))
not_too_noisy = pyregion.open(paths.rpath('not_too_noisy_box.reg'))

for line, restfreq, velocity_res, spw in line_to_image_list:
    fn = paths.dpath('merge/moments/W51_b6_7M_12M.{0}.image.pbcor_medsub_moment0.fits'.format(line))
    #noisefn = paths.dpath('merge/moments/W51_b6_7M_12M.{0}.image.pbcor_medsub_madstd.fits'.format(line))
    if os.path.exists(fn):
        fh = fits.open(fn)
        # integrate over 15 km/s...
        #noise = fits.getdata(noisefn) * 15
        mask = threecore_reg.get_mask(fh[0])
        mask_12m = twelvem_ptgs.get_mask(fh[0])
        not_too_noisy_mask = not_too_noisy.get_mask(fh[0])
        pos = fh[0].data > np.nanstd(fh[0].data) * 2

        total_pos = fh[0].data[pos & mask_12m & not_too_noisy_mask].sum()
        total_threecore = fh[0].data[pos & mask & mask_12m & not_too_noisy_mask].sum()
        print("{0} has 3-core fraction = {1}".format(line, total_threecore/total_pos))
    else:
        print("Skipped {0}".format(line))

"""
H2CO303_202 has 3-core fraction = 0.058371385227015314
H2CO321_220 has 3-core fraction = 0.0751399004994529
H2CO322_221 has 3-core fraction = 0.05633021072527247
12CO2-1 has 3-core fraction = 0.015078962977107289
HC3N24-23 has 3-core fraction = 0.20604911069498635
HC3Nv7=124-23 has 3-core fraction = 0.6891917088368432
HC3Nv7=1_24-23 has 3-core fraction = 0.6653100531030567
Skipped HC3Nv7=2_24-23
OCS18-17 has 3-core fraction = 0.24658824270891316
OCS19-18 has 3-core fraction = 0.22620151022279547
SO65-54 has 3-core fraction = 0.09535002004924674
HNCO10110-919 has 3-core fraction = 0.6960424610898156
HNCO1028-927 has 3-core fraction = 0.7015274678037113
HNCO10010-909 has 3-core fraction = 0.3173586566020047
HNCO1055-954 has 3-core fraction = 0.7877860655254458
HNCO1046-945 has 3-core fraction = 0.2605853128887476
HNCO1038-937 has 3-core fraction = 0.80287297173315
HNCO28128-29029 has 3-core fraction = 0.01819461521342965
Skipped HNCO1019-918
HNCO1046-945 has 3-core fraction = 0.2605853128887476
HNCO1038-937 has 3-core fraction = 0.80287297173315
CH3OH422-312 has 3-core fraction = 0.06995864405207089
CH3OH423-514 has 3-core fraction = 0.21526706112038796
CH3OH5m42-6m43 has 3-core fraction = 0.2705469330280903
CH3OH808-716 has 3-core fraction = 0.1533925943145102
CH3OH1029-936 has 3-core fraction = 0.29473757520098254
CH3OH18315-17414 has 3-core fraction = 0.39680471636537057
CH3OH25322-24420 has 3-core fraction = 0.18415541722286188
CH3OH23519-22617 has 3-core fraction = 0.3198050225238298
/Users/adam/work/w51/alma/analysis/molecular_spatial_distribution.py:33: RuntimeWarning: invalid value encountered in double_scalars
13CH3OH515-414 has 3-core fraction = nan
Skipped CH3OH_10m38-11m210
Skipped CH3OH_18316-17413
Skipped CH3OH_514-422
Skipped CH3OH_3m22-220
Skipped CH3OH_10m55-11m48
13CS5-4 has 3-core fraction = 0.09273134444151157
PN5-4 has 3-core fraction = 0.06298571389458033
NH2CHO11210-1029 has 3-core fraction = 0.4966767016090414
NH2CHO1156-1055 has 3-core fraction = 0.7142704980674808
H30alpha has 3-core fraction = 0.05653341828339364
C18O2-1 has 3-core fraction = 0.05506450547094455
H2CCO11-10 has 3-core fraction = 0.14393541683796685
HCOOH431-524 has 3-core fraction = 0.6455287493204379
CH3OCHO17314-16313E has 3-core fraction = 0.37732119171479406
CH3OCHO17413-16412A has 3-core fraction = 0.40666350270932705
CH3OCHO17314-16313A has 3-core fraction = 0.38879233902978094
Skipped CH3OCHO1249-1138E
CH3CH2CN24321-23320 has 3-core fraction = 0.5676066081469056
CH3CH2CN24222-23221 has 3-core fraction = 0.40874714377888793
Acetone21120-20219AE has 3-core fraction = 0.7006274658445651
Acetone21120-20119EE has 3-core fraction = 0.35966030078337624
Skipped Acetone871-744AE
Skipped Acetone1294-1183EE
Skipped Acetone12112-11101AE
H213CO312-211 has 3-core fraction = 0.14749507232385622
CH3CH2OH550-541 has 3-core fraction = 0.5848653091461419
O13CS18-17 has 3-core fraction = 0.3221345008638999
CH3OCH323321-23222AA has 3-core fraction = 0.32472861872960557
CH3OCH313013-12112AA has 3-core fraction = 0.20750696403437224
CH3OCH323321-23222EE has 3-core fraction = 0.40989863243854646
N2D+_3-2 has 3-core fraction = 0.35475762208569334
Skipped 13CO2-1
Skipped CH3CNv8=1_123-113l-1
Skipped CH3CNv8=1_124-114l-1
Skipped CH3CNv8=1_123-113l+1
Skipped CH3CNv8=1_122-112l-1
Skipped CH3CNv8=1_125-115l+1
Skipped CH3CNv8=1_126-116l+1
Skipped CH3CNv8=1_127-117l+1
Skipped CH3CNv8=1_125-115l-1
Skipped NH2CHO1138-1037
Skipped NH2CHO1037-11110
Skipped NH2CHO12012-11111
Skipped NH2CHO1019-918
Skipped CH3CN_128-118
Skipped CH3CN_127-117
Skipped CH3CN_126-116
Skipped CH3CN_125-115
Skipped CH3CN_124-114
Skipped CH3CN_123-113
Skipped CH3CN_122-112
Skipped CH3OCHO_954-845E
Skipped CH3OCHO_990-981E
Skipped CH3OCHO_991-982A
Skipped SO2_422-313
Skipped SO2_28325-28226
SO2_16610-17513 has 3-core fraction = 0.5690400542434736
Skipped g-CH3CH2OH_13211-12210
Skipped g-CH3CH2OH_651-541
Skipped t-CH3CH2OH_16511-16412
Skipped t-CH3CH2OH_14014-13113
SO2_22715-23618 has 3-core fraction = 0.7109607931246041
SO2_16610-17513 has 3-core fraction = 0.5690400542434736
SO2v2=1_20218-19317 has 3-core fraction = 0.2809554596991058
SO2v2=1_22220-22121 has 3-core fraction = 0.6807151868461042
SO2v2=1_16313-16214 has 3-core fraction = 0.4690871788104719
SO2v2=1_642-735 has 3-core fraction = 0.00576215337260039
Skipped CH3NCO_25223-24222
Skipped CH3NCO_2500-2400
Skipped CH3NCO_2520-2420
CH3NCO_25124-24123 has 3-core fraction = 0.6741348590729396
/Users/adam/work/w51/alma/analysis/molecular_spatial_distribution.py:33: RuntimeWarning: invalid value encountered in double_scalars
CH3NCO_27226-26225 has 3-core fraction = nan
Skipped CH3SH_130-121
CH3SH_152-151 has 3-core fraction = 0.5840554755086186
CH3SH_162-161 has 3-core fraction = 0.5908614401559547
CH3SH_73-82 has 3-core fraction = 0.1879621230542317
CH3SH_232-231 has 3-core fraction = 0.660239351834896
Skipped SiO_5-4
"""
