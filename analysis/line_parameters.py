"""
Lines to examine as possible brightest lines in bad for determining the source
velocity.  This is a more restrictive list than the lines_to_image_list; it
excludes things that are likely to be wrong like h30a and C18O (which will
trace ionized/diffuse gas)
"""
from astropy import units as u

#lines_to_overlay = ['OCS','H2CO', 'HNCO', 'SO']
frequencies = {'H2CO303_202': 218.22219*u.GHz,
               'H2CO321_220': 218.76007*u.GHz,
               'H2CO322_221': 218.47564*u.GHz,
               'CH3OH422-312': 218.44005*u.GHz,
               'HC3N24-23': 218.32471*u.GHz,
               'OCS18-17': 218.90336*u.GHz,
               'OCS19-18': 231.06099*u.GHz,
               'SO65-54': 219.94944*u.GHz,
               'HNCO10110-919': 218.98102*u.GHz,
               'HNCO1028-927': 219.73719*u.GHz,
               'CH3OH423-514': 234.68345*u.GHz,
               'CH3OH5m42-6m43': 234.69847*u.GHz,
               'CH3OH808-716': 220.07849*u.GHz,
               '13CS5-4': 231.22069*u.GHz,
               'CH3OCH3_13013-12112': 231.98792*u.GHz,
               'NH2CHO11210-1029': 232.27363*u.GHz,
               'NH2CHO1156-1055': 233.59451*u.GHz,
               'HC3Nv7=124-23': 219.17358*u.GHz,
               # offset too much 'N2D+_3-2': 231.32183*u.GHz,
               #'H30alpha': 231.90093*u.GHz,
               #'C18O2-1': 219.56036*u.GHz,
              }
freq_name_mapping = {v:k for k,v in frequencies.items()}
yoffset = {'H2CO303_202': 0,
           'H2CO321_220': 1,
           'H2CO322_221': 2,
           'OCS18-17': 3,
           'SO65-54': 4,
           'CH3OH423-514': 5,
           'CH3OH5m42-6m43': 6,
           'OCS19-18': 7,
           '13CS5-4': 8,
           'CH3OCH3_13013-12112': 9,
           'HNCO1028-927': 10,
           'HNCO10110-919': 11,
           'HC3N24-23': 12,
           'HC3Nv7=124-23': 13,
           'NH2CHO11210-1029': 14,
           'NH2CHO1156-1055': 15,
           'CH3OH422-312': 16,
           'CH3OH808-716': 17,
           # offset too much 'N2D+_3-2': 18,
           #'H30alpha': 4.5,
           #'C18O2-1': 3.5,

          }
