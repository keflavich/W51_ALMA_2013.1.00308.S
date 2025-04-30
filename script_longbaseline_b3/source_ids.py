"""
A dictionary of source names and their locations
"""

from astropy import coordinates, units as u

sources = {#'w51e2e': (290.9331907, 14.50958333),
           'w51e8': (290.9329337, 14.50784566),
           #'w51e2se': (290.9336503, 14.50930947),
           #'w51e2w': (290.93295, 14.509592),
           #'w51e2nw': (290.93282, 14.510005),
           'w51n': (290.9168748, 14.51818126),
           #'w51d2': (290.9159067, 14.51800749),
           'w51n_mm24': (290.9164773, 14.51814722),
           #'w51n_123': (290.9171007, 14.51825111),
           'w51e5': (290.92426433, 14.51525949),
           'w51sw': (290.92427542, 14.51522778),
           'w51e_n': (290.93291768, 14.51648581),
           'w51main': (290.92400257, 14.50718514),
          }

source_field_mapping = {'w51e2e': 'w51e2',
                        'w51e8': 'w51e2',
                        'w51e2se': 'w51e2',
                        'w51e2w': 'w51e2',
                        'w51e2nw': 'w51e2',
                        'w51n': 'w51n',
                        'w51d2': 'w51n',
                        'w51n_mm24': 'w51n',
                        'w51n_123': 'w51n',
                        'w51e5': 'w51n',
                        'w51sw': 'w51n',
                        'w51e_n': 'w51e2',
                        'w51main': 'w51e2',
                       }

def casafmt(crd):
    center = coordinates.SkyCoord(*crd, frame='icrs', unit=(u.deg, u.deg))

    return 'ICRS {} {}'.format(center.ra.to_string(unit=u.hour),
                               center.dec.to_string())

sources_fmtd = {key: casafmt(val) for key, val in sources.items()}
