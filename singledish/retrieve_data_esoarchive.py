"""
Download files for our project from the ESO archive using astroquery
"""
from astroquery.eso import Eso

def retrieve_ESO_files(rawpath='/Volumes/passport/w51-apex/raw/',
                       username='aginsburg', projids=['098.C-0421',],
                       download=False):

    # don't touch the rest
    Eso.cache_location = rawpath

    Eso.login(username)

    Eso.ROW_LIMIT = 1000000
    #tbl = Eso.query_instrument('apex', pi_coi='ginsburg', cache=False)
    #stbl = tbl[np.char.startswith(tbl['Object'], 'Map') & (tbl['Scan Mode']=='OTF') & (tbl['Number of subscans'] > 10)]
    #programs = set(tbl['ProgId'])
    #projids = set(tbl['APEX Project ID'])

    all_files = []
    for proj in projids:
        tbl = Eso.query_apex_quicklooks(proj, cache=False)
        print(tbl)

        if download:
            files = Eso.retrieve_data(tbl['Product ID'], continuation=True)
            all_files.append(files)

    return all_files
