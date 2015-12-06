from astroquery.alma import Alma

a = Alma()
a.login('keflavich')

# Find my project
result = a.query_object('W51', payload=dict(pi_name='ginsburg'), public=False)
_7m = result['Spatial resolution'] > 3
_12m = result['Spatial resolution'] < 3

# unfortunately, ALMA doesn't support querying for MOUSes
uid_url_table = a.stage_data(result['Member ous id'])

filelist = a.download_and_extract_files(uid_url_table['URL'],
                                        include_asdm=True, regex='.*',
                                        delete=True, verbose=True)
