raise("spectral_extraction is meant for bigger machines...")
def execfile(fn):
    with open(fn, 'r') as fh:
        contents = fh.read()
    exec(contents)

execfile('spectral_extraction.py')
execfile('core_velocities.py')
execfile('core_photometry.py')
execfile('outflow_photometry.py')
execfile('merge_outflow_core_tables.py')
