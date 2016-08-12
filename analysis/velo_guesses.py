from astropy.table import Table
import paths

myvtbl = Table.read(paths.tpath('core_velocities.txt'),
                    # write as ascii.fixed_width
                    format='ascii.fixed_width', delimiter='|')

guesses = {row['source']:row['velocity'] for row in myvtbl}
