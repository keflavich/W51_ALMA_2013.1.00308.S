import numpy as np
def mass_conversion_factor(TK = 20, d=5.4):
    return 14.30 * (np.exp(13.01/TK) - 1)*d**2
def col_conversion_factor(TK = 20):
    return 2.19e22 * (np.exp(13.01/TK - 1))

