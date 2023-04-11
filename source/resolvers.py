import numpy as np

def deltaR(rap1, rap2, phi1, phi2):
    return np.sqrt(np.power(rap1 - rap2, 2) + np.power(phi1 - phi2, 2))