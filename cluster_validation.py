import numpy as np
import fitsio
import matplotlib.pyplot as plt
import matplotlib.ticker
import sys
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

# for cosomology dependence
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo 
from astropy.coordinates import SkyCoord
from copy import copy
import matplotlib.pyplot as plt
import numpy.lib.recfunctions as rfn

import matplotlib.patches as mpatches
import scipy
from scipy.optimize import curve_fit

from angle_calc import compute_lensing_angles_astropy
from angle_calc import _compute_lensing_angles_astropy2

class ClusterLensing:
    def __init__(self, phot_cat_path, dgamma, z_cl, lensband, ra_c, dec_c):
        self.phot_cat = fitsio.read(phot_cat_path)
        self.DGAMMA = dgamma
        self.z_cl = z_cl
        self.lensband = lensband
        self.RA_C = ra_c
        self.DEC_C = dec_c
        
        # Compute lensing angles
        self.r_rad, self.phi = compute_lensing_angles_astropy(self.RA_C, self.DEC_C, self.phot_cat["ra"], self.phot_cat["dec"])
        self.r_rad2, self.phi2 = _compute_lensing_angles_astropy2(self.RA_C, self.DEC_C, self.phot_cat["ra"], self.phot_cat["dec"], coordinate_system="euclidean")
    
    def plot_radius(self):
        plt.scatter(phi,phi2)
        plt.title("RXJ2129")
        plt.xlabel(r'$\phi_{old}$ [rad]')
        plt.ylabel(r'$\phi_{new}$ [rad]')

        plt.show()


# Example usage:
rxj2129_lensing = ClusterLensing(phot_cat_path='/sdf/home/h/holden/from_lucie_2024/RXJ2129/RXJ2129_run001_v2_adam.fits',
                           dgamma=2*0.01, z_cl=0.234, lensband='W-C-RC', ra_c=322.41650000, dec_c=0.10588889)

macs1115_lensing = ClusterLensing(phot_cat_path='/sdf/home/h/holden/from_lucie_2024/MACS1115/MACS1115_run002_v2_adam_2.fits',
                           dgamma=2*0.01, z_cl=0.355, lensband='W-C-RC', ra_c=168.96666667, dec_c=1.49861111)

rxj2129_lensing.plot_radius()
macs1115_lensing.plot_radius()