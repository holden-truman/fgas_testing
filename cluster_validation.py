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
from angle_calc import _compute_lensing_angles_flatsky

class ClusterLensing:
    def __init__(self, name, phot_cat_path, dgamma, z_cl, lensband, ra_c, dec_c):
        self.name = name
        self.phot_cat = fitsio.read(phot_cat_path)
        self.DGAMMA = dgamma
        self.z_cl = z_cl
        self.lensband = lensband
        self.RA_C = ra_c
        self.DEC_C = dec_c
        
        # Compute lensing angles
        self.r_rad, self.phi = compute_lensing_angles_astropy(self.RA_C, self.DEC_C, self.phot_cat["ra"], self.phot_cat["dec"])
        self.r_rad2, self.phi2 = _compute_lensing_angles_astropy2(self.RA_C, self.DEC_C, self.phot_cat["ra"], self.phot_cat["dec"], coordinate_system="celestial")
        self.r_rad2, self.phi2 = _compute_lensing_angles_astropy3(self.RA_C, self.DEC_C, self.phot_cat["ra"], self.phot_cat["dec"], coordinate_system="celestial")
    
    def plot_radius(self):
        plt.scatter(self.phi2,self.phi3)
        plt.title(self.name)
        plt.xlabel(r'$\phi_{old}$ [rad]')
        plt.ylabel(r'$\phi_{new}$ [rad]')

        plt.show()

    def plot_ra_dec(self):
        # Extract RA and DEC
        ra = self.phot_cat['ra']
        dec = self.phot_cat['dec']
        
        # Create the plot
        plt.figure(figsize=(8, 6))
        plt.scatter(ra, dec, s=1, color='blue', alpha=0.5)
        
        # Set labels and title
        plt.xlabel('Right Ascension (deg)')
        plt.ylabel('Declination (deg)')
        plt.title('RA vs DEC from FITS Catalog')
        
        # Show grid
        plt.grid(True)
        
        # Show the plot
        plt.show()

# Example usage:
rxj2129_lensing = ClusterLensing(name="RXJ2129", phot_cat_path='/sdf/home/h/holden/from_lucie_2024/RXJ2129/RXJ2129_run001_v2_adam.fits',
                           dgamma=2*0.01, z_cl=0.234, lensband='W-C-RC', ra_c=322.41650000, dec_c=0.10588889)

macs1115_lensing = ClusterLensing(name="MACS1115", phot_cat_path='/sdf/home/h/holden/from_lucie_2024/MACS1115/MACS1115_run002_v2_adam_2.fits',
                           dgamma=2*0.01, z_cl=0.355, lensband='W-C-RC', ra_c=168.96666667, dec_c=1.49861111)

rxj2129_lensing.plot_radius()
macs1115_lensing.plot_radius()

#rxj2129_lensing.plot_ra_dec()
#macs1115_lensing.plot_ra_dec()

'''
sk_lens = SkyCoord(10 * u.deg, 0 * u.deg, frame="icrs")
sk_src = SkyCoord(10 * u.deg, 2 * u.deg, frame="icrs")
angsep, phi = sk_lens.separation(sk_src).rad, sk_lens.position_angle(sk_src).rad

print(sk_lens)
print(sk_src)
print(phi)
'''
