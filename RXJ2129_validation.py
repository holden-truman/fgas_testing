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
from copy import copy
import matplotlib.pyplot as plt
import numpy.lib.recfunctions as rfn

import matplotlib.patches as mpatches
import scipy
from scipy.optimize import curve_fit

#load shear and photometry catalogs
phot_cat = fitsio.read('/sdf/home/h/holden/from_lucie_2024/RXJ2129/RXJ2129_run001_v2_adam.fits')
#shear_cat=fitsio.read('/sps/lsst/users/lbaumont/data/SUBARU/RXJ2129/lensing/ngmix/RXJ2129_lensing_run004.fits')


def compute_lensing_angles_astropy(ra_lens, dec_lens, ra_source_list, dec_source_list):
    #from clmm
    r"""Compute the angular separation between the lens and the source and the azimuthal
    angle from the lens to the source in radians.

    Parameters
    ----------
    ra_lens: float
        Right ascension of the lensing cluster in degrees
    dec_lens: float
        Declination of the lensing cluster in degrees
    ra_source_list: array
        Right ascensions of each source galaxy in degrees
    dec_source_list: array
        Declinations of each source galaxy in degrees

    Returns
    -------
    angsep: array
        Angular separation between the lens and the source in radians
    phi: array
        Azimuthal angle from the lens to the source in radians
    """
    if not -360. <= ra_lens <= 360.:
        raise ValueError(f"ra = {ra_lens} of lens if out of domain")
    if not -90. <= dec_lens <= 90.:
        raise ValueError(f"dec = {dec_lens} of lens if out of domain")
    if not all(-360. <= x_ <= 360. for x_ in ra_source_list):
        raise ValueError("Cluster has an invalid ra in source catalog")
    if not all(-90. <= x_ <= 90 for x_ in dec_source_list):
        raise ValueError("Cluster has an invalid dec in the source catalog")
    sk_lens = SkyCoord(ra_lens*u.deg, dec_lens*u.deg, frame='icrs')
    sk_src = SkyCoord(ra_source_list*u.deg,
                      dec_source_list*u.deg, frame='icrs')
    angsep, phi = sk_lens.separation(
        sk_src).rad, sk_lens.position_angle(sk_src).rad
    # Transformations for phi to have same orientation as _compute_lensing_angles_flatsky
    phi += 0.5*np.pi
    phi[phi > np.pi] -= 2*np.pi
    return angsep, phi



def _compute_lensing_angles_astropy2(
    ra_lens, dec_lens, ra_source_list, dec_source_list, coordinate_system="euclidean"
):
    r"""Compute the angular separation between the lens and the source and the azimuthal
    angle from the lens to the source in radians.

    Parameters
    ----------
    ra_lens: float
        Right ascension of the lensing cluster in degrees
    dec_lens: float
        Declination of the lensing cluster in degrees
    ra_source_list: array
        Right ascensions of each source galaxy in degrees
    dec_source_list: array
        Declinations of each source galaxy in degrees
    coordinate_system: str, optional
        Coordinate system of the ellipticity components. Must be either 'celestial' or
        euclidean'. See https://doi.org/10.48550/arXiv.1407.7676 section 5.1 for more details.
        Default is 'euclidean'.

    Returns
    -------
    angsep: array
        Angular separation between the lens and the source in radians
    phi: array
        Azimuthal angle from the lens to the source in radians
    """
    sk_lens = SkyCoord(ra_lens * u.deg, dec_lens * u.deg, frame="icrs")
    sk_src = SkyCoord(ra_source_list * u.deg, dec_source_list * u.deg, frame="icrs")
    angsep, phi = sk_lens.separation(sk_src).rad, sk_lens.position_angle(sk_src).rad
    # Transformations for phi to have same orientation as _compute_lensing_angles_flatsky
    phi += 0.5 * np.pi
    if np.iterable(phi):
        phi[phi > np.pi] -= 2 * np.pi
        phi[angsep == 0] = 0
    else:
        phi -= 2 * np.pi if phi > np.pi else 0
        phi = 0 if angsep == 0 else phi
    if coordinate_system == "celestial":
        phi = np.pi - phi
    return angsep, phi

r_rad,phi = compute_lensing_angles_astropy(RA_C, DEC_C, phot_cat["ra"], phot_cat["dec"])
r_rad2,phi2 = _compute_lensing_angles_astropy2(RA_C, DEC_C, phot_cat["ra"], phot_cat["dec"],coordinate_system="euclidean"