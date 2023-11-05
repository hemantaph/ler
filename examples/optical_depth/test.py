# https://git.ligo.org/narola.bharatbhai/gw-lensing/-/blob/main/multi_messenger_tools/sh/tabulated_lum_functions.py

from scipy.special import gamma
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
from scipy.misc import derivative
import astropy.units as u
from astropy import units as u, constants as c
import numpy as np

def phi(s, z):
    return s**4 * philocbernardi(s) * pdf_phi_z_div_0(s, z) * pdfz(z)

def philocbernardi(sigma):
    # Very different sigmastars depending on sample selection.
    alpha_sch = 0.94
    beta_sch = 1.85
    phistar = 2.099e-2 * (cosmo.h / 0.7) ** 3  # Mpc**-3
    sigmastar = 113.78
    philoc_ = phistar*(sigma/sigmastar)**alpha_sch * np.exp(-(sigma/sigmastar)**beta_sch) * beta_sch/gamma(alpha_sch/beta_sch)/sigma
    return philoc_

def pdf_phi_z_div_0(s, z):
    phi_sim_z = 10 ** cvdf_fit(np.log10(s), z) / s * derivative(lambda ss: cvdf_fit(ss, z), np.log10(s), dx=1e-8)
    phi_sim_0 = 10 ** cvdf_fit(np.log10(s), 0) / s * derivative(lambda ss: cvdf_fit(ss, 0), np.log10(s), dx=1e-8)
    # phi_loc_0 = derivative(lambda ss: num_d.cvdf_fit(np.log10(ss),0), s, dx=1e-8)
    return phi_sim_z / phi_sim_0

def pdfz(z):
    # Over whole sky, so equal to dV/dz Mpc^3
    Dl = cosmo.angular_diameter_distance(z).to_value(u.Mpc)
    Dh = (c.c / cosmo.H0).to_value(u.Mpc)
    return 4 * np.pi * Dh * ((1 + z) ** 2 * Dl ** 2 / cosmo.efunc(z))

def cvdf_fit(log_vd, redshift):
    # from https://github.com/ptorrey/torrey_cmf/blob/master/torrey_cmf.py
    this_vars = np.array([
        [7.39149763, 5.72940031, -1.12055245],
        [-6.86339338, -5.27327109, 1.10411386],
        [2.85208259, 1.25569600, -0.28663846],
        [0.06703215, -0.04868317, 0.00764841]])
    coeffs = [this_vars[i][0] + this_vars[i][1] * redshift + this_vars[i][2] * redshift ** 2 for i in range(4)]
    mstar = log_vd - coeffs[3]
    return coeffs[0] + coeffs[1] * mstar + coeffs[2] * mstar ** 2 - np.exp(mstar)

def conditional_rejection_sample(pdf, xmin, xmax, y_fixed, size=100, chunk_size=10000):
    
    chunk_size = 10000
    
    x = np.linspace(xmin, xmax, chunk_size)
    y = np.ones(chunk_size) * y_fixed
    z = pdf(x,y)
    zmax = np.max(z)
    
    
    # Rejection sample in chunks
    x_sample = []
    while len(x_sample) < size:
        x_try = np.random.uniform(xmin, xmax, size=chunk_size)
        y_try = np.ones(chunk_size) * y_fixed
        
        z_try = np.random.uniform(0, zmax, size=chunk_size)
        zmax = max(zmax, np.max(z_try))

        x_sample += list(x_try[z_try < pdf(x_try, y_try)])
        
    # Transform the samples to a 1D numpy array
    x_sample = np.array(x_sample).flatten()
    # Return the correct number of samples
    return x_sample[:size]
