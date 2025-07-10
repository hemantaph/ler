import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
# from scipy.interpolate import interp1d
from .jit_functions import sfr_madau_fragos2017

def sfr_with_time_delay(input_args):
    """
    Compute the star formation rate at redshift z, given parameters a, b, c, and d, 
    and cosmological parameters H0, Omega_M, and Omega_Lambda. 
    The star formation rate is time-delayed relative to the observed redshift, 
    with a time delay uniformly distributed between td_min and td_max. 
    The time delay is computed using the cosmology provided by astropy.

    Parameters
    ----------
    input_args : list
        z : float
            observed redshift
        idx : int
            index of the galaxy
        td_min : float
            minimum time delay in Gyr
        td_max : float
            maximum time delay in Gyr
        H0 : float
            Hubble constant in km/s/Mpc
        Omega_M : float
            matter density parameter
        Omega_Lambda : float
            dark energy density parameter
        a : float
            parameter of the Madau-Fragos star formation rate
        b : float
            parameter of the Madau-Fragos star formation rate
        c : float
            parameter of the Madau-Fragos star formation rate
        d : float
            parameter of the Madau-Fragos star formation rate

    Returns
    -------
    idx : int
        index of the galaxy
    result : float
        star formation rate at observed redshift z
    """
    z = input_args[0]
    idx = input_args[1]
    td_min = input_args[2]
    td_max = input_args[3]
    H0 = input_args[4]
    Omega_M = input_args[5]
    Omega_Lambda = input_args[6]
    a = input_args[7]
    b = input_args[8]
    c = input_args[9]
    d = input_args[10]

    def E(z_prime):
        """Hubble parameter as a function of redshift."""
        return np.sqrt(Omega_M * (1 + z_prime)**3 + Omega_Lambda)

    def integrand(z_prime):
        """Integrand for the equation."""
        return 1 / (H0* (1 + z_prime) * E(z_prime)) * 977.813 

    def time_delay(zform, z):
        """Time delay between formation and observation."""
        integral, _ = quad(integrand, z, zform)
        return integral

    def equation_to_solve(zform, z, td):
        """Equation to solve for zform."""
        return td - time_delay(zform, z)
    
    def find_zform(z, td):
        """Find zform using grid search and linear interpolation."""
        zform_solution = fsolve(equation_to_solve, z, args=(z, td))
        return zform_solution
    
    # montecalo integration
    def integrand_rates(z, size=100000, zform_max=1000.):

        td = np.random.uniform(td_min, td_max, size)
        td_max_allowed = time_delay(zform_max, z)
        # idx = td < td_max_allowed
        idx = np.where(td < td_max_allowed)[0]
        P_td = np.zeros_like(td)
        # print('idx', idx)
        # print('td', td)
        # print('td_max_allowed', td_max_allowed)
        P_td[idx] = 1/(np.log(td_max/td_min) * td[idx]) 

        zform = np.zeros_like(td)
        for idx_ in idx:
            zform[idx_] = find_zform(z, td[idx_])

        psi = sfr_madau_fragos2017(zform, a, b, c, d)

        integral = 1/(td_max - td_min) * np.sum(P_td * psi)
        return integral

    result = integrand_rates(z)
    return int(idx), result