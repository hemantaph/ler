import numpy as np
from numba import njit

from .prior_functions import powerlaw_pdf, _gaussian_pdf, _smoothing_S
from ..utils import inverse_transform_sampler, cumulative_trapezoid

# ------------------------------
# powerlaw + gaussian
# ------------------------------
@njit(fastmath=True)
def powerlaw_plus_peak_unnormalized(m, mminbh=4.98, mmaxbh=112.5, alpha=3.78, mu_g=32.27, sigma_g=3.88, lambda_peak=0.003, delta_m=4.8):
    """
        Compute unnormalized PDF for power-law + Gaussian model.

        Evaluates the unnormalized mixture of power-law and Gaussian distributions
        with low-mass smoothing applied, prior to normalization.

        Parameters
        ----------
        m : ``numpy.ndarray``
            Mass values (solar masses).
        mminbh : ``float``
            Minimum BH mass (solar masses).
        mmaxbh : ``float``
            Maximum BH mass (solar masses).
        alpha : ``float``
            Power-law spectral index.
        mu_g : ``float``
            Mean of Gaussian peak (solar masses).
        sigma_g : ``float``
            Standard deviation of Gaussian peak (solar masses).
        lambda_peak : ``float``
            Fraction in Gaussian component (0-1).
        delta_m : ``float``
            Smoothing width (solar masses).

        Returns
        -------
        density : ``numpy.ndarray``
            Unnormalized probability density values.
    """

    pdf_unnormalized = ((1-lambda_peak)*powerlaw_pdf(m, alpha, mminbh, mmaxbh) + (lambda_peak * _gaussian_pdf(m, mu_g, sigma_g)))* _smoothing_S(m, mminbh, delta_m)

    return pdf_unnormalized

@njit(fastmath=True)
def powerlaw_plus_peak_pdf(m, mminbh=4.98, mmaxbh=112.5, alpha=3.78, mu_g=32.27, sigma_g=3.88, lambda_peak=0.003, delta_m=4.8, normalization_size=500):
    """
        Evaluate fully normalized primary-mass probability density (power-law + Gaussian model).

        Computes the normalized PDF combining a power-law distribution with a single
        Gaussian peak for black hole primary masses.

        Parameters
        ----------
        m : ``numpy.ndarray``
            Mass values to evaluate (solar masses).
        mminbh : ``float``
            Minimum BH mass (solar masses).
        mmaxbh : ``float``
            Maximum BH mass (solar masses).
        alpha : ``float``
            Power-law spectral index for the power-law component.
        mu_g : ``float``
            Mean of the Gaussian peak (solar masses).
        sigma_g : ``float``
            Standard deviation of the Gaussian peak (solar masses).
        lambda_peak : ``float``
            Fraction of sources in the Gaussian component (0-1).
        delta_m : ``float``
            Low-mass smoothing width (solar masses).
        normalization_size : ``int``
            Grid size for normalization. 
            default: 500

        Returns
        -------
        pdf : ``numpy.ndarray``
            Normalized probability density values.

        Examples
        --------
        >>> from ler.gw_source_population.prior_functions import powerlaw_plus_peak_pdf
        >>> import numpy as np
        >>> m = np.array([20.0, 40.0, 60.0])
        >>> pdf_values = powerlaw_plus_peak_pdf(m, mminbh=5, mmaxbh=100, alpha=2.3, mu_g=50, sigma_g=10, lambda_peak=0.1, delta_m=5)
        >>> print(pdf_values)  # doctest: +SKIP
    """

    m1_grid = np.geomspace(mminbh, mmaxbh, normalization_size)
    pdf_unnormalized = powerlaw_plus_peak_unnormalized(
        m1_grid, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m
    )
    _, _, normalization = cumulative_trapezoid(y=pdf_unnormalized, x=m1_grid, initial=0.0)

    pdf_m = powerlaw_plus_peak_unnormalized(
        m, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m
    )

    return pdf_m / normalization

@njit(fastmath=True)
def powerlaw_plus_peak_cdf(size, mminbh=4.98, mmaxbh=112.5, alpha=3.78, mu_g=32.27, sigma_g=3.88, lambda_peak=0.003, delta_m=4.8):
    """
        Compute cumulative distribution function for power-law + Gaussian mass model.

        Evaluates the CDF using numerical integration for primary mass distribution.

        Parameters
        ----------
        size : ``int``
            Number of grid points for CDF computation.
        mminbh : ``float``
            Minimum BH mass (solar masses).
        mmaxbh : ``float``
            Maximum BH mass (solar masses).
        alpha : ``float``
            Power-law spectral index.
        mu_g : ``float``
            Mean of the Gaussian peak (solar masses).
        sigma_g : ``float``
            Standard deviation of the Gaussian peak (solar masses).
        lambda_peak : ``float``
            Fraction in Gaussian component (0-1).
        delta_m : ``float``
            Low-mass smoothing width (solar masses).

        Returns
        -------
        cdf_values : ``numpy.ndarray``
            Normalized CDF values at each grid point.
    """

    m_try = np.geomspace(mminbh, mmaxbh, size)
    pdf_unnormalized = powerlaw_plus_peak_unnormalized(
        m_try, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m
    )
    cdf_values, m_try, _ = cumulative_trapezoid(y=pdf_unnormalized, x=m_try, initial=0.0)

    return cdf_values, m_try

@njit(fastmath=True)
def powerlaw_plus_peak_function(m, mminbh=4.98, mmaxbh=112.5, alpha=3.78, mu_g=32.27, sigma_g=3.88, lambda_peak=0.003, delta_m=4.8, normalization_size=500, R0=16.158, kappa=3.166, z_eval=0.2):
    """
        Evaluate the redshift-dependent merger rate density for the power-law + Gaussian mass model.
        Combines the normalized mass distribution with a redshift evolution model to compute the merger rate density as a function of mass and redshift.

        Parameters
        ----------
        m : ``numpy.ndarray``
            Mass values to evaluate (solar masses).
        mminbh : ``float``
            Minimum BH mass (solar masses).
        mmaxbh : ``float``  
            Maximum BH mass (solar masses).
        alpha : ``float``
            Power-law spectral index for the mass distribution.
        mu_g : ``float``
            Mean of the Gaussian peak (solar masses).
        sigma_g : ``float``
            Standard deviation of the Gaussian peak (solar masses).
        lambda_peak : ``float``
            Fraction of sources in the Gaussian component (0-1).
        delta_m : ``float``
            Low-mass smoothing width (solar masses).
        normalization_size : ``int``
            Grid size for normalization of the mass distribution. 
            default: 500
        R0 : ``float``
            Local merger rate density at z=0 (Gpc^-3 yr^-1). 
            default: 16.158
        kappa : ``float``   
            Power-law index for redshift evolution of the merger rate. 
            default: 3.166  
        z_eval : ``float``
            Redshift at which to evaluate the merger rate density.
            default: 0.2

        Returns
        -------
        rate_density : ``numpy.ndarray``
            Merger rate density as a function of mass at the specified redshift (Gpc^-3 yr^-1 Msun^-1).
    """

    # get pdf value
    pdf_m1 = powerlaw_plus_peak_pdf(m, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, normalization_size)

    # rate
    R_z = R0 * (1.0 + z_eval) ** kappa

    return R_z * pdf_m1

@njit(fastmath=True)
def powerlaw_plus_peak_rvs(size, mminbh=4.98, mmaxbh=112.5, alpha=3.78, mu_g=32.27, sigma_g=3.88, lambda_peak=0.003, delta_m=4.8, normalization_size=500):
    """
        Generate random samples from the power-law + Gaussian mass distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        mminbh : ``float``
            Minimum BH mass (solar masses).
        mmaxbh : ``float``
            Maximum BH mass (solar masses).
        alpha : ``float``
            Power-law spectral index for the power-law component.
        mu_g : ``float``
            Mean of the Gaussian peak (solar masses).
        sigma_g : ``float``
            Standard deviation of the Gaussian peak (solar masses).
        lambda_peak : ``float``
            Fraction of sources in the Gaussian component (0-1).
        delta_m : ``float``
            Low-mass smoothing width (solar masses).
        normalization_size : ``int``
            Grid size for CDF computation. 
            default: 500

        Returns
        -------
        samples : ``numpy.ndarray``
            Mass samples drawn from the specified distribution (solar masses).
    """

    cdf_values, m_grid = powerlaw_plus_peak_cdf(size=normalization_size, mminbh=mminbh, mmaxbh=mmaxbh, alpha=alpha, mu_g=mu_g, sigma_g=sigma_g, lambda_peak=lambda_peak, delta_m=delta_m)

    samples = inverse_transform_sampler(size=size, cdf=cdf_values, x=m_grid)

    return samples
