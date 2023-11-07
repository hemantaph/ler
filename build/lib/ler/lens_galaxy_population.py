# -*- coding: utf-8 -*-
"""
This module contains the LensGalaxyPopulation class, which is used to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed, image properties, and lensed SNRs. \n
The class inherits from the CompactBinaryPopulation class, which is used to sample source parameters. \n
"""
import warnings

warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import gengamma, rayleigh, norm
from scipy.interpolate import interp1d
from scipy.integrate import quad
from lenstronomy.Util.param_util import phi_q2_ellipticity

# for redshift to luminosity distance conversion
from astropy.cosmology import Planck18
import astropy.units as u
from astropy import constants as const

# the following .py file will be called if they are not given in the class initialization
from ler.source_population import CompactBinaryPopulation
from ler.helperroutines import add_dictionaries_together, trim_dictionary

# for multiprocessing
from multiprocessing import Pool
from tqdm import tqdm

# multiprocessing routines
import ler.multiprocessing_routine as mp


class LensGalaxyPopulation:
    """
    Class to sample lens galaxy parameters

    Parameters
    ----------
    CompactBinaryPopulation_ : CompactBinaryPopulation class
        This is an already initialized class that contains a function (CompactBinaryPopulation.sample_gw_parameters) that actually samples the source parameters. \n
        :class:`~ler.source_population.CompactBinaryPopulation`

    Examples
    --------
    >>> from ler.lens_galaxy_population import LensGalaxyPopulation
    >>> lens_pop = LensGalaxyPopulation()
    >>> # list all the methods of the class
    >>> print([method for method in dir(lens_pop) if method.startswith('__') is False])
    ['Dc_to_z', 'angular_diameter_distance', 'angular_diameter_distance_z1z2', 'cbc_pop', 'compute_einstein_radii', 'create_lookup_table', 'differential_comoving_volume', 'get_image_properties', 'get_lensed_snrs', 'lens_redshift_sampler_helper_function', 'm_max', 'm_min', 'normalization_pdf_z', 'rejection_sample_lensing_probability', 'sample_axis_ratio_angle_phi', 'sample_galaxy_shear', 'sample_gamma', 'sample_lens_parameters', 'sample_lens_parameters_routine', 'sample_lens_redshifts', 'sample_strongly_lensed_source_parameters', 'sample_velocity_dispersion_axis_ratio', 'strong_lensing_optical_depth', 'z_max', 'z_min', 'z_to_Dc', 'z_to_luminosity_distance']
    >>> # sample lens parameters
    >>> lens_parameters = lens_pop.sample_lens_parameters(size=1000)
    >>> lens_parameters.keys()
    dict_keys(['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a2', 'tilt1', 'tilt2', 'phi12', 'phi_jl'])
    >>> # get image properties
    >>> lens_parameters = lens_pop.get_image_properties(lens_parameters, n_min_images=2, n_max_images=4, lensModelList=['EPL_NUMBA', 'SHEAR'], npool=4)
    solving lens equations...
    100%|█████████████████████████████████████████████████████████| 1000/1000 [00:00<00:00, 1258.38it/s]
    >>> lens_parameters.keys()
    dict_keys(['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a2', 'tilt1', 'tilt2', 'phi12', 'phi_jl', 'n_images', 'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'image_type', 'weights'])
    >>> # get lensed SNRs
    >>> from gwsnr import GWSNR
    >>> snr = GWSNR()
    Given: IMR waveform
    psds not given. Choosing bilby's default psds
    given psds:  {'L1': 'aLIGO_O4_high_asd.txt', 'H1': 'aLIGO_O4_high_asd.txt', 'V1': 'AdV_asd.txt'}
    Interpolator will be generated for L1 detector at ./interpolator_pickle/L1/halfSNR_dict_0.pickle
    Interpolator will be generated for H1 detector at ./interpolator_pickle/H1/halfSNR_dict_0.pickle
    Interpolator will be generated for V1 detector at ./interpolator_pickle/V1/halfSNR_dict_0.pickle
    Generating interpolator for ['L1', 'H1', 'V1'] detectors
    interpolation for each mass_ratios: 100%|███████████████████████████| 50/50 [00:23<00:00,  2.10it/s]
    interpolator generated
    >>> lens_snrs = lens_pop.get_lensed_snrs(snr, lens_parameters, n_max_images=4)
    >>> lens_snrs.keys()
    dict_keys(['opt_snr_net', 'L1', 'H1', 'V1'])

    Instance Attributes
    ----------
    LensGalaxyPopulation class has the following instance attributes:\n
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`~cbc_pop`                     | CompactBinaryPopulation class    |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_min`                       | float                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_max`                       | float                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~m_min`                       | float                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~m_max`                       | float                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~normalization_pdf_z`         | float                            |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------
    LensGalaxyPopulation class has the following instance methods:\n
    +-------------------------------------+----------------------------------+
    | Methods                             | Type                             |
    +=====================================+==================================+
    |:meth:`~create_lookup_table`         | Function to create a lookup      |
    |                                     | table for the differential       |
    |                                     | comoving volume and luminosity   |
    |                                     | distance wrt redshift            |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_lens_parameters`      | Function to sample lens galaxy   |
    |                                     | parameters                       |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_lens_parameters_routine`                                 |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample lens galaxy   |
    |                                     | parameters                       |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_strongly_lensed_source_parameters`                       |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample source        |
    |                                     | parameters conditioned on the    |
    |                                     | source being strongly lensed     |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_lens_redshifts`       | Function to sample lens redshifts|
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_velocity_dispersion_axis_ratio`                          |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample velocity      |
    |                                     | dispersion and axis ratio of the |
    |                                     | lens galaxy                      |
    +-------------------------------------+----------------------------------+
    |:meth:`~compute_einstein_radii`      | Function to compute the Einstein |
    |                                     | radii of the lens galaxies       |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_axis_ratio_angle_phi` | Function to sample the axis      |
    |                                     | rotation angle of the elliptical |
    |                                     | lens galaxy                      |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_galaxy_shear`         | Function to sample the lens      |
    |                                     | galaxy shear                     |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_gamma`                | Function to sample the lens      |
    |                                     | galaxy spectral index of the     |
    |                                     | density profile                  |
    +-------------------------------------+----------------------------------+
    |:meth:`~rejection_sample_lensing_probability`                           |
    +-------------------------------------+----------------------------------+
    |                                     | Function to conduct rejection    |
    |                                     | sampling wrt einstein radius     |
    +-------------------------------------+----------------------------------+
    |:meth:`~strong_lensing_optical_depth`| Function to compute the strong   |
    |                                     | lensing optical depth            |
    +-------------------------------------+----------------------------------+
    |:meth:`~get_image_properties`        | Function to get the image        |
    |                                     | properties e.g. image positions, |
    |                                     | magnifications, time delays, etc.|
    +-------------------------------------+----------------------------------+
    |:meth:`~get_lensed_snrs`             | Function to get the lensed SNRs  |
    +-------------------------------------+----------------------------------+

    """

    # Attributes
    cbc_pop = None
    """:class:`~CompactBinaryPopulation` class\n
    This is an already initialized class that contains a function (CompactBinaryPopulation.sample_gw_parameters) that actually samples the source parameters. 
    """

    z_min = None
    """`float`\n
    minimum redshift
    """
    z_max = None
    """`float`\n
    maximum redshift
    """

    m_min = None
    """`float`\n
    minimum mass in detector frame
    """

    m_max = None
    """`float`\n
    maximum mass in detector frame
    """

    normalization_pdf_z = None
    """`float`\n
    normalization constant of the pdf p(z)
    """

    def __init__(self, CompactBinaryPopulation_=False):
        if CompactBinaryPopulation_ == False:
            # initialization of clasess
            # CompactBinaryPopulation already inherits from Source_Galaxy_Population_Model class form source_population.py
            self.cbc_pop = CompactBinaryPopulation(
                z_min=0.0001,
                z_max=10.0,
                m_min=4.59,
                m_max=86.22,
                event_type="popI_II",
                merger_rate_density_param=None,
                src_model_params=None,
            )
        else:
            # if the classes are already initialized, then just use them
            self.cbc_pop = CompactBinaryPopulation_

        self.z_min = self.cbc_pop.z_min
        self.z_max = self.cbc_pop.z_max
        self.m_min = self.cbc_pop.m_min
        self.m_max = self.cbc_pop.m_max
        self.create_lookup_table(z_min=self.z_min, z_max=self.z_max)

        # To find the normalization constant of the pdf p(z)
        # this under the assumption that the event is strongly lensed
        # Define the merger-rate density function
        merger_rate_density_detector_frame = lambda z: self.cbc_pop.merger_rate_density(
            z
        ) / (1 + z)
        # Define the pdf p(z)
        pdf_unnormalized = (
            lambda z: merger_rate_density_detector_frame(z)
            * self.differential_comoving_volume(z)
            * self.strong_lensing_optical_depth_SIS(z)
        )
        # Normalize the pdf
        # this normalization factor is common no matter what you choose for z_min and z_max
        self.normalization_pdf_z = quad(pdf_unnormalized, self.z_min, self.z_max)[0]
        return None

    def create_lookup_table(self, z_min, z_max):
        """
        Functions to create lookup tables
        1. Redshift to co-moving distance.
        2. Co-moving distance to redshift.
        3. Redshift to luminosity distance
        4. Redshift to angular diameter distance.
        5. Lens redshift sampler helper function.
        6. Redshift to differential comoving volume.

        Parameters
        ----------
            z_min : `float`
                minimum redshift
            z_max : `float`
                maximum redshift

        """
        # initialing cosmological functions for fast calculation through interpolation
        # z_min is fixed to 0.0, as the lens redshifts are drawn between 0.0 and z_max
        z = np.linspace(0.0, z_max, 500)  # red-shift
        Dc = Planck18.comoving_distance(z).value  # co-moving distance in Mpc
        self.z_to_Dc = interp1d(z, Dc, kind="cubic")
        self.Dc_to_z = interp1d(Dc, z, kind="cubic")
        self.z_to_luminosity_distance = self.cbc_pop.z_to_luminosity_distance

        # for angular diameter distance
        quad_ = []  # refers to integration (with quad algorithm) from scipy
        for ii in range(len(z)):
            quad_.append(
                quad(
                    Planck18._inv_efunc_scalar,
                    0.0,
                    z[ii],
                    args=Planck18._inv_efunc_scalar_args,
                )[0]
            )
        quad_ = np.array(quad_)
        quad_int = interp1d(z, np.array(quad_), kind="cubic")

        H0d = Planck18._hubble_distance.value
        self.angular_diameter_distance_z1z2 = (
            lambda zl0, zs0: H0d * (quad_int(zs0) - quad_int(zl0)) / (zs0 + 1.0)
        )
        self.angular_diameter_distance = lambda zs0: H0d * quad_int(zs0) / (zs0 + 1.0)

        # create a lookup table for the lens redshift draws
        r = np.linspace(0, 1, num=100)
        u = (
            10 * r**3 - 15 * r**4 + 6 * r**5
        )  # See the integral of Eq. A7 of https://arxiv.org/pdf/1807.07062.pdf (cdf)
        self.lens_redshift_sampler_helper_function = interp1d(
            u, r, kind="cubic"
        )  # Computes r(u)

        # Create a lookup table for the differential comoving volume
        self.differential_comoving_volume = self.cbc_pop.differential_comoving_volume
        return None

    def sample_lens_parameters(
        self, size=1000, lens_parameters_input={}, verbose=False
    ):
        """
        Function to sample galaxy lens parameters

        Parameters
        ----------
            size : `int`
                number of lens parameters to sample
            lens_parameters_input : `dict`
                dictionary of lens parameters to sample

        Returns
        -------
            lens_parameters : `dict`
                dictionary of lens parameters and source parameters (lens conditions applied)
                e.g. dictionary keys:\n
                lensing related=>['zl':redshift of lens, 'zs': redshift of source, 'sigma':velocity dispersion, 'q':axis ratios, 'e1':ellipticity, 'e2':ellipticity, 'gamma1':external-shear, 'gamma2':external-shear, 'Dl':angular diameter distance of lens, 'Ds':angular diameter distance of source, 'Dls':angular diameter distance between lens and source, 'theta_E': einstein radius in radian, 'gamma':spectral index of mass density distribution] \n
                source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'iota': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination, 'a_1': spin magnitude of the more massive black hole, 'a2': spin magnitude of the less massive black hole, 'tilt_1': tilt angle of the more massive black hole, 'tilt_2': tilt angle of the less massive black hole, 'phi_12': azimuthal angle between the two spins, 'phi_jl': azimuthal angle between the total angular momentum and the orbital angular momentum]

        """
        if verbose:
            print(
                f"refer to https://arxiv.org/abs/1908.06068 for parameter definitions, priors and sampling methods"
            )

        return self.sample_lens_parameters_routine(
            size=size, lens_parameters_input=lens_parameters_input
        )

    def sample_lens_parameters_routine(self, size=1000, lens_parameters_input={}):
        """
        Function to sample galaxy lens parameters

        Parameters
        ----------
            size : `int`
                number of lens parameters to sample
            lens_parameters_input : `dict`
                dictionary of lens parameters to sample

        Returns
        -------
            lens_parameters : `dict`
                dictionary of lens parameters and source parameters (lens conditions applied)
                e.g. dictionary keys:\n
                lensing related=>['zl':redshift of lens, 'zs': redshift of source, 'sigma':velocity dispersion, 'q':axis ratios, 'e1':ellipticity, 'e2':ellipticity, 'gamma1':external-shear, 'gamma2':external-shear, 'Dl':angular diameter distance of lens, 'Ds':angular diameter distance of source, 'Dls':angular diameter distance between lens and source, 'theta_E': einstein radius in radian, 'gamma':spectral index of mass density distribution] \n
                source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'iota': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination, 'a_1': spin magnitude of the more massive black hole, 'a2': spin magnitude of the less massive black hole, 'tilt_1': tilt angle of the more massive black hole, 'tilt_2': tilt angle of the less massive black hole, 'phi_12': azimuthal angle between the two spins, 'phi_jl': azimuthal angle between the total angular momentum and the orbital angular momentum]

        """

        # Sample source redshifts from the source population
        # rejection sampled with optical depth
        # print('sampling source parameters...')
        source_params_strongly_lensed = self.sample_strongly_lensed_source_parameters(
            size=size
        )
        zs = source_params_strongly_lensed["zs"]

        # Sample lens redshifts
        # print("sampling lens's redshifts...")
        zl = self.sample_lens_redshifts(zs)
        # Sample velocity dispersions and axis ratios (note: these should be sampled together because the lensing probability depends on the combination of these two parameters)
        sigma, q = self.sample_velocity_dispersion_axis_ratio(zs)

        # Compute the Einstein radii
        # print("sampling einstein radius...")
        theta_E = self.compute_einstein_radii(sigma, zl, zs)

        # Sample the axis ratio angle
        # print("sampling axis ratio and angle...")
        axis_ratio_angle_phi = self.sample_axis_ratio_angle_phi(size=size)

        # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
        e1, e2 = phi_q2_ellipticity(axis_ratio_angle_phi, q)

        # Sample shears
        # print("sampling external shears...")
        gamma1, gamma2 = self.sample_galaxy_shear(size=size)

        # Sample the spectral index of the mass density distribution
        # print("sampling spectral index...")
        gamma = self.sample_gamma(size=size)

        # Compute the angular diameter distances
        Dl = self.angular_diameter_distance(zl)  # for the lens
        Ds = self.angular_diameter_distance(zs)  # for the source
        # Compute Dls also
        Dls = self.angular_diameter_distance_z1z2(zl, zs)  # for the lens-source pair

        # Create a dictionary of the lens parameters
        lens_parameters = {
            "zl": zl,
            "zs": zs,
            "sigma": sigma,
            "q": q,
            "e1": e1,
            "e2": e2,
            "gamma1": gamma1,
            "gamma2": gamma2,
            "Dl": Dl,
            "Ds": Ds,
            "Dls": Dls,
            "theta_E": theta_E,
            "gamma": gamma,
        }

        # Add source params strongly lensed to the lens params
        lens_parameters.update(source_params_strongly_lensed)

        # Rejection sample based on the lensing probability, that is, rejection sample wrt theta_E
        mask = self.rejection_sample_lensing_probability(
            theta_E
        )  # proportional to pi theta_E^2
        lens_parameters = {key: val[mask] for key, val in lens_parameters.items()}

        # Add the lensing parameter dictionaries together
        # Note: This is done after rejection sampling, so that we don't have to worry about the size of the dictionary
        lens_parameters = add_dictionaries_together(
            lens_parameters, lens_parameters_input
        )

        # Check if the lens are larger than requested size
        # len_ = len(lens_parameters['zl'])
        # print(f'len(zs) = {len_}')
        if len(lens_parameters["zl"]) >= size:
            # Trim dicitionary to right size
            lens_parameters = trim_dictionary(lens_parameters, size)
            return lens_parameters
        else:
            # Run iteratively until we have the right number of lensing parmaeters
            # print("current sampled size", len(lens_parameters['zl']))
            return self.sample_lens_parameters_routine(
                size=size, lens_parameters_input=lens_parameters
            )

    def sample_strongly_lensed_source_parameters(self, size=1000):
        """
        Function to sample source redshifts and other parameters, conditioned on the source being strongly lensed.

        Parameters
        ----------
            size : `int`
                number of lens parameters to sample

        Returns
        -------
            gw_param_strongly_lensed : `dict`
                dictionary of source parameters. `zs` is sampled considering the merger rate density at source frame, comoving volume and strong lensing optical depth. \n
                e.g. gw_param_strongly_lensed.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a2', 'tilt1', 'tilt2', 'phi12', 'phi_jl']

        """
        z_max = self.z_max
        z_min = self.z_min

        def zs_strongly_lensed_function(zs_strongly_lensed):
            # get zs
            zs = self.cbc_pop.sample_source_redshifts(
                size=size, z_min=z_min, z_max=z_max
            )

            # put strong lensing condition with optical depth
            tau = self.strong_lensing_optical_depth_SIS(zs)
            tau_max = np.max(self.strong_lensing_optical_depth_SIS(z_max))
            r = np.random.uniform(0, tau_max, size=len(zs))
            pick_strongly_lensed = r < tau  # pick strongly lensed sources
            # Add the strongly lensed source redshifts to the list
            zs_strongly_lensed += list(zs[pick_strongly_lensed])  # list concatenation

            # Check if the zs_strongly_lensed are larger than requested size
            if len(zs_strongly_lensed) >= size:
                # Trim dicitionary to right size
                zs_strongly_lensed = zs_strongly_lensed[:size]
                return zs_strongly_lensed
            else:
                # Run iteratively until we have the right number of lensing parmaeters
                # print("current sampled size", len(lens_parameters['zl']))
                return zs_strongly_lensed_function(zs_strongly_lensed)

        zs_strongly_lensed = []
        zs_ = zs_strongly_lensed_function(zs_strongly_lensed)
        # gravitional waves source parameter sampling
        gw_param_strongly_lensed = self.cbc_pop.sample_gw_parameters(
            zs=np.array(zs_), nsamples=size, verbose=False
        )

        return gw_param_strongly_lensed

    def sample_lens_redshifts(self, zs):
        """
        Function to sample lens redshifts, conditioned on the lens being strongly lensed
        Input parameters:
            zs : source redshifts
        Output parameters:
            zl : lens redshifts
        """
        # lens redshift distribution
        r = self.lens_redshift_sampler_helper_function(
            np.random.uniform(0, 1, size=len(zs))
        )
        # comoing distance to the lens galaxy
        # on the condition that lens lie between the source and the observer
        lens_galaxy_Dc = (
            self.z_to_Dc(zs) * r
        )  # corresponding element-wise multiplication between 2 arrays
        # lens redshift
        zl = self.Dc_to_z(lens_galaxy_Dc)  # 2D array
        return zl

    def sample_velocity_dispersion_axis_ratio(self, zs):
        """
        Function to sample velocity dispersion and axis ratio of the lens galaxy

        Parameters
        ----------
            zs : `float`
                source redshifts

        Returns
        -------
            sigma : `float`
                velocity dispersion of the lens galaxy
            q : `float`
                axis ratio of the lens galaxy

        """
        size = len(zs)
        sigma = []
        q = []

        # Draw the velocity dispersions and axis ratios in 'chunks'
        size_ = size
        while size_ != 0:
            # Draw the velocity dispersion
            a = gengamma.rvs(2.32 / 2.67, 2.67, size=size_)
            sigma_ = 161.0 * a

            # Draw the axis ratio see Appendix of https://arxiv.org/pdf/1807.07062.pdf
            s = 0.38 - 0.09177 * a
            b = rayleigh.rvs(scale=s, size=size_)
            q_ = 1.0 - b

            # Weed out sigmas and axis ratios that have axis ratio below 0.2
            idx = q_ > 0.2
            sigma_ = sigma_[idx]
            q_ = q_[idx]

            # Append the velocity dispersions and axis ratios
            sigma += list(sigma_)
            q += list(q_)

            size_ = abs(size - len(sigma))

        # Transform to an array
        sigma = np.array(sigma)
        q = np.array(q)
        return sigma, q

    def compute_einstein_radii(self, sigma, zl, zs):
        """
        Function to compute the Einstein radii of the lens galaxies

        Parameters
        ----------
            sigma : `float`
                velocity dispersion of the lens galaxy
            zl : `float`
                lens redshifts
            zs : `float`
                source redshifts

        Returns
        -------
            theta_E : `float`
                Einstein radii of the lens galaxies in radian

        """
        # Compute the angular diameter distances
        Ds = self.angular_diameter_distance(zs)
        Dls = self.angular_diameter_distance_z1z2(zl, zs)
        # Compute the Einstein radii
        theta_E = (
            4.0 * np.pi * (sigma / const.c.to("km/s").value) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc
        return theta_E

    def sample_axis_ratio_angle_phi(self, size=1000):
        """
        Function to sample the axis rotation angle of the elliptical lens galaxy

        Parameters
        ----------
            size : `int`
                number of lens parameters to sample

        Returns
        -------
            phi : `float`
                axis rotation angle of the elliptical lens galaxy

        """
        # Draw the angles
        phi = np.random.uniform(0, 2 * np.pi, size=size)
        return phi

    def sample_galaxy_shear(self, size):
        """
        Function to sample the lens galaxy shear

        Parameters
        ----------
            size : `int`
                number of lens parameters to sample

        Returns
        -------
            gamma_1 : `float`
                shear component in the x-direction
            gamma_2 : `float`
                shear component in the y-direction

        """
        # Draw an external shear
        gamma_1 = norm.rvs(size=size, scale=0.05)
        gamma_2 = norm.rvs(size=size, scale=0.05)
        return gamma_1, gamma_2

    def sample_gamma(self, size=1000):
        """
        Function to sample the lens galaxy spectral index of the density profile

        Parameters
        ----------
            size : `int`
                number of lens parameters to sample

        Returns
        -------
            gamma : `float`
                spectral index of the density profile

        """
        self.gamma_mean = 2
        self.gamma_std = 0.2
        return np.random.normal(loc=self.gamma_mean, scale=self.gamma_std, size=size)

    def rejection_sample_lensing_probability(self, theta_E):
        """
        Function to conduct rejection sampling wrt einstein radius

        Parameters
        ----------
            theta_E : `float`
                Einstein radii of the lens galaxies

        Returns
        -------
            idx : `bool`
                boolean array of size len(theta_E) indicating whether the sample is accepted or not

        """
        size = len(theta_E)
        theta_E_max = np.max(theta_E)
        u = np.random.uniform(0, theta_E_max**2, size=size)
        idx = u < theta_E**2
        return idx

    def strong_lensing_optical_depth_SIE(self, zs):
        """
        Function to compute the strong lensing optical depth SIE

        Parameters
        ----------
            zs : `float`
                source redshifts

        Returns
        -------
            tau : `float`
                strong lensing optical depth

        """
        # for SIE model
        # dedine a function to calculate the cross section number
        # here we already assume that we have initialized the crossection spline 
        def getcrosssect_num(theta_E, q):
            fid_b_I = 10.0  # constant
            idx = q > 0.999
            q[idx] = 0.999
            idx = q < 0.1
            q[idx] = 0.1
            b_I = theta_E * np.sqrt(q)
            return self.cross_sect_spl(q)*(b_I/fid_b_I)**2

        # theta_E=bsis=einstein_radius_SIS  # refer to compute_einstein_radii
        
        return getcrosssect_num(theta_E, q)/(4*np.pi)

    def strong_lensing_optical_depth_SIS(self, zs):
        """
        Function to compute the strong lensing optical depth (SIS)

        Parameters
        ----------
            zs : `float`
                source redshifts

        Returns
        -------
            tau : `float`
                strong lensing optical depth

        """
        # For SIS model
        # z to luminosity_distance (luminosity_distance) conversion
        Dc = self.z_to_Dc(zs) * 1e-3  # 1e-3 converts Mpc to Gpc

        return (Dc / 62.2) ** 3

    def get_image_properties(
        self,
        lens_parameters,
        n_min_images=int(2),
        n_max_images=int(4),
        lensModelList=["EPL_NUMBA", "SHEAR"],
        npool=4,
    ):
        """
        Function to get the image properties e.g. image positions, magnifications, time delays, etc.

        Parameters
        ----------
            lens_parameters : `dict`
                dictionary of lens parameters
                e.g. lens_parameters.keys() = ['zs', 'zl', 'gamma1', 'gamma2', 'e1', 'e2', 'gamma', 'theta_E']
            n_min_images : `int`
                minimum number of images to consider
                default: 2
            n_max_images : `int`
                maximum number of images to consider
                default: 4
            lensModelList : `list`
                list of lens models
                default: ['EPL_NUMBA', 'SHEAR']
            npool : `int`
                number of processes to use
                default: 4

        Returns
        -------
            lens_parameters : `dict`
                dictionary of lens parameters and image properties
                e.g. lens_parameters contains the following keys:\n
                lens related=>['zs': source redshift, 'zl': lens redshift, 'gamma1': shear component in the x-direction, 'gamma2': shear component in the y-direction, 'e1': ellipticity component in the x-direction, 'e2': ellipticity component in the y-direction, 'gamma': spectral index of the mass density distribution, 'theta_E': einstein radius in radian]\n
                source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'iota': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination, 'a_1': spin magnitude of the more massive black hole, 'a2': spin magnitude of the less massive black hole, 'tilt_1': tilt angle of the more massive black hole, 'tilt_2': tilt angle of the less massive black hole, 'phi_12': azimuthal angle between the two spins, 'phi_jl': azimuthal angle between the total angular momentum and the orbital angular momentum]\n
                image related=>['x_source': source position in the x-direction, 'y_source': source position in the y-direction, 'x0_image_position': image position in the x-direction, 'x1_image_position': image position in the y-direction, 'magnifications': magnifications, 'time_delays': time delays, 'n_images': number of images formed, 'determinant': determinants, 'trace': traces, 'iteration': to keep track of the iteration number, 'weights': weights for the caustic considered]

        """
        zs = lens_parameters["zs"]
        size = len(zs)
        zl = lens_parameters["zl"]
        # external shear params to the 'PEMD' galaxy lens
        gamma1, gamma2 = lens_parameters["gamma1"], lens_parameters["gamma2"]
        # ellipticity of the galaxy lens
        e1, e2 = lens_parameters["e1"], lens_parameters["e2"]
        gamma = lens_parameters["gamma"]
        einstein_radius = lens_parameters["theta_E"]
        # Create the lens model list (note: can be a different lens model for different samples)
        lensModelList = np.array(lensModelList) * np.ones(
            (size, len(lensModelList)), dtype=object
        )
        min_img_arr = n_min_images * np.ones((size), dtype=int)

        # get image properties (with Multiprocessing)
        iterations = np.arange(size)
        input_arguments = np.array(
            [
                min_img_arr,
                e1,
                e2,
                gamma,
                gamma1,
                gamma2,
                zl,
                zs,
                einstein_radius,
                iterations,
            ],
            dtype=object,
        ).T
        input_arguments = np.concatenate((input_arguments, lensModelList), axis=1)
        # Initialize the image positions and lens argument list here.
        x0_image_positions = np.ones((size, n_max_images)) * np.nan
        x1_image_positions = np.ones((size, n_max_images)) * np.nan
        magnifications = np.ones((size, n_max_images)) * np.nan
        time_delays = np.ones((size, n_max_images)) * np.nan
        determinants = np.ones((size, n_max_images)) * np.nan
        traces = np.ones((size, n_max_images)) * np.nan
        n_images = np.ones(size, dtype=int)
        x_source, y_source = np.ones(size) * np.nan, np.ones(size) * np.nan
        weights = np.ones(size) * np.nan

        # Solve the lens equation
        print("solving lens equations...")
        if n_min_images == 2:
            solve_lens_equation = mp.solve_lens_equation1
        elif n_min_images > 2:
            solve_lens_equation = mp.solve_lens_equation2
        else:
            raise ValueError("n_min_images should be greater than 1")
        with Pool(processes=npool) as pool:
            # call the same function with different data in parallel
            # imap->retain order in the list, while map->doesn't
            for result in tqdm(
                pool.imap(solve_lens_equation, input_arguments),
                total=len(input_arguments),
                ncols=100,
                disable=False,
            ):
                # print(result)
                """
                for i in tqdm(range(size)):
                    result = self.solve_lens_equation(input_arguments[i])
                """
                (
                    x_source_i,
                    y_source_i,
                    x0_image_position_i,
                    x1_image_position_i,
                    magnifications_i,
                    time_delays_i,
                    n_image_i,
                    determinant_i,
                    trace_i,
                    iter_i,
                    weights_i,
                ) = result

                n_image_i = min(n_image_i, n_max_images)
                n_images[iter_i] = n_image_i
                x0_image_position = np.ones(n_max_images) * np.nan
                x1_image_position = np.ones(n_max_images) * np.nan
                x0_image_position[:n_image_i] = x0_image_position_i[:n_image_i]
                x1_image_position[:n_image_i] = x1_image_position_i[:n_image_i]
                x0_image_positions[
                    iter_i
                ] = x0_image_position  # shape = (size, n_max_images)
                x1_image_positions[
                    iter_i
                ] = x1_image_position  # shape = (size, n_max_images)
                magnification = np.ones(n_max_images) * np.nan
                time_delay = np.ones(n_max_images) * np.nan
                determinant = np.ones(n_max_images) * np.nan
                trace = np.ones(n_max_images) * np.nan
                magnification[:n_image_i] = magnifications_i[:n_image_i]
                time_delay[:n_image_i] = time_delays_i[:n_image_i]
                determinant[:n_image_i] = determinant_i[:n_image_i]
                trace[:n_image_i] = trace_i[:n_image_i]
                # Add the magnifications, time delays, determinants, and traces to their respective arrays
                magnifications[iter_i] = magnification
                time_delays[iter_i] = time_delay
                determinants[iter_i] = determinant
                traces[iter_i] = trace
                x_source[iter_i] = x_source_i
                y_source[iter_i] = y_source_i
                weights[iter_i] = weights_i

        # time-delays: convert to positive values
        # time-delays will be relative to the first arrived signal of an lensed event
        time_delays = time_delays - np.array([np.sort(time_delays, axis=1)[:, 0]]).T

        # select only strongly lensed events are selected
        assert np.all(n_images >= 2), "There are events with no images!"

        # image type classification (morse phase)
        number_of_lensed_events = size
        image_type = np.zeros((number_of_lensed_events, n_max_images))
        image_type[traces < 0] = 3
        image_type[traces > 0] = 1
        image_type[determinants < 0] = 2

        # Return a dictionary with all of the lens information but also the BBH parameters from gw_param
        image_parameters = {
            "n_images": n_images,
            "x0_image_positions": x0_image_positions,
            "x1_image_positions": x1_image_positions,
            "magnifications": magnifications,
            "time_delays": time_delays,
            "image_type": image_type,
            "weights": weights,
        }
        lens_parameters.update(image_parameters)

        return lens_parameters

    def get_lensed_snrs(self, snr_calculator, lensed_param, n_max_images=4):
        """
        Function to calculate the signal to noise ratio for each image in each event.

        Parameters
        ----------
            snr_calculator : `class`
                snr_calculator class
                this is an already initialized class that contains a function (snr_calculator.snr) that actually calculates snr with the given gw_params.\n
                Luminosity distance and time delay are modified to be effective luminosity distance and effective time delay, respectively, for each image using the magnifications and time delays.\n
            lensed_param : `dict`
                dictionary containing the both already lensed source paramters and image parameters.
                e.g. lensed_param.keys() = ['mass_1', 'mass_2', 'zs', 'luminosity_distance', 'iota', 'psi', 'phi', 'ra', 'dec', 'geocent_time', 'phase', 'a_1', 'a2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'magnifications', 'time_delays']
            n_max_images : `int`
                maximum number of images to consider
                default: 4

        Returns
        -------
            snrs : `dict`
                signal to noise ratio for each image in each event.
                (dictionary containing 'H1', 'L1', ..., and 'opt_snr_net', which is the network snr, for each image as an array with dimensions (number_of_lensed_events,n_max_images) )

        """
        # needed to calculate effective luminosity distance and effective time delay
        magnifications = lensed_param["magnifications"]
        time_delays = lensed_param["time_delays"]

        # Get the binary parameters
        number_of_lensed_events = len(magnifications)
        mass_1, mass_2, luminosity_distance, iota, psi, ra, dec, geocent_time, phase, a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl  = (
            lensed_param["mass_1"],
            lensed_param["mass_2"],
            lensed_param["luminosity_distance"],
            lensed_param["iota"],
            lensed_param["psi"],
            lensed_param["ra"],
            lensed_param["dec"],
            lensed_param["geocent_time"],
            lensed_param["phase"],
            lensed_param["a_1"],
            lensed_param["a_2"],
            lensed_param["tilt_1"],
            lensed_param["tilt_2"],
            lensed_param["phi_12"],
            lensed_param["phi_jl"],
        )

        # setting up snr dictionary
        detectors = snr_calculator.list_of_detectors
        optimal_snrs = dict()
        optimal_snrs["opt_snr_net"] = (
            np.ones((number_of_lensed_events, n_max_images)) * np.nan
        )
        for detector in detectors:
            optimal_snrs[detector] = (
                np.ones((number_of_lensed_events, n_max_images)) * np.nan
            )

        # LALSimulation cannot handle NaN
        if snr_calculator.snr_type == "inner_product":
            print("There will be {} progress bar iteration".format(n_max_images))

        for i in range(n_max_images):
            # Get the optimal signal to noise ratios for each image
            buffer = magnifications[:, i]
            idx = ~np.isnan(buffer)  # index of not-nan
            effective_luminosity_distance = luminosity_distance[idx] / np.sqrt(
                np.abs(buffer[idx])
            )
            effective_geocent_time = geocent_time[idx] + time_delays[idx, i]
            # if GPS time is negative, shift it
            # by a year until it is positive
            effective_geocent_time[
                effective_geocent_time < 0
            ] += 31556952  # number of seconds in a year

            # Each image has their own effective luminosity distance and effective geocent time
            if len(effective_luminosity_distance) != 0:
                # Returns a dictionary
                optimal_snr = snr_calculator.snr(
                    mass_1[idx],
                    mass_2[idx],
                    effective_luminosity_distance,
                    iota[idx],
                    psi[idx],
                    phase[idx],
                    effective_geocent_time,
                    ra[idx],
                    dec[idx],
                    a_1[idx],
                    a_2[idx],
                    tilt_1[idx],
                    tilt_2[idx],
                    phi_12[idx],
                    phi_jl[idx],
                    jsonFile=False,
                )

                optimal_snrs["opt_snr_net"][idx, i] = optimal_snr["opt_snr_net"]
                for detector in detectors:
                    optimal_snrs[detector][idx, i] = optimal_snr[detector]

        self.optimal_snrs = optimal_snrs
        return optimal_snrs
