from numba import njit
from multiprocessing import Pool
import numpy as np
from scipy.integrate import quad
from scipy.stats import gengamma
from astropy.cosmology import LambdaCDM
from ..utils import  interpolator_from_pickle, cubic_spline_interpolator, inverse_transform_sampler
from .jit_functions import phi_cut_SIE, axis_ratio_rayleigh, axis_ratio_SIS, phi, phi_loc_bernardi

class OpticalDepth():
    """
    Class to calculate the optical depth, velocity dispersion and axis-ratio of a lens galaxy population.

    Parameters
    ----------
    npool : `int`
        number of processors to use for multiprocessing
    z_min : `float`
        minimum redshift of the lens galaxy population
    z_max : `float`
        maximum redshift of the lens galaxy population
    optical_depth_function : `str` or `callable`
        Function or function name to calculate optical depth.
        Check for default/available optical depth functions by running,
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> print(OpticalDepth().available_optical_depth_list_and_its_params)
    sampler_priors, sampler_priors_params : `dict`, `dict`
        dictionary of sampler functions and it's parameters to sample velocity dispersion and axis-ratio.
        Check for default/available sampler priors and corresponding input parameters by running,
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> print(OpticalDepth().available_velocity_dispersion_list_and_its_params)
        >>> print(OpticalDepth().available_axis_ratio_list_and_its_params)
    cosmology : `astropy.cosmology`
        Cosmology to use
        default: None/astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)
    directory : `str`
        directory to store interpolator pickle files
        default: "./interpolator_pickle"
    create_new_interpolator : `dict`
        dictionary to create new interpolator for velocity dispersion and optical depth.

    Examples
    --------
    >>> from ler.lens_galaxy_population import OpticalDepth
    >>> od = OpticalDepth()
    >>> print(od.strong_lensing_optical_depth(0.5))

    Instance Attributes
    -------------------
    OpticalDepth class has the following attributes:
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`~npool`                       | `int`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_min`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_max`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~functions`                   | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~sampler_priors`              | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~sampler_priors_params`       | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~cosmo`                       | `astropy.cosmology`              |
    +-------------------------------------+----------------------------------+
    |:attr:`~directory`                   | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~c_n_i`                       | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~vd_inv_cdf`                  | `numpy.ndarray`                  |
    +-------------------------------------+----------------------------------+
    |:attr:`~available_velocity_dispersion_list_and_its_params`              |
    +-------------------------------------+----------------------------------+
    |                                     | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~available_axis_ratio_list_and_its_params`                       |
    +-------------------------------------+----------------------------------+
    |                                     | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~available_optical_depth_list_and_its_params`                    |
    +-------------------------------------+----------------------------------+
    |                                     | `dict`                           |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------------
    OpticalDepth class has the following methods:
    +-------------------------------------+----------------------------------+
    | Methods                             | Type                             |
    +=====================================+==================================+
    |:meth:`~create_lookup_table`         | Function to create a lookup      |
    |                                     | table for redshift dependent     |
    |                                     | cosmological functions           |
    +-------------------------------------+----------------------------------+
    |:meth:`~initialize_velocity_dispersion_sampler`                         |
    +-------------------------------------+----------------------------------+
    |                                     | Function to initialize velocity  |
    |                                     | dispersion sampler               |
    +-------------------------------------+----------------------------------+
    |:meth:`~initialize_optical_depth_function`                              |
    +-------------------------------------+----------------------------------+
    |                                     | Function to initialize optical   |
    |                                     | depth function                   |
    +-------------------------------------+----------------------------------+
    |:meth:`~velocity_dispersion_gengamma`| Function to sample velocity      |
    |                                     | dispersion from gengamma         |
    |                                     | distribution                     |
    +-------------------------------------+----------------------------------+
    |:meth:`~velocity_dispersion_bernardi`| Function to sample velocity      |
    |                                     | dispersion from Bernardi et al.  |
    |                                     | (2010)                           |
    +-------------------------------------+----------------------------------+
    |:meth:`~velocity_dispersion_ewoud`   | Function to sample velocity      |
    |                                     | dispersion from Wempe et al.     |
    |                                     | (2022)                           |
    +-------------------------------------+----------------------------------+
    |:meth:`~cross_section_SIS`           | Function to compute the SIS      |
    |                                     | cross-section                    |
    +-------------------------------------+----------------------------------+
    |:meth:`~tau_zl_zs`               | Function to compute the optical  |
    |                                     | depth integrand                  |
    +-------------------------------------+----------------------------------+
    |:meth:`~optical_depth_calculator`    | Function to compute the optical  |
    |                                     | depth without multiprocessing    |
    +-------------------------------------+----------------------------------+
    |:meth:`~optical_depth_multiprocessing`                                  |
    +-------------------------------------+----------------------------------+
    |                                     | Function to compute the optical  |
    |                                     | depth with multiprocessing       |
    +-------------------------------------+----------------------------------+
    |:meth:`~optical_depth_SIS_haris`     | Function to compute the strong   |
    |                                     | lensing optical depth (SIS)      |
    |                                     | Haris et al. (2018) (analytic)   |
    +-------------------------------------+----------------------------------+
    |:meth:`~optical_depth_SIS_hemanta`   | Function to compute the strong   |
    |                                     | lensing optical depth (SIS)      |
    |                                     | (numerical)                      |
    +-------------------------------------+----------------------------------+
    |:meth:`~optical_depth_SIE_hemanta`   | Function to compute the strong   |
    |                                     | lensing optical depth (SIE)      |
    |                                     | (numerical)                      |
    +-------------------------------------+----------------------------------+
    |:meth:`~strong_lensing_optical_depth`| Function to compute the strong   |
    |                                     | lensing optical depth            |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_velocity_dispersion`  | Function to sample velocity      |
    |                                     | dispersion                       |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_axis_ratio`           | Function to sample axis ratio    |
    +-------------------------------------+----------------------------------+
    """

    def __init__(self,
        npool=4,
        z_min=0.001,
        z_max=10.0,
        optical_depth_function=None,
        sampler_priors=None,
        sampler_priors_params=None,
        cosmology=None,
        directory="./interpolator_pickle",
        create_new_interpolator=False,
    ):

        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

        self.sampler_priors = dict(
            velocity_dispersion="velocity_dispersion_bernardi",
            axis_ratio="axis_ratio_rayleigh",
        )
        if sampler_priors:
            self.sampler_priors.update(sampler_priors)
        self.sampler_priors_params = dict(
            velocity_dispersion=dict(vd_min=0., vd_max=600.),
            axis_ratio=dict(q_min=0.2, q_max=1.),
        )
        if sampler_priors_params:
            self.sampler_priors_params.update(sampler_priors_params)

        optical_depth_function = optical_depth_function if optical_depth_function else "optical_depth_SIE_hemanta"

        self.directory = directory
        self.c_n_i = dict(
            velocity_dispersion=dict(create_new=False, resolution=500),
            axis_ratio=dict(create_new=False, resolution=500),
            optical_depth=dict(create_new=False, resolution=100),
            z_to_Dc=dict(create_new=False, resolution=500),
            Dc_to_z=dict(create_new=False, resolution=500),
            angular_diameter_distance=dict(create_new=False, resolution=500),
            differential_comoving_volume=dict(create_new=False, resolution=500),
        )
        if isinstance(create_new_interpolator, dict):
            self.c_n_i.update(create_new_interpolator)
        elif create_new_interpolator is True:
            self.c_n_i = dict(
                velocity_dispersion=dict(create_new=True, resolution=500),
                axis_ratio=dict(create_new=True, resolution=500),
                optical_depth=dict(create_new=True, resolution=100),
                z_to_Dc=dict(create_new=True, resolution=500),
                Dc_to_z=dict(create_new=True, resolution=500),
                angular_diameter_distance=dict(create_new=True, resolution=500),
                differential_comoving_volume=dict(create_new=True, resolution=500),
            )

        vd_name = self.sampler_priors["velocity_dispersion"]  # velocity dispersion sampler name
        tau_name = optical_depth_function  # optical depth function name

        self.create_lookup_table_fuction(self.z_max)

        # setting velocity dispersion sampler method
        self.initialize_velocity_dispersion_sampler(vd_name);

        # setting optical depth method
        if callable(optical_depth_function):
            self.strong_lensing_optical_depth = optical_depth_function
        else:
            self.initialize_optical_depth_function(tau_name, vd_name);
        
    def initialize_velocity_dispersion_sampler(self, vd_name):
        """
        Function to initialize velocity dispersion sampler

        Parameters
        ----------
        vd_name : `str`
            name of velocity dispersion sampler
        """

        # setting up input parameters for interpolation
        # check vd_min and vd_max exists in sampler_priors_params
        try:
            self.vd_min = self.sampler_priors_params['velocity_dispersion']['vd_min']
            self.vd_max = self.sampler_priors_params['velocity_dispersion']['vd_max']
        except:
            self.vd_min = 0.
            self.vd_max = 600.

        # generating inverse cdf interpolator for velocity dispersion
        param_dict_given_ = dict(z_min=self.z_min, z_max=self.z_max, vd_min=self.vd_min, vd_max=self.vd_max, cosmology=self.cosmo, name=vd_name, resolution=self.c_n_i["velocity_dispersion"]["resolution"])
        sub_directory_ = vd_name
        x_ = np.linspace(self.vd_min, self.vd_max, self.c_n_i["velocity_dispersion"]["resolution"])
        category_ = "inv_cdf"
        create_new_ = self.c_n_i["velocity_dispersion"]["create_new"]
        
        if vd_name == "velocity_dispersion_gengamma":
            # setting up input parameters for interpolation
            pdf_func_ = lambda vd_: gengamma.pdf(vd_/161., a=2.32 / 2.67, c=2.67)   # gengamma pdf
            conditioned_y_ = None
            dimension_ = 1
     
        elif vd_name == "velocity_dispersion_bernardi":
            # setting up input parameters for interpolation
            pdf_func_ = lambda vd_: vd_**4*phi_loc_bernardi(sigma=vd_, cosmology_h=self.cosmo.h)
            conditioned_y_ = None
            dimension_ = 1

        elif vd_name == "velocity_dispersion_ewoud":
            if self.z_min==0.0:
                self.zl_list = np.linspace(self.z_min+0.001, self.z_max, 100)
            else:
                self.zl_list = np.linspace(self.z_min, self.z_max, 100)
            # setting up input parameters for interpolation
            pdf_func_ = lambda vd_, zl_: phi(vd_, zl_, cosmology_h=self.cosmo.h)*self.differential_comoving_volume(np.array([zl_]))
            conditioned_y_ = self.zl_list
            dimension_ = 2

        else:
            # this is to allow user to pass in their own sampler
            if callable(self.sampler_priors["velocity_dispersion"]):
                pass
            else:
                # raise error and terminate
                raise ValueError("velocity dispersion name not in `ler` library. And sampler_priors['velocity_dispersion'] is not a callable function.")
            
        self.vd_inv_cdf = interpolator_from_pickle(
            param_dict_given = param_dict_given_,
            directory=self.directory,
            sub_directory=sub_directory_,
            name=vd_name,
            x = x_,
            pdf_func= pdf_func_,
            conditioned_y=conditioned_y_,
            dimension=dimension_,
            category=category_,
            create_new=create_new_,
        )

        self.sample_velocity_dispersion = self.sampler_priors["velocity_dispersion"]

    def initialize_optical_depth_function(self, tau_name, vd_name):
        """
        Function to initialize optical depth function.

        Parameters
        ----------
        tau_name : `str`
            name of optical depth function
        vd_name : `str`
            name of velocity dispersion sampler
        """

        if tau_name=="optical_depth_SIS_haris":
            self.sample_axis_ratio = axis_ratio_SIS
            # no interpolation needed
            optical_depth_setter = getattr(self, tau_name)

        elif callable(tau_name):
            self.sample_axis_ratio = axis_ratio_SIS
            optical_depth_setter = tau_name

        else:
            # setting up input parameters for interpolation
            resolution = self.c_n_i["optical_depth"]["resolution"]
            param_dict_given_ = dict(z_min=self.z_min, z_max=self.z_max, vd_min=self.vd_min, vd_max=self.vd_max, cosmology=self.cosmo, tau_name=tau_name, vd_name=vd_name, q_name=self.sampler_priors["axis_ratio"], resolution=resolution)
            #self.param_dict_given_ = param_dict_given_
            sub_directory_ = tau_name
            if self.z_min==0.0:
                x_ = np.geomspace(self.z_min+0.001, self.z_max, resolution)
            else:
                x_ = np.geomspace(self.z_min, self.z_max, resolution)
            pdf_func_ = self.optical_depth_multiprocessing
            dimension_ = 1
            category_ = "function"
            create_new_ = self.c_n_i["optical_depth"]["create_new"]

            if tau_name == "optical_depth_SIS_hemanta" or tau_name == "SIS":
                from .mp import optical_depth_sis_mp
                self.sample_axis_ratio = axis_ratio_SIS
                self.tau_mp_routine = optical_depth_sis_mp

            elif tau_name=="optical_depth_SIE_hemanta" or tau_name == "SIE":
                # axis-ratio sampler
                if self.sampler_priors["axis_ratio"]=="axis_ratio_rayleigh":
                    self.sample_axis_ratio = axis_ratio_rayleigh
                else:
                    self.sample_axis_ratio = self.sampler_priors["axis_ratio"]

                if vd_name == "velocity_dispersion_ewoud":
                    from .mp import optical_depth_sie2_mp
                    self.tau_mp_routine = optical_depth_sie2_mp
                else:
                    from .mp import optical_depth_sie1_mp
                    self.tau_mp_routine = optical_depth_sie1_mp

            # this will initialize the interpolator
            optical_depth_setter = interpolator_from_pickle(
                param_dict_given = param_dict_given_,
                directory=self.directory,
                sub_directory=sub_directory_,
                name=tau_name,
                x = x_,
                pdf_func= pdf_func_,
                conditioned_y=None,
                dimension=dimension_,
                category=category_,
                create_new=create_new_,
            )

        #self.optical_depth_setter = optical_depth_setter
        self.strong_lensing_optical_depth = optical_depth_setter
        self.sample_axis_ratio = self.sampler_priors["axis_ratio"]

    def axis_ratio_rayleigh(self, sigma, q_min=0.2, q_max=1.0, get_attribute=False, param=None, **kwargs):
        """
        Function to sample axis ratio from rayleigh distribution with given velocity dispersion.

        Parameters
        ----------
        sigma : `float: array`
            velocity dispersion of the lens galaxy
        q_min, q_max : `float`
            minimum and maximum axis ratio
        get_attribute : `bool`
            if True, returns a function that can be used to sample axis ratio

        Returns
        -------
        q : `float: array`
            axis ratio of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(sampler_priors=dict(axis_ratio="axis_ratio_rayleigh"))
        >>> print(od.sample_axis_ratio(sigma=200.))
        """

        if param:
            q_min = param["q_min"]
            q_max = param["q_max"]
            
        if get_attribute:
            return njit(lambda sigma: axis_ratio_rayleigh(sigma, q_min, q_max))
        else:
            return axis_ratio_rayleigh(sigma, q_min, q_max)
        
    def axis_ratio_padilla_strauss(self, size=1000, q_min=0.2, q_max=1.0, get_attribute=False, param=None, **kwargs):
        """
        Function to sample axis ratio using Padilla and Strauss 2008 distribution for axis ratio

        Parameters
        ----------
        size : `int`
            sample size
        q_min, q_max : `float`
            minimum and maximum axis ratio
        get_attribute : `bool`
            if True, returns a function that can be used to sample axis ratio

        Returns
        -------
        q : `float: array`
            axis ratio of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(sampler_priors=dict(axis_ratio="axis_ratio_padilla_strauss"))
        >>> print(od.sample_axis_ratio(size=10))
        """

        try:
            size = len(kwargs["sigma"])
        except:
            pass
        if param:
            q_min = param["q_min"]
            q_max = param["q_max"]

        # Using Padilla and Strauss 2008 distribution for axis ratio
        q = np.array([0.04903276402927845, 0.09210526315789469, 0.13596491228070173, 0.20789473684210524, 0.2899703729522482, 0.3230132450331126, 0.35350877192982455, 0.37946148483792264, 0.4219298245614036, 0.4689525967235971, 0.5075026141512723, 0.5226472638550018, 0.5640350877192983, 0.6096491228070177, 0.6500000000000001, 0.6864848379226213, 0.7377192982456142, 0.7787295224817011, 0.8007581038689441, 0.822786685256187, 0.8668438480306729, 0.8973684210526317, 0.9254385964912283])
        pdf = np.array([0.04185262687135349, 0.06114520695141845, 0.096997499638376, 0.1932510900336828, 0.39547914337673706, 0.49569751276216234, 0.6154609137685201, 0.7182049959882812, 0.920153741243567, 1.1573982157399754, 1.3353263628106684, 1.413149656448315, 1.5790713532948977, 1.7280185150744938, 1.8132994441344819, 1.8365803753840484, 1.8178662203211204, 1.748929843583365, 1.688182592496342, 1.6274353414093188, 1.4948487090314488, 1.402785526832393, 1.321844068356993])

        spline_coeff = interpolator_from_pickle(
            param_dict_given = dict(parameter="axis_ratio", pdf="Padilla and Strauss 2008 distribution for axis ratio"),
            directory=self.directory,
            sub_directory="axis_ratio",
            name="axis_ratio_spline_coeff",
            x = q,
            pdf_func=None,
            y=pdf,
            conditioned_y=None,
            dimension=1,
            category="function",
            create_new=True,
        )

        resolution = self.c_n_i["axis_ratio"]["resolution"]
        create_new = self.c_n_i["axis_ratio"]["create_new"]
        q_array = np.linspace(q_min, q_max, resolution)
        pdf_array = cubic_spline_interpolator(q_array, spline_coeff[0], spline_coeff[1])

        q_inv_cdf = interpolator_from_pickle(
            param_dict_given = dict(parameter="axis_ratio", pdf="Padilla and Strauss 2008 distribution for axis ratio", q_min=q_min, q_max=q_max, resolution=resolution),
            directory=self.directory,
            sub_directory="axis_ratio",
            name="axis_ratio",
            x = q_array,
            pdf_func=None,
            y=pdf_array,
            conditioned_y=None,
            dimension=1,
            category="inv_cdf",
            create_new=create_new,
        )

        # def sampler(sigma):
        #     size = len(sigma)
        #     return inverse_transform_sampler(size, q_inv_cdf[0], q_inv_cdf[1])
            
        if get_attribute:
            return njit(lambda sigma: inverse_transform_sampler(len(sigma), q_inv_cdf[0], q_inv_cdf[1]))
        else:
            return inverse_transform_sampler(size, q_inv_cdf[0], q_inv_cdf[1])

    def velocity_dispersion_gengamma(self, size, a=2.32 / 2.67, c=2.67, get_attribute=False, param=None, **kwargs):
        """
        Function to sample velocity dispersion from gengamma distribution

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        a,c : `float`
            parameters of gengamma distribution
            refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gengamma.html
        get_attribute : `bool`
            if True, returns a function that can be used to sample velocity dispersion
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(a=2.32 / 2.67, c=2.67)

        Returns
        -------
        sigma : `numpy.ndarray` (1D array of floats)
            velocity dispersion of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(sampler_priors=dict(velocity_dispersion="velocity_dispersion_gengamma"), sampler_priors_params=dict(velocity_dispersion=dict(a=2.32 / 2.67, c=2.67)))
        >>> print(od.sample_velocity_dispersion(size=10))

        """

        if param:
            a = param["a"]
            c = param["c"]

        if get_attribute:
            return lambda size: 161.*gengamma.rvs(a, c, size=size)
        else:
            # sample velocity dispersion from gengamma distribution
            return 161.*gengamma.rvs(a, c, size=size)  # km/s
    
    def velocity_dispersion_bernardi(self, size, get_attribute=False, **kwargs):
        """
        Function to sample velocity dispersion from Bernardi et al. (2010). This uses inverse transform sampling.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        get_attribute : `bool`
            if True, returns a function that can be used to sample velocity dispersion

        Returns
        -------
        sigma : `numpy.ndarray` (1D array of floats)
            velocity dispersion of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(sampler_priors=dict(velocity_dispersion="velocity_dispersion_bernardi"))
        >>> print(od.sample_velocity_dispersion(size=10))
        """

        vd_inv_cdf = self.vd_inv_cdf.copy()
        # get the interpolator (inverse cdf) and sample
        if get_attribute:
            return njit(lambda size: inverse_transform_sampler(size, vd_inv_cdf[0], vd_inv_cdf[1]))
        else:
            return inverse_transform_sampler(size, vd_inv_cdf[0], vd_inv_cdf[1])

    def velocity_dispersion_ewoud(self, size, zl, get_attribute=False, **kwargs):
        """
        Function to sample velocity dispersion (redshift dependent) from Wempe et al. (2022). This uses inverse transform sampling.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        zl : `float`
            redshift of the lens galaxy
        get_attribute : `bool`
            if True, returns a function that can be used to sample velocity dispersion

        Returns
        -------
        sigma : `numpy.ndarray` (1D array of floats)
            velocity dispersion of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(sampler_priors=dict(velocity_dispersion="velocity_dispersion_ewoud"))
        >>> print(od.sample_velocity_dispersion(size=10, zl=0.5))
        """

        z_max = self.z_max
        z_min = self.z_min
        zlist = np.linspace(z_min, z_max, 100)
        vd_inv_cdf = self.vd_inv_cdf.copy()

        @njit
        def sampler(size, zl, inv_cdf, zlist):
            
            idx = np.searchsorted(zlist, zl)
            return inverse_transform_sampler(size, inv_cdf[idx][0], inv_cdf[idx][1])
        
        if get_attribute:
            return njit(lambda size, zl: sampler(size, zl, vd_inv_cdf, zlist))
        else:
            return sampler(size, zl, vd_inv_cdf, zlist)

    def cross_section_SIS(self, sigma, zl, zs):
        """
        Function to compute the SIS cross-section

        Parameters
        ----------
        sigma : `float`
            velocity dispersion of the lens galaxy
        zl : `float`
            redshift of the lens galaxy
        zs : `float`
            redshift of the source galaxy

        Returns
        -------
        cross_section : `float`
            SIS cross-section

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.cross_section_SIS(sigma=200., zl=0.5, zs=1.0))
        """
        zl = np.array([zl]).reshape(-1)
        zs = np.array([zs]).reshape(-1)
        Ds = self.angular_diameter_distance(zs)
        Dls = self.angular_diameter_distance_z1z2(zl, zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        return np.pi * theta_E**2
    
    def tau_zl_zs(self, zl, zs):
        """
        Function to compute the optical depth for a given lens redshift and source redshift

        Parameters
        ----------
        zl : `float`
            redshift of the lens galaxy
        zs : `float`
            redshift of the source galaxy

        Returns
        -------
        tau : `float`
            optical depth

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.tau_zl_zs(zl=0.5, zs=1.0))
        """
        zl = np.array([zl]).reshape(-1)
        zs = np.array([zs]).reshape(-1)
        # size=5000 will take ~ 48s to run, with ewoud vd sampler
        try:
            sigma = self.sample_velocity_dispersion(size=5000)
        except:
            sigma = self.sample_velocity_dispersion(size=5000, zl=zl)

        q = self.sample_axis_ratio(sigma)  # if SIS, q=array of 1.0
        no = 8*1e-3*self.cosmo.h**3
        test = phi_cut_SIE(q)*self.cross_section_SIS(sigma=sigma, zl=zl, zs=zs)/(4*np.pi)*no*self.differential_comoving_volume(zl)
        # average
        return np.mean(test)
    
    def optical_depth_calculator(self, zs):
        """
        Function to compute the optical depth without multiprocessing. This is the integrated version of tau_zl_zs from z=0 to z=zs.

        Parameters
        ----------
        zs : `float`
            source redshifts

        Returns
        -------
        tau : `float`
            optical depth

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.optical_depth_calculator(zs=1.0))
        """

        zs = np.array([zs]).reshape(-1)
        tau_list = []
        for z in zs:
            tau_list.append(quad(self.tau_zl_zs, 0, z, args=(z))[0])

        return np.array(tau_list)
    
    def optical_depth_multiprocessing(self, zs):
        """
        Function to compute the optical depth with multiprocessing. This is the integrated version of optical depth from z=0 to z=zs.

        Parameters
        ----------
        zs : `float`
            source redshifts

        Returns
        -------
        tau : `float`
            optical depth
            
        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.optical_depth_multiprocessing(zs=1.0))
        """

        zs = np.array([zs]).reshape(-1)
        no = 8*1e-3*self.cosmo.h**3
        vd_inv_cdf = self.vd_inv_cdf
        splinedVcdz = self.splinedVcdz
        splineDa = self.splineDa
        idx = np.arange(len(zs))
        try:
            zl_list = self.zl_list
            input_params = [(zs[i], no, vd_inv_cdf, splinedVcdz, splineDa, idx[i], zl_list) for i in range(len(zs))]
        except:
            input_params = [(zs[i], no, vd_inv_cdf, splinedVcdz, splineDa, idx[i]) for i in range(len(zs))]

        # Create a pool of workers and parallelize the integration
        with Pool(processes=self.npool) as pool:
            result = list(pool.map(self.tau_mp_routine, input_params))

        result = np.array(result)
        tau_list = result[:,1][np.array(result[:,0], dtype=int)]

        return tau_list

    def optical_depth_SIS_haris(self, zs):
        """
        Function to compute the strong lensing optical depth (SIS). \n
        LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.

        Parameters
        ----------
        zs : `float`
            source redshifts

        Returns
        -------
        tau : `float`
            strong lensing optical depth

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.optical_depth_SIS_haris(zs=1.0))
        """

        # z to luminosity_distance (luminosity_distance) conversion
        Dc = self.z_to_Dc(zs) * 1e-3  # 1e-3 converts Mpc to Gpc

        return (Dc / 62.2) ** 3  # 62.2 is the critical density in Gpc
    
    def create_lookup_table_fuction(self, z_max):
        """
        Functions to create lookup tables
        1. Redshift to co-moving distance.
        2. Co-moving distance to redshift.
        3. Redshift to angular diameter distance.

        Parameters
        ----------
        z_max : `float`
            maximum redshift of the lens galaxy population
        """

        Dc = lambda z_: self.cosmo.comoving_distance(z_).value  # co-moving distance in Mpc
        resolution = self.c_n_i["z_to_Dc"]["resolution"]
        create_new = self.c_n_i["z_to_Dc"]["create_new"]
        splineDc = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo, resolution=resolution),
            directory=self.directory,
            sub_directory="z_to_Dc", 
            name="z_to_Dc",
            x = np.geomspace(0.001, z_max, resolution),
            pdf_func= Dc, 
            conditioned_y=None,
            dimension=1,
            category="function",
            create_new= create_new,
        )
        self.splineDc = splineDc
        self.z_to_Dc = njit(lambda z_: cubic_spline_interpolator(z_, splineDc[0], splineDc[1]))

        # for co-moving distance to redshift
        resolution = self.c_n_i["Dc_to_z"]["resolution"]
        create_new = self.c_n_i["Dc_to_z"]["create_new"]
        splineDcInv = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo, resolution=resolution), 
            directory=self.directory,
            sub_directory="Dc_to_z", 
            name="Dc_to_z",
            x = np.geomspace(0.001, z_max, resolution),
            pdf_func= Dc, 
            conditioned_y=None, 
            dimension=1,
            category="function_inverse",
            create_new=create_new,
        )
        self.splineDcInv = splineDcInv
        self.Dc_to_z = njit(lambda Dc_: cubic_spline_interpolator(Dc_, splineDcInv[0], splineDcInv[1]))

        # for angular diameter distance
        resolution = self.c_n_i["angular_diameter_distance"]["resolution"]
        create_new = self.c_n_i["angular_diameter_distance"]["create_new"]
        Da = lambda z_: self.cosmo.angular_diameter_distance(z_).value
        splineDa = interpolator_from_pickle(
                    param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo, resolution=resolution), 
                    directory=self.directory,
                    sub_directory="angular_diameter_distance", 
                    name="angular_diameter_distance",
                    x = np.geomspace(0.001, z_max, resolution),
                    pdf_func= Da, 
                    conditioned_y=None, 
                    dimension=1,
                    category="function",
                    create_new=create_new,
                )
        self.splineDa = splineDa
        angular_diameter_distance = njit(lambda z_: cubic_spline_interpolator(z_, splineDa[0], splineDa[1]))
        self.angular_diameter_distance = angular_diameter_distance

        # for angular diameter distance between two redshifts
        self.angular_diameter_distance_z1z2 = njit(lambda zl0, zs0: (angular_diameter_distance(zs0)*(1.+zs0) - angular_diameter_distance(zl0)*(1.+zl0))/(1.+zs0))

        # get differential co-moving volume interpolator
        resolution = self.c_n_i["differential_comoving_volume"]["resolution"]
        create_new = self.c_n_i["differential_comoving_volume"]["create_new"]
        splinedVcdz = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo, resolution=resolution), 
            directory=self.directory,
            sub_directory="differential_comoving_volume", 
            name="differential_comoving_volume",
            x = np.geomspace(0.001, z_max, resolution),
            pdf_func= lambda z_: self.cosmo.differential_comoving_volume(z_).value * 4 * np.pi,
            conditioned_y=None, 
            dimension=1,
            category="function",
            create_new=create_new,
        )
        self.splinedVcdz = splinedVcdz
        self.differential_comoving_volume = njit(lambda z_: cubic_spline_interpolator(z_, splinedVcdz[0], splinedVcdz[1]))

    @property
    def strong_lensing_optical_depth(self):
        """
        Function to compute the strong lensing optical depth.

        Parameters
        ----------
        zs : `numpy.ndarray` (1D array of floats)
            source redshifts

        Returns
        -------
        tau : `numpy.ndarray` (1D array of floats)
            strong lensing optical depth

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.strong_lensing_optical_depth(np.array([0.1,0.2,0.3])))
        """

        return self._strong_lensing_optical_depth

    @strong_lensing_optical_depth.setter
    def strong_lensing_optical_depth(self, input_):
        if callable(input_):
            self._strong_lensing_optical_depth = input_
        else:
            # input can be function or spline interpolator array
            try:
                self._strong_lensing_optical_depth = njit(lambda zs: cubic_spline_interpolator(zs, input_[0], input_[1]))
            except:
                raise ValueError("strong_lensing_optical_depth must be a callable function or spline interpolator array.")

    @property
    def sample_velocity_dispersion(self):
        """
        Function to sample velocity dispersion. `zl` is required only if velocity dispersion sampler is redshift dependent.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        zl : `float`
            redshift of the lens galaxy

        Returns
        -------
        sigma : `numpy.ndarray` (1D array of floats)
            velocity dispersion of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.sample_velocity_dispersion(size=10))
        """

        return self._sample_velocity_dispersion

    @sample_velocity_dispersion.setter
    def sample_velocity_dispersion(self, prior):
        try:
            self._sample_velocity_dispersion = getattr(self, prior)(size=None, zl=None, get_attribute=True)
        except:
            self._sample_velocity_dispersion = prior

    @property
    def sample_axis_ratio(self):
        """
        Function to sample axis ratio from rayleigh distribution with given velocity dispersion.

        Parameters
        ----------
        sigma : `numpy.ndarray` (1D array of floats)
            velocity dispersion of the lens galaxy

        Returns
        -------
        q : `numpy.ndarray` (1D array of floats)
            axis ratio of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.sample_axis_ratio(sigma=200.))
        """

        return self._sample_axis_ratio

    @sample_axis_ratio.setter
    def sample_axis_ratio(self, prior):
        try:
            args = self.sampler_priors_params["axis_ratio"]
            self._sample_axis_ratio = getattr(self, prior)(sigma=None, get_attribute=True, param=args)
        except:
            self._sample_axis_ratio = prior

    @property
    def available_velocity_dispersion_list_and_its_params(self):
        """
        Function to list all available velocity dispersion sampler and its parameters.
        """
        self._available_velocity_dispersion_list_and_its_params = dict(
            velocity_dispersion_haris=dict(a=2.32 / 2.67, c=2.67),
            velocity_dispersion_gengamma=dict(a=2.32 / 2.67, c=2.67),
            velocity_dispersion_bernardi=None,
            velocity_dispersion_ewoud=None,
        )

        return self._available_velocity_dispersion_list_and_its_params
    
    @property
    def available_axis_ratio_list_and_its_params(self):
        """
        Function to list all available axis ratio sampler.
        """
        self._available_axis_ratio_list_and_its_params = dict(
            axis_ratio_rayleigh=None,
            axis_ratio_SIS=None,
        )

        return self._available_axis_ratio_list_and_its_params
    
    @property
    def available_optical_depth_list_and_its_params(self):
        """
        Function to list all available optical depth sampler.
        """
        self._available_optical_depth_list_and_its_params = dict(
            optical_depth_SIS_haris=None,
            optical_depth_SIS_hemanta=None,
            optical_depth_SIE_hemanta=None,
        )
        return self._available_optical_depth_list_and_its_params








