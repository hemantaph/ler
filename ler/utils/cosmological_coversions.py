import numpy as np
from numba import njit
from .function_interpolation import FunctionConditioning
# for redshift to luminosity distance conversion
from astropy.cosmology import LambdaCDM


def generate_mixed_grid(
    x_min,
    x_max,
    resolution,
    power_law_part='lower',
    spacing_trend='increasing',
    power=2.3,
    value_transition_fraction=0.6,
    num_transition_fraction=0.8,
    auto_match_slope=True,
):
    """
    Generalized mixed spacing grid generator. Safely handles negative ranges.

    Parameters
    ----------
    x_min, x_max : float
        The physical domain boundaries (can be negative, crossing zero, etc.).
    power_law_part : str
        'lower' to apply power-law to the first segment.
        'upper' to apply power-law to the second segment.
    spacing_trend : str
        'increasing': Step sizes grow as x increases.
        'decreasing': Step sizes shrink as x increases.
    power : float
        Exponent (>= 1.0) to control curvature. Overridden if auto_match_slope=True.
    auto_match_slope : bool
        If True, analytically calculates the required power AND trend to guarantee 
        a smooth transition without a kink.

    Returns
    -------
    numpy.ndarray
        The generated grid points.

    Examples
    --------
    from ler.utils.cosmological_conversions import generate_mixed_grid

    resolution=20
    # linear+power-law with power-law in the upper segment and decreasing step sizes
    x = generate_mixed_grid(
        x_min=0.0, x_max=10.0, resolution=resolution,
        power_law_part='upper',
        spacing_trend='decreasing',  # Forces largest steps near z_trans
        power=2.5,
        value_transition_fraction=0.6,
        num_transition_fraction=0.3,
        auto_match_slope=True       # We accept the kink to control the exact power
    )
    # powerlaw+linear with power-law in the lower segment and increasing step sizes
    x = generate_mixed_grid(
        x_min=0.0, x_max=10.0, resolution=resolution,
        power_law_part='lower',
        spacing_trend='increasing',  # Forces largest steps near z_trans
        power=2.5,
        value_transition_fraction=0.3,
        num_transition_fraction=0.6,
        auto_match_slope=True       # We accept the kink to control the exact power
    )
    """
    if x_max <= x_min:
        return np.linspace(x_min, x_max, resolution)

    # 1. Normalized transition parameter
    u_trans = float(np.clip(value_transition_fraction, 0.0, 1.0))
    
    if u_trans <= 0.0 or u_trans >= 1.0:
        return np.linspace(x_min, x_max, resolution)

    if num_transition_fraction is None:
        num_transition_fraction = u_trans

    n_low = max(2, int(resolution * num_transition_fraction))
    n_low = min(n_low, resolution - 1)
    n_high = resolution - n_low

    # 2. Build the grid in normalized [0, 1] space
    if power_law_part == 'lower':
        du_lin = (1.0 - u_trans) / n_high
        ratio = du_lin / u_trans
        N_int = n_low - 1

        if auto_match_slope and N_int > 1:
            if ratio >= 1.0:
                power, spacing_trend = 1.0, 'increasing'
            elif ratio > 1.0 / N_int:
                spacing_trend = 'increasing'
                power = np.log(1.0 - ratio) / np.log((N_int - 1) / N_int)
            else:
                spacing_trend = 'decreasing'
                power = np.log(ratio) / np.log(1.0 / N_int)

        t = np.linspace(0.0, 1.0, n_low)
        if spacing_trend == 'increasing':
            u_low = u_trans * (t ** power)
        else:
            u_low = u_trans * (1.0 - (1.0 - t) ** power)
        
        u_high = np.linspace(u_trans, 1.0, n_high + 1)[1:]
        u_grid = np.concatenate([u_low, u_high])

    elif power_law_part == 'upper':
        du_lin = u_trans / (n_low - 1)
        ratio = du_lin / (1.0 - u_trans)
        N_int = n_high

        if auto_match_slope and N_int > 1:
            if ratio >= 1.0:
                power, spacing_trend = 1.0, 'increasing'
            elif ratio < 1.0 / N_int:
                spacing_trend = 'increasing'
                power = np.log(ratio) / np.log(1.0 / N_int)
            else:
                spacing_trend = 'decreasing'
                power = np.log(1.0 - ratio) / np.log((N_int - 1) / N_int)

        t = np.linspace(0.0, 1.0, n_high + 1)[1:]
        if spacing_trend == 'increasing':
            u_high = u_trans + (1.0 - u_trans) * (t ** power)
        else:
            u_high = u_trans + (1.0 - u_trans) * (1.0 - (1.0 - t) ** power)

        u_low = np.linspace(0.0, u_trans, n_low)
        u_grid = np.concatenate([u_low, u_high])

    else:
        raise ValueError("power_law_part must be 'lower' or 'upper'")

    # 3. Map the normalized grid [0, 1] to the physical domain [x_min, x_max]
    return x_min + u_grid * (x_max - x_min)

def luminosity_distance(z=None, z_min=0.001, z_max=10., cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory="./interpolator_json", create_new=False, resolution=500, get_attribute=True):
    """
    Function to create a lookup table for the luminosity distance wrt redshift.

    Parameters
    ----------
    z : `numpy.ndarray` or `float`
        Source redshifts
    z_min : `float`
        Minimum redshift of the source population
    z_max : `float`
        Maximum redshift of the source population

    Attributes
    ----------
    z_to_luminosity_distance : `ler.utils.FunctionConditioning`
        Object of FunctionConditioning class containing the luminosity distance wrt redshift
    """

    z_min = 0.001 if z_min == 0. else z_min

    zs = generate_mixed_grid(z_min, z_max, resolution)
    Dl = cosmo.luminosity_distance(zs).value
    
    luminosity_distance_object = FunctionConditioning(
        function=Dl,
        x_array=zs,
        conditioned_y_array=None,
        identifier_dict=dict(z_min=z_min, z_max=z_max, cosmology=cosmo, resolution=resolution, details="luminosity_distance from astropy.cosmology"),
        directory=directory,
        sub_directory="luminosity_distance",
        name="luminosity_distance",
        create_new=create_new,
        create_function_inverse=True,
        create_function=True,
        create_pdf=False,
        create_rvs=False,
        callback='function',
    )
    luminosity_distance_object.__doc__ = """
    Redshift to luminosity distance conversion.

    Parameters
    ----------
    zs : `numpy.ndarray` or `float`
        Source redshifts

    Returns
    ----------
    luminosity_distance : `numpy.ndarray`
        luminosity distance in Mpc

    Examples
    ----------
    >>> from ler.gw_source_population import SourceGalaxyPopulationModel
    >>> ler = SourceGalaxyPopulationModel()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    >>> luminosity_distance = ler.luminosity_distance(1.)
    >>> luminosity_distance = ler.luminosity_distance.function(np.array([1., 2.]))
    >>> redshift = ler.luminosity_distance.function_inverse(np.array([100., 200.]))
    """

    return luminosity_distance_object if get_attribute else luminosity_distance_object(z)


def differential_comoving_volume(z=None, z_min=0.001, z_max=10., cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory="./interpolator_json", create_new=False, resolution=500, get_attribute=True):
        
    z_min = 0.001 if z_min == 0. else z_min

    # get differential co-moving volume interpolator
    zs = generate_mixed_grid(z_min, z_max, resolution)
    dVcdz = cosmo.differential_comoving_volume(zs).value * 4 * np.pi  # volume of shell in Mpc^3
    differential_comoving_volume_object = FunctionConditioning(
        function=dVcdz,
        x_array=zs,
        conditioned_y_array=None,
        identifier_dict=dict(z_min=z_min, z_max=z_max, cosmology=cosmo, resolution=resolution, details="differential_comoving_volume from astropy.cosmology"),
        directory=directory,
        sub_directory="differential_comoving_volume",
        name="differential_comoving_volume",
        create_new=create_new,
        create_function_inverse=False,
        create_function=True,
        create_pdf=False,
        create_rvs=False,
        callback='function',
    )

    differential_comoving_volume_object.__doc__ = """
    Redshift to differential comoving volume conversion.

    Parameters
    ----------
    zs : `numpy.ndarray` or `float`
        Source redshifts

    Returns
    ----------
    differential_comoving_volume : `numpy.ndarray`
        differential comoving volume in Mpc^3

    Examples
    ----------
    >>> from ler.len_galaxy_population import OpticalDepth
    >>> ler = OpticalDepth()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    >>> differential_comoving_volume = ler.differential_comoving_volume(1.)
    >>> differential_comoving_volume = ler.differential_comoving_volume.function(np.array([1., 2.]))
    """

    return differential_comoving_volume_object if get_attribute else differential_comoving_volume_object(z)

def comoving_distance(z=None, z_min=0.001, z_max=10., cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory="./interpolator_json", create_new=False, resolution=500, get_attribute=True):

    z_min = 0.001 if z_min == 0. else z_min
    zs = generate_mixed_grid(z_min, z_max, resolution)

    Dc = cosmo.comoving_distance(zs).value  # co-moving distance in Mpc
    comoving_distance_object = FunctionConditioning(
        function=Dc,
        x_array=zs,
        conditioned_y_array=None,
        identifier_dict=dict(
            z_min=z_min,
            z_max=z_max,
            cosmology=cosmo,
            resolution=resolution,
            details="comoving_distance from astropy.cosmology",
        ),
        directory=directory,
        sub_directory="comoving_distance",
        name="comoving_distance",
        create_new=create_new,
        create_function_inverse=True,
        create_function=True,
        create_pdf=False,
        create_rvs=False,
        callback="function",
    )
    comoving_distance_object.__doc__ = """
    Redshift to comoving distance conversion.

    Parameters
    ----------
    zs : `numpy.ndarray` or `float`
        Source redshifts

    Returns
    ----------
    comoving_distance : `numpy.ndarray`
        comoving distance in Mpc

    Examples
    ----------
    >>> from ler.len_galaxy_population import OpticalDepth
    >>> ler = OpticalDepth()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    >>> comoving_distance = ler.comoving_distance(1.)
    >>> comoving_distance = ler.comoving_distance.function(np.array([1., 2.]))
    >>> redshift = ler.comoving_distance.function_inverse(np.array([100., 200.]))
    """

    return comoving_distance_object if get_attribute else comoving_distance_object(z)

def angular_diameter_distance(z=None, z_min=0.001, z_max=10., cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory="./interpolator_json", create_new=False, resolution=500, get_attribute=True):    

    z_min = 0.001 if z_min == 0. else z_min
    zs = generate_mixed_grid(z_min, z_max, resolution)

    Da = cosmo.angular_diameter_distance(zs).value
    angular_diameter_distance_object = FunctionConditioning(
        function=Da,
        x_array=zs,
        conditioned_y_array=None,
        identifier_dict=dict(
            z_min=z_min,
            z_max=z_max,
            cosmology=cosmo,
            resolution=resolution,
            details="angular_diameter_distance from astropy.cosmology",
        ),
        directory=directory,
        sub_directory="angular_diameter_distance",
        name="angular_diameter_distance",
        create_new=create_new,
        create_function_inverse=False,
        create_function=True,
        create_pdf=False,
        create_rvs=False,
        callback="function",
    )

    angular_diameter_distance_object.__doc__ = """
    Redshift to angular diameter distance conversion.

    Parameters
    ----------
    zs : `numpy.ndarray` or `float`
        Source redshifts

    Returns
    ----------
    angular_diameter_distance : `numpy.ndarray`
        angular diameter distance in Mpc

    Examples
    ----------
    >>> from ler.len_galaxy_population import OpticalDepth
    >>> ler = OpticalDepth()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    >>> angular_diameter_distance = ler.angular_diameter_distance(1.)
    >>> angular_diameter_distance = ler.angular_diameter_distance.function(np.array([1., 2.]))
    >>> redshift = ler.angular_diameter_distance.function_inverse(np.array([100., 200.]))
    """

    return angular_diameter_distance_object if get_attribute else angular_diameter_distance_object(z)

def angular_diameter_distance_z1z2(z1=None, z2=None, z_min=0.001, z_max=10., cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory="./interpolator_json", create_new=False, resolution=500, get_attribute=True):

    z_min = 0.001 if z_min == 0. else z_min

    angular_diameter_distance_object = angular_diameter_distance(z_min=z_min, z_max=z_max, cosmo=cosmo, directory=directory, create_new=create_new, resolution=resolution, get_attribute=get_attribute)

    # for angular diameter distance between two redshifts
    _Da = angular_diameter_distance_object.function
    @njit(cache=True)
    def angular_diameter_distance_z1z2(zl0, zs0):
        return (_Da(zs0) * (1.0 + zs0) - _Da(zl0) * (1.0 + zl0)) / (1.0 + zs0)

    angular_diameter_distance_z1z2_object = FunctionConditioning(
        function=None,
        x_array=None,
        identifier_dict=dict(
            z_min=z_min,
            z_max=z_max,
            cosmology=cosmo,
            resolution=resolution,
            details="angular_diameter_distance_z1z2 from astropy.cosmology",
        ),
        create_function=angular_diameter_distance_z1z2,
        callback="function",
    )

    angular_diameter_distance_z1z2_object.__doc__ = """
    Redshift to angular diameter distance conversion.

    Parameters
    ----------
    zl0 : `numpy.ndarray` or `float`
        Lens redshifts
    zs0 : `numpy.ndarray` or `float`
        Source redshifts

    Returns
    ----------
    angular_diameter_distance_z1z2 : `numpy.ndarray`
        angular diameter distance in Mpc

    Examples
    ----------
    >>> from ler.len_galaxy_population import OpticalDepth
    >>> ler = OpticalDepth()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    >>> angular_diameter_distance_z1z2 = ler.angular_diameter_distance_z1z2(1., 2.)
    >>> angular_diameter_distance_z1z2 = ler.angular_diameter_distance_z1z2.function(np.array([1., 2.]), np.array([1., 2.]))
    """

    return angular_diameter_distance_z1z2_object if get_attribute else angular_diameter_distance_z1z2_object(z1, z2)

