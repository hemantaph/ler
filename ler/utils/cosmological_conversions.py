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
    geomspace_part=False,
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
    x_min : float
        Minimum value of the grid.
    x_max : float
        Maximum value of the grid.
    resolution : int
        Total number of grid points.
    power_law_part : str, optional
        Which part of the grid should follow the power-law spacing. Options: 'lower' or 'upper'. Default is 'lower'.
    geomspace_part : bool or str, optional
        If `False`, keep the existing linear + power-law behavior. If `'lower'` or `'upper'`,
        replace that segment with geometric spacing while keeping the other segment linear.
        Geometric spacing is only used when the selected segment endpoints are strictly positive;
        otherwise the function falls back to the standard mixed-grid construction. Default is `False`.
    spacing_trend : str, optional
        Whether the power-law spacing should be increasing or decreasing. Options: 'increasing' or 'decreasing'. Default is 'increasing'.
    power : float, optional
        The power-law exponent. Higher values lead to more extreme spacing. Default is 2.3.
    value_transition_fraction : float, optional
        The fraction of the total value range at which to transition from linear to power-law spacing. Must be between 0 and 1. Default is 0.6.
    num_transition_fraction : float, optional
        The fraction of the total number of points at which to transition from linear to power-law spacing. Must be between 0 and 1. Default is 0.8.
    auto_match_slope : bool, optional
        Whether to automatically adjust the power-law exponent to match the slope of the linear spacing at the transition point. Default is True.
        This is ignored for the geometric-spacing segment when `geomspace_part` is used.

    Returns
    -------
    numpy.ndarray
        The generated grid points.

    Examples
    --------
    from ler.utils import generate_mixed_grid

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
    # linear+geomspace with geometric spacing in the upper segment
    x = generate_mixed_grid(
        x_min=0.1, x_max=10.0, resolution=resolution,
        geomspace_part='lower',
        value_transition_fraction=0.3,
        num_transition_fraction=0.6,
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
    x_trans = x_min + u_trans * (x_max - x_min)

    if geomspace_part not in (False, None, 'lower', 'upper'):
        raise ValueError("geomspace_part must be False, 'lower' or 'upper'")

    if geomspace_part == 'lower' and x_min > 0.0 and x_trans > 0.0:
        x_low = np.geomspace(x_min, x_trans, n_low)
        x_high = np.linspace(x_trans, x_max, n_high + 1)[1:]
        return np.concatenate([x_low, x_high])

    if geomspace_part == 'upper' and x_trans > 0.0 and x_max > 0.0:
        x_low = np.linspace(x_min, x_trans, n_low)
        x_high = np.geomspace(x_trans, x_max, n_high + 1)[1:]
        return np.concatenate([x_low, x_high])

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
    Function to create a lookup table for luminosity distance as a function of redshift.

    The interpolated quantity is

    .. math::

        D_L(z) = (1+z) D_C(z),

    as returned by ``astropy.cosmology`` for the supplied cosmology.

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
        Object of FunctionConditioning class containing luminosity distance as a
        function of redshift.
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

    .. math::

        D_L(z) = (1+z) D_C(z).

    Parameters
    ----------
    zs : `numpy.ndarray` or `float`
        Source redshifts

    Returns
    -------
    luminosity_distance : `numpy.ndarray`
        luminosity distance in Mpc

    Examples
    --------
    >>> import numpy as np
    >>> from ler.utils import luminosity_distance
    >>> dl = luminosity_distance(get_attribute=True)
    >>> distances = dl.function(np.array([1., 2.]))
    >>> redshift = dl.function_inverse(np.array([100., 200.]))
    """

    return luminosity_distance_object if get_attribute else luminosity_distance_object(z)


def differential_comoving_volume(z=None, z_min=0.001, z_max=10., cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory="./interpolator_json", create_new=False, resolution=500, get_attribute=True):
    """
    Create a FunctionConditioning object for the differential comoving volume dVc/dz.

    The stored table is full-sky:

    .. math::

        \\frac{dV_c}{dz} = 4\\pi \\frac{dV_c}{dz\\,d\\Omega}.

    Parameters
    ----------
    z : ``float`` or ``numpy.ndarray`` or ``None``
        Redshift(s) at which to evaluate. If None, returns the FunctionConditioning object.
    z_min : ``float``
        Minimum redshift. default: 0.001
    z_max : ``float``
        Maximum redshift. default: 10.0
    cosmo : ``astropy.cosmology``
        Cosmology object. default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    directory : ``str``
        Directory for storing interpolator JSON files. default: './interpolator_json'
    create_new : ``bool``
        If True, create new interpolator. default: False
    resolution : ``int``
        Number of grid points for the interpolator. default: 500
    get_attribute : ``bool``
        If True, return the FunctionConditioning object. default: True

    Returns
    -------
    differential_comoving_volume : ``FunctionConditioning`` or ``numpy.ndarray``
        dVc/dz in Mpc^3 sr^-1 (multiplied by 4*pi for full sky).
    """

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

    The returned values are full-sky ``dVc/dz`` in ``Mpc^3`` per unit redshift.

    Parameters
    ----------
    zs : `numpy.ndarray` or `float`
        Source redshifts

    Returns
    -------
    differential_comoving_volume : `numpy.ndarray`
        differential comoving volume in Mpc^3 per unit redshift.

    Examples
    --------
    >>> import numpy as np
    >>> from ler.utils import differential_comoving_volume
    >>> dvc_dz = differential_comoving_volume(get_attribute=True)
    >>> values = dvc_dz.function(np.array([1., 2.]))
    """

    return differential_comoving_volume_object if get_attribute else differential_comoving_volume_object(z)

def comoving_distance(z=None, z_min=0.001, z_max=10., cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory="./interpolator_json", create_new=False, resolution=500, get_attribute=True):
    """
    Create a FunctionConditioning object for the comoving distance.

    .. math::

        D_C(z) = c \\int_0^z \\frac{dz'}{H(z')}.

    Parameters
    ----------
    z : ``float`` or ``numpy.ndarray`` or ``None``
        Redshift(s) at which to evaluate. If None, returns the FunctionConditioning object.
    z_min : ``float``
        Minimum redshift. default: 0.001
    z_max : ``float``
        Maximum redshift. default: 10.0
    cosmo : ``astropy.cosmology``
        Cosmology object. default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    directory : ``str``
        Directory for storing interpolator JSON files. default: './interpolator_json'
    create_new : ``bool``
        If True, create new interpolator. default: False
    resolution : ``int``
        Number of grid points for the interpolator. default: 500
    get_attribute : ``bool``
        If True, return the FunctionConditioning object. default: True

    Returns
    -------
    comoving_distance : ``FunctionConditioning`` or ``numpy.ndarray``
        Comoving distance in Mpc.
    """

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

    .. math::

        D_C(z) = c \\int_0^z \\frac{dz'}{H(z')}.

    Parameters
    ----------
    zs : `numpy.ndarray` or `float`
        Source redshifts

    Returns
    -------
    comoving_distance : `numpy.ndarray`
        comoving distance in Mpc

    Examples
    --------
    >>> import numpy as np
    >>> from ler.utils import comoving_distance
    >>> dc = comoving_distance(get_attribute=True)
    >>> distances = dc.function(np.array([1., 2.]))
    >>> redshift = dc.function_inverse(np.array([100., 200.]))
    """

    return comoving_distance_object if get_attribute else comoving_distance_object(z)

def angular_diameter_distance(z=None, z_min=0.001, z_max=10., cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory="./interpolator_json", create_new=False, resolution=500, get_attribute=True):    
    """
    Create a FunctionConditioning object for the angular diameter distance.

    .. math::

        D_A(z) = \\frac{D_C(z)}{1+z}.

    Parameters
    ----------
    z : ``float`` or ``numpy.ndarray`` or ``None``
        Redshift(s) at which to evaluate. If None, returns the FunctionConditioning object.
    z_min : ``float``
        Minimum redshift. default: 0.001
    z_max : ``float``
        Maximum redshift. default: 10.0
    cosmo : ``astropy.cosmology``
        Cosmology object. default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    directory : ``str``
        Directory for storing interpolator JSON files. default: './interpolator_json'
    create_new : ``bool``
        If True, create new interpolator. default: False
    resolution : ``int``
        Number of grid points for the interpolator. default: 500
    get_attribute : ``bool``
        If True, return the FunctionConditioning object. default: True

    Returns
    -------
    angular_diameter_distance : ``FunctionConditioning`` or ``numpy.ndarray``
        Angular diameter distance in Mpc.
    """

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

    .. math::

        D_A(z) = \\frac{D_C(z)}{1+z}.

    Parameters
    ----------
    zs : `numpy.ndarray` or `float`
        Source redshifts

    Returns
    -------
    angular_diameter_distance : `numpy.ndarray`
        angular diameter distance in Mpc

    Examples
    --------
    >>> import numpy as np
    >>> from ler.utils import angular_diameter_distance
    >>> da = angular_diameter_distance(get_attribute=True)
    >>> distances = da.function(np.array([1., 2.]))
    """

    return angular_diameter_distance_object if get_attribute else angular_diameter_distance_object(z)

def angular_diameter_distance_z1z2(z1=None, z2=None, z_min=0.001, z_max=10., cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory="./interpolator_json", create_new=False, resolution=500, get_attribute=True):
    """
    Create a FunctionConditioning object for the angular diameter distance between two redshifts.

    Uses the relation

    .. math::

        D_A(z_1, z_2) =
        \\frac{D_A(z_2)(1+z_2) - D_A(z_1)(1+z_1)}{1+z_2}.

    Parameters
    ----------
    z1 : ``float`` or ``numpy.ndarray`` or ``None``
        Lens redshift(s).
    z2 : ``float`` or ``numpy.ndarray`` or ``None``
        Source redshift(s).
    z_min : ``float``
        Minimum redshift. default: 0.001
    z_max : ``float``
        Maximum redshift. default: 10.0
    cosmo : ``astropy.cosmology``
        Cosmology object. default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    directory : ``str``
        Directory for storing interpolator JSON files. default: './interpolator_json'
    create_new : ``bool``
        If True, create new interpolator. default: False
    resolution : ``int``
        Number of grid points for the interpolator. default: 500
    get_attribute : ``bool``
        If True, return the FunctionConditioning object. default: True

    Returns
    -------
    angular_diameter_distance_z1z2 : ``FunctionConditioning`` or ``numpy.ndarray``
        Angular diameter distance between z1 and z2 in Mpc.
    """

    z_min = 0.001 if z_min == 0. else z_min

    angular_diameter_distance_object = angular_diameter_distance(z_min=z_min, z_max=z_max, cosmo=cosmo, directory=directory, create_new=create_new, resolution=resolution, get_attribute=get_attribute)

    # for angular diameter distance between two redshifts
    _Da = angular_diameter_distance_object.function
    @njit()
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
    Angular diameter distance between two redshifts.

    .. math::

        D_A(z_l, z_s) =
        \\frac{D_A(z_s)(1+z_s) - D_A(z_l)(1+z_l)}{1+z_s}.

    Parameters
    ----------
    zl0 : `numpy.ndarray` or `float`
        Lens redshifts
    zs0 : `numpy.ndarray` or `float`
        Source redshifts

    Returns
    -------
    angular_diameter_distance_z1z2 : `numpy.ndarray`
        angular diameter distance between ``zl0`` and ``zs0`` in Mpc.

    Examples
    --------
    >>> import numpy as np
    >>> from ler.utils import angular_diameter_distance_z1z2
    >>> da_z1z2 = angular_diameter_distance_z1z2(get_attribute=True)
    >>> distances = da_z1z2.function(np.array([0.5, 1.0]), np.array([1.0, 2.0]))
    """

    return angular_diameter_distance_z1z2_object if get_attribute else angular_diameter_distance_z1z2_object(z1, z2)
