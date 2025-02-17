# This is a short example to simulate lensed and unlensed binary black hole mergers and calculate their rates ($yr^{-1}$) and finally compare the results.
# LeR settings and important results (like simulated events file names and location, rates and rate ratios, etc.) will be saved in the `./ler_data/ler_params.json` file.
# The unlensed events will be saved in the `./ler_data/unlensed_param.json` file and the detectable unlensed events will be saved in the `./ler_data/unlensed_param_detectable.json` file.
# The lensed events will be saved in the `./ler_data/lensed_param.json` file and the detectable lensed events will be saved in the `./ler_data/lensed_param_detectable.json` file.
# A plot of the redshift distributions of lensed and unlensed events will be saved in the current directory.

# call the LeR class
from ler.rates import LeR

# other necessary imports
from astropy.cosmology import LambdaCDM

# define the main function to avoid multiprocessing issues
def main():
    """
    Function to simulate lensed and unlensed binary black hole mergers and calculate their rates ($yr^{-1}$) and finally compare the results. The outputs are saved in the `ler_data` directory by default.
    """

    # if you want to initialize with the default values, you can simply call the class without any arguments; ler = LeR()
    # below I am showing an example with all the arguments, including the arguments of the parent classes (CBCSourceParameterDistribution, LensGalaxyParameterDistribution, ImageProperties) and the arguments of the gwsnr package.
    ler = LeR(
        # LeR setup arguments
        npool=4, # number of processors to use
        z_min=0.0, # minimum redshift
        z_max=10.0, # maximum redshift
        event_type='BBH', # event type
        size=100000, # number of events to simulate
        batch_size=50000, # batch size
        cosmology=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7), # cosmology
        snr_finder=None, # snr calculator from 'gwsnr' package will be used
        pdet_finder=None,  # will not be consider unless specified
        list_of_detectors=None, # list of detectors that will be considered when calculating snr or pdet for lensed events. if None, all the detectors from 'gwsnr' will be considered
        json_file_names=dict(
            ler_params="ler_params.json", # to store initialization parameters and important results
            unlensed_param="unlensed_param.json", # to store all unlensed events
            unlensed_param_detectable="unlensed_param_detectable.json", # to store only detectable unlensed events
            lensed_param="lensed_param.json", # to store all lensed events 
            lensed_param_detectable="lensed_param_detectable.json"), # to store only detectable lensed events
        interpolator_directory='./interpolator_pickle', # directory to store the interpolator pickle files. 'ler' uses interpolation to get values of various functions to speed up the calculations (relying on numba njit).
        ler_directory='./ler_data', # directory to store all the outputs
        verbose=True, # if True, will print all information at initialization

        # CBCSourceParameterDistribution class arguments
        source_priors= {'merger_rate_density': 'merger_rate_density_bbh_popI_II_oguri2018', 'source_frame_masses': 'binary_masses_BBH_popI_II_powerlaw_gaussian', 'zs': 'sample_source_redshift', 'geocent_time': 'sampler_uniform', 'ra': 'sampler_uniform', 'dec': 'sampler_cosine', 'phase': 'sampler_uniform', 'psi': 'sampler_uniform', 'theta_jn': 'sampler_sine'},
        source_priors_params= {'merger_rate_density': {'R0': 2.39e-08, 'b2': 1.6, 'b3': 2.0, 'b4': 30}, 'source_frame_masses': {'mminbh': 4.98, 'mmaxbh': 112.5, 'alpha': 3.78, 'mu_g': 32.27, 'sigma_g': 3.88, 'lambda_peak': 0.03, 'delta_m': 4.8, 'beta': 0.81}, 'zs': None, 'geocent_time': {'min_': 1238166018, 'max_': 1269702018}, 'ra': {'min_': 0.0, 'max_': 6.283185307179586}, 'dec': None, 'phase': {'min_': 0.0, 'max_': 6.283185307179586}, 'psi': {'min_': 0.0, 'max_': 3.141592653589793}, 'theta_jn': None},
        spin_zero= True, # if True, spins will be set to zero
        spin_precession= False, # if True, spins will be precessing

        # LensGalaxyParameterDistribution class arguments
        lens_type = 'epl_shear_galaxy',
        lens_functions =  {'strong_lensing_condition': 'rjs_with_cross_section_SIE', 'optical_depth': 'optical_depth_SIE_hemanta', 'param_sampler_type': 'sample_all_routine'},
        lens_priors =  {'source_redshift_sl': 'strongly_lensed_source_redshifts', 'lens_redshift': 'lens_redshift_SDSS_catalogue', 'velocity_dispersion': 'velocity_dispersion_ewoud', 'axis_ratio': 'axis_ratio_rayleigh', 'axis_rotation_angle': 'axis_rotation_angle_uniform', 'external_shear': 'shear_norm', 'density_profile_slope': 'density_profile_slope_normal', 'source_parameters': 'sample_gw_parameters'},
        lens_priors_params =  {'source_redshift_sl': None, 'lens_redshift': None, 'velocity_dispersion': None, 'axis_ratio': {'q_min': 0.2, 'q_max': 1.0}, 'axis_rotation_angle': {'phi_min': 0.0, 'phi_max': 6.283185307179586}, 'external_shear': {'scale': 0.05}, 'density_profile_slope': {'mean': 2.0, 'std': 0.2}, 'source_parameters': None},

        # ImageProperties class arguments
        n_min_images = 2,
        n_max_images = 4,
        geocent_time_min = 1126259462.4,
        geocent_time_max = 1756979462.4,
        lens_model_list = ['EPL_NUMBA', 'SHEAR'],

        # gwsnr package arguments
        mtot_min = 2.0,
        mtot_max = 184.98599853446768,
        ratio_min = 0.1,
        ratio_max = 1.0,
        mtot_resolution = 500,
        ratio_resolution = 50,
        sampling_frequency = 2048.0,
        waveform_approximant = 'IMRPhenomD',
        minimum_frequency = 20.0,
        snr_type = 'interpolation',
        psds = {'L1':'aLIGO_O4_high_asd.txt','H1':'aLIGO_O4_high_asd.txt', 'V1':'AdV_asd.txt', 'K1':'KAGRA_design_asd.txt'},
        ifos = ['L1', 'H1', 'V1'],
        interpolator_dir = './interpolator_pickle',
        gwsnr_verbose = False,
        multiprocessing_verbose = True,
        mtot_cut = True,

        # common arguments, to generate interpolator
        create_new_interpolator = dict(
            redshift_distribution=dict(create_new=False, resolution=1000),
            z_to_luminosity_distance=dict(create_new=False, resolution=1000),
            velocity_dispersion=dict(create_new=False, resolution=1000),
            axis_ratio=dict(create_new=False, resolution=1000),
            optical_depth=dict(create_new=False, resolution=200),
            z_to_Dc=dict(create_new=False, resolution=1000),
            Dc_to_z=dict(create_new=False, resolution=1000),
            angular_diameter_distance=dict(create_new=False, resolution=1000),
            differential_comoving_volume=dict(create_new=False, resolution=1000),
            Dl_to_z=dict(create_new=False, resolution=1000),
        )
    )

    # run the simulation

    # unlensed events
    unlensed_param = ler.unlensed_cbc_statistics(
        size=100000, # number of events to simulate
        resume=False, # if True, will resume the simulation from the last saved batch
        save_batch=False, # if True, will save the batch in each iteration
        output_jsonfile=None, # if not None, file name from ler.json_file_names['unlensed_param'] will be used.
    )

    # lensed events
    lensed_param = ler.lensed_cbc_statistics(
        size=100000, # number of events to simulate
        resume=False, # if True, will resume the simulation from the last saved batch
        save_batch=False, # if True, will save the batch in each iteration
        output_jsonfile=None, # if not None, file name from ler.json_file_names['lensed_param'] will be used.
    )

    # select detectable events and save them.
    # calculate rates.
    # compare the rates.
    rate_ratio, unlensed_param_detectable, lensed_param_detectable = ler.rate_comparison_with_rate_calculation(
            unlensed_param=None,  # if None, will get get param wrt the unlensed_param file name from ler.json_file_names['unlensed_param']
            snr_threshold_unlensed=8.0,  # snr threshold for unlensed events
            output_jsonfile_unlensed=None, # if not None, file name from ler.json_file_names['unlensed_param_detectable'] will be used.
            lensed_param=None, # if None, will get get param wrt the lensed_param file name from ler.json_file_names['lensed_param']
            snr_threshold_lensed=[8.0,8.0], # snr threshold for lensed images
            num_img=[1,1], # number of images corresponding to the snr_threshold_lensed
            output_jsonfile_lensed=None, # if not None, file name from ler.json_file_names['lensed_param_detectable'] will be used.
            nan_to_num=True, # if True, will replace nan values with 0
            detectability_condition="step_function",
        )

def plot():
    """
    Function to plot the redshift distribution of lensed and unlensed events. The plot will be saved in the current directory.
    """

    # plot the results
    # redshift distribution only
    # quick plot
    from ler.utils import plots as lerplt
    import matplotlib.pyplot as plt

    # plotting the distribution of event parameters
    # comparision of redshift distribution for lensed and unlensed events
    # param_dict can be either a dictionary or a json file name that contains the parameters
    plt.figure(figsize=(6, 4))
    # for unlensed case
    lerplt.param_plot(
        param_name='zs',
        param_dict='ler_data/unlensed_param_detectable.json',
        plot_label='zs (bbh-detectable)',
        histogram=False,
        kde=True,
        kde_bandwidth=0.5,
    )
    lerplt.param_plot(
        param_name='zs',
        param_dict='ler_data/unlensed_param.json',
        plot_label='zs (bbh-all)',
        histogram=False,
        kde=True,
    )
    # for lensed case
    lerplt.param_plot(
        param_name='zs',
        param_dict='ler_data/lensed_param_detectable.json',
        plot_label='zs (bbh-lensed-detectable)',
        histogram=False,
        kde=True,
        kde_bandwidth=0.5,
    )
    lerplt.param_plot(
        param_name='zs',
        param_dict='ler_data/lensed_param.json',
        plot_label='zs (bbh-lensed-all)',
        histogram=False,
        kde=True,
    )
    plt.xlim(0.001,8)
    plt.grid(alpha=0.4)
    plt.xlabel('Source redshift (zs)')
    plt.ylabel('Probability Density')
    # save the plot
    print('Saving the plot...')
    plt.savefig('zs_distribution.png', dpi=300)

if __name__ == "__main__":
    main();

    # plot the results
    print('calling plot function to plot redshift distribution...')
    plot();