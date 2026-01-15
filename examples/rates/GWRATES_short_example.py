# This is a short example to simulate binary black hole mergers and calculate their rates ($yr^{-1}$) and finally compare the results.
# GWRATES settings and important results (like simulated events file names and location, rates etc.) will be saved in the `./ler_data/gwrates_params.json` file.
# The GW events will be saved in the `./ler_data/gw_param.json` file and the detectable unlensed events will be saved in the `./ler_data/gw_param_detectable.json` file.
# A plot of the redshift distributions of GW events will be saved in the current directory.

# call the GWRATES class
from ler.rates import GWRATES
# other necessary imports
from astropy.cosmology import LambdaCDM

ler = GWRATES(
    # GWRATES setup arguments
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
    interpolator_directory='./interpolator_json', # directory to store the interpolator pickle files. 'ler' uses interpolation to get values of various functions to speed up the calculations (relying on numba njit).
    ler_directory='./ler_data', # directory to store all the outputs
    verbose=True, # if True, will print all information at initialization

    # CBCSourceParameterDistribution class arguments
    source_priors= {'merger_rate_density': 'merger_rate_density_bbh_oguri2018', 'source_frame_masses': 'binary_masses_BBH_powerlaw_gaussian', 'zs': 'sample_source_redshift', 'geocent_time': 'sampler_uniform', 'ra': 'sampler_uniform', 'dec': 'sampler_cosine', 'phase': 'sampler_uniform', 'psi': 'sampler_uniform', 'theta_jn': 'sampler_sine'},
    source_priors_params= {'merger_rate_density': {'R0': 2.39e-08, 'b2': 1.6, 'b3': 2.0, 'b4': 30}, 'source_frame_masses': {'mminbh': 4.98, 'mmaxbh': 112.5, 'alpha': 3.78, 'mu_g': 32.27, 'sigma_g': 3.88, 'lambda_peak': 0.03, 'delta_m': 4.8, 'beta': 0.81}, 'zs': None, 'geocent_time': {'min_': 1238166018, 'max_': 1269702018}, 'ra': {'min_': 0.0, 'max_': 6.283185307179586}, 'dec': None, 'phase': {'min_': 0.0, 'max_': 6.283185307179586}, 'psi': {'min_': 0.0, 'max_': 3.141592653589793}, 'theta_jn': None},
    spin_zero= True, # if True, spins will be set to zero
    spin_precession= False, # if True, spins will be precessing

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
    interpolator_dir = './interpolator_json',
    gwsnr_verbose = False,
    multiprocessing_verbose = True,
    mtot_cut = True,

    # common arguments, to generate interpolator
    create_new_interpolator = dict(
        redshift_distribution=dict(create_new=False, resolution=1000),
        z_to_luminosity_distance=dict(create_new=False, resolution=1000),
        differential_comoving_volume=dict(create_new=False, resolution=1000),
        Dl_to_z=dict(create_new=False, resolution=1000),
    )
)

# simulate events, save it
gw_params = ler.gw_cbc_statistics(size=100000, resume=False)
# select detectable events and save it.
# calculate rates.
rate, param_detectable = ler.gw_rate()


# plot the redshift distribution of GW events
import matplotlib.pyplot as plt
from ler.utils import plots as lerplt


# input param_dict can be either a dictionary or a json file name that contains the parameters
plt.figure(figsize=(6, 4))
lerplt.param_plot(
    param_name='zs',
    param_dict='ler_data/gw_param_detectable.json',
    plot_label='zs (detectable)',
)
lerplt.param_plot(
    param_name='zs',
    param_dict='ler_data/gw_param.json',
    plot_label='zs (all)',
)
plt.xlim(0.001,10)
plt.grid(alpha=0.4)
plt.xlabel(r'$z_s$')
plt.ylabel('pdf')
print('\n Saving the plot...')
plt.savefig('redshift_distribution_gwrates.png')