import numpy as np
from quintet import Quintet
from scipy.stats import norm
import pylab as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad
from astropy.cosmology import Planck18
import json
from lens_galaxy_population import LensGalaxyPopulation
from source_population import CompactBinaryPopulation

# Conversions from SI units to CGS units
C = 299792458.
G = 6.67408*1e-11
Mo = 1.989*1e30
Mpc = 3.086*1e22
# Define the maximum number of images
max_images = 5
min_images = 2

class LeR():
    def __init__(self, npool=int(4), z_min=0.0001, z_max=10., m_min=4.59, m_max=86.22, \
                 event_type = "StellarBBH", equal_mass=False, \
                 SNR_finder='quintet', **kwargs_snr_function):
        '''
        class for rate calculation
        class initialization
        Input parameters:
            npool: number of cores to use
            z_min: minimum redshift
            z_max: maximum redshift
            m_min: minimum mass of the binary binary black hole
            m_max: maximum mass of the binary binary black hole
            event_type: type of event, either "BBH" or "BNS"
            equal_mass: if True, the binary black holes will have equal masses
            SNR_finder: if 'quintet', the SNR will be calculated using the quintet package
                        if 'custom', the SNR will be calculated using a custom function
            kwargs_snr_function: keyword arguments for the SNR function
        Output parameters:
            None
        '''
        self.z_min      = z_min
        self.z_max      = z_max
        self.m_min      = m_min
        self.m_max      = m_max
        self.event_type = event_type
        self.gw_param   = False          # this needed not repeat sampling
        self.lensed_param   = False          # this needed not repeat sampling
        self.npool      = npool
        # initializing function for fast SNR calculation
        self.equal_mass = equal_mass
        
        # initializing function for fast SNR calculation (in case of quintet)
        if SNR_finder=='quintet':
            # default
            self.snr_ = self.quintet_intialization(kwargs_snr_function)
        else:
            # custom SNR function
            self.snr_ = None
        # extra note on how to change snr finder function
        # self.snr_ = custom_snr_finder_function()     

        # Create lookup tables
        self.create_lookup_tables(z_min, z_max)
        # initialization of clasess
        # CompactBinaryPopulation already inherits from Source_Galaxy_Population_Model class form source_population.py
        self.CompactBinaryPopulation = CompactBinaryPopulation(z_min=z_min, z_max=z_max, m_min=m_min, m_max=m_max,\
                                               event_type = event_type)
        
        self.LensGalaxyPopulation = LensGalaxyPopulation(self.CompactBinaryPopulation)
        
        # both UnlensedCBCStatistics and LensedCBCStatistics sampler needs this
        '''
        self.gw_param_sampler_dict = {'zs':False, 'm_min':m_min, 'm_max':m_max, 'z_min':z_min, 'z_max':z_max, \
                         'event_type':event_type,'model_pars_gwcosmo': {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, \
                                                                   'mmin': m_min, 'mmax': m_max, 'lambda_peak': 0.08, \
                                                                    'mu_g': 33.07, 'sigma_g': 5.69} }
        '''
        # only LensedCBCStatistics sampler needs this
        self.lens_param_sampler_dict = { 'min_lensed_images':min_images, 'max_lensed_images':max_images, \
                                        'lensModelList':['EPL_NUMBA', 'SHEAR']}


    def quintet_intialization(self, kwargs_snr_function):  
        '''
        function for initializing the quintet package
        Input parameters:
            kwargs_snr_function: keyword arguments for the SNR function
        Output Parameters:
            snr_: the quintet object
                    quintet object is used to calculate the SNR and pdet (probability of detection)
        '''
        quintet_param_dict = dict(nsamples_mtot=100, nsamples_mass_ratio=50,\
                                      sampling_frequency=4096.,\
                                      waveform_approximant = "IMRPhenomD", minimum_frequency = 20., \
                                      snr_type = 'interpolation', waveform_inspiral_must_be_above_fmin=True, \
                                      psds=False, psd_file=False)

        for key, value in kwargs_snr_function.items():
            quintet_param_dict[key] = value

        snr_ = Quintet(npool=self.npool, mtot_min=self.m_min*2, mtot_max=self.m_max*2, nsamples_mtot=quintet_param_dict['nsamples_mtot'],\
                            nsamples_mass_ratio=quintet_param_dict['nsamples_mass_ratio'],\
                            sampling_frequency=quintet_param_dict['sampling_frequency'],\
                            waveform_approximant = quintet_param_dict['waveform_approximant'], \
                            minimum_frequency = quintet_param_dict['minimum_frequency'], \
                            snr_type = quintet_param_dict['snr_type'],\
                            waveform_inspiral_must_be_above_fmin=quintet_param_dict['waveform_inspiral_must_be_above_fmin'],\
                            psds = quintet_param_dict['psds'],\
                            psd_file = quintet_param_dict['psd_file'])
        return snr_
            

    def create_lookup_tables(self, z_min, z_max):
        '''
        function for creating lookup tables for fast calculation
        Intput Parameters:
            z_min: minimum redshift
            z_max: maximum redshift
        Output Parameters:
            None
        '''
        self.z_min      = z_min
        self.z_max      = z_max
        self.gw_param   = False          # this needed not repeat sampling 
        
        # initialing cosmological functions for fast calculation through interpolation
        z               = np.linspace(0,z_max,500) # red-shift
        Dc              = Planck18.comoving_distance(z).value # co-moving distance in Mpc
        luminosity_distance              = Planck18.luminosity_distance(z).value # luminosity distance in Mpc
        self.z_to_Dc    = interp1d( z, Dc, kind = 'cubic')
        self.Dc_to_z    = interp1d( Dc, z, kind = 'cubic')
        self.z_to_luminosity_distance    = interp1d( z, luminosity_distance, kind = 'cubic')
        # for angular diameter distance
        quad_ = [] # refer to quad integradtion from scipy 
        for ii in range(len(z)):
            quad_.append(quad(Planck18._inv_efunc_scalar, 0., z[ii], args=Planck18._inv_efunc_scalar_args)[0])
        quad_ = np.array(quad_)
        self.quad_int = interp1d(z, np.array(quad_), kind = 'cubic')
        # Lookup table for differential comoving distance
        differential_comoving_volume = Planck18.differential_comoving_volume(z).value*4*np.pi # differential comoving volume in Mpc^3
        self.differential_comoving_volume = interp1d(z, differential_comoving_volume, kind = 'cubic')


    def UnlensedCBCStatistics(self, nsamples=1000, jsonFile=True, **kwargs):
        '''
        function to generate unlensed GW source parameters
        Intput Parameters:
            nsamples: number of samples
            snr_threshold: snr threshold of detection
            jsonFile: if True, store all gravitational waves source parameters in json file 
                        (for all sources, detected and undetected)
            kwargs: if new paramteres are provided, it will be used for sampling source parameters
        Output Parameters:
            unlensed_gw_params: dictionary of unlensed GW source parameters
        '''                
        # gw parameter sampling
        gw_param = self.CompactBinaryPopulation.sample_gw_parameters(nsamples=nsamples)

        # with quintet
        # calculate signal to noise ratio and probability of detection for each gw parameters set
        print('snr calculation started...')
        snrs = self.snr_.snr(GWparam_dict = gw_param)
        gw_param.update(snrs)
        print('snr calculation ended...')
        
        self.gw_param = gw_param
        # store all params in json file
        if jsonFile:
            file_name = './gw_params.json'
            json_dump = json.dumps(gw_param, cls=NumpyEncoder)
            with open(file_name, "w") as write_file:
                json.dump(json.loads(json_dump), write_file, indent=4)
        
        return None
    
    def unlensed_rate(self, size=1000, snr_threshold=8., jsonFile=True):
        '''
        function to calculate unlensed merger rate
        Intput Parameters:
            size: number of samples
            snr_threshold: threshold for detection signal to noise ratio
            jsonFile: if True, store all gravitational waves source parameters in json file 
                        (for detected sources)
        Output Parameters:
            unlensed_rate: unlensed merger rate in yr^-1
        '''
        if self.gw_param == False:   
            # sample the source parameters
            print('gw_param not sampled beforehand. Sampling now...')
            self.UnlensedCBCStatistics(nsamples=size, snr_threshold=snr_threshold, jsonFile=jsonFile)
        else:
            print('already sampled gw_param found.')
            size = len(self.gw_param['zs'])
            print('sample size will be taken as that gw_param, size=',size)
        
        gw_param = self.gw_param.copy()
        snr = gw_param['opt_snr_net']
        pdet = 1 - norm.cdf(snr_threshold - snr)
        gw_param['pdet_net'] = pdet
        
        # selecting only detectable
        idx_detectable = snr > snr_threshold
        
        # store all detectable params in json file
        for key, value in gw_param.items():
            gw_param[key] = value[idx_detectable]
        if jsonFile:
            file_name = './gw_params_detectable.json'
            json_dump = json.dumps(gw_param, cls=NumpyEncoder)
            with open(file_name, "w") as write_file:
                json.dump(json.loads(json_dump), write_file, indent=4)
        
        # montecarlo integration
        # The total rate is Eq. A4 of https://arxiv.org/pdf/2106.06303.pdf
        # R = C0 int Theta(rho-rhoc) p(z) p(theta) dtheta dz_s, where C0 = int R(zs)/(1+zs) dVc/dzs dzs is the normalization constant for p(z)
        # Thus R = C0 <Theta(rho-rhoc)>
        c0 = self.CompactBinaryPopulation.normalization_pdf_z
        total_rate_step = c0 * np.mean(idx_detectable)
        print("total unlensed rate with step function: {}".format(total_rate_step))
        
        # with pdet 
        total_rate_pdet = c0 * np.mean(pdet)
        print("total unlensed rate with pdet function: {}".format(total_rate_pdet))
        
        return(total_rate_step,total_rate_pdet)
    
    def LensedCBCStatistics(self, nsamples=1000, jsonFile=True, **kwargs):
        '''
        function to generate lensed GW source parameters, lens parameters and image parameters
        Intput Parameters:
            nsamples: number of samples
            snr_threshold: threshold for detection signal to noise ratio
            jsonFile: if True, store lensed GW source parameters, lens parameters and image parameters in json file
                        (both for detected and undetected sources)
            **kwargs: if new parameters are provided, it will be used for sampling
        Output Parameters:
            lensed_param: dictionary of lensed GW source parameters, lens parameters and image parameters
        '''
        lens_sampler_dict = self.lens_param_sampler_dict
        
        # if new paramteres are provided
        for key, value in kwargs.items():
            if key in lens_sampler_dict:
                lens_sampler_dict[key] = value
        
        # gw_param will not be kept same as that of unlensed case. So, it is sampled newly
        lensed_param = self.LensGalaxyPopulation.sample_lens_parameters( size=nsamples )
        # now get (strongly lensed) image paramters along with lens parameters 
        lensed_param = self.LensGalaxyPopulation.get_image_properties(lens_parameters=lensed_param, \
                                            lensModelList=lens_sampler_dict['lensModelList'], npool=self.npool)
        
        #---------------------------------------------------#
        #      Get all of the SNRs and Pdet 
        #---------------------------------------------------#
        # Get all of the signal to noise ratios
        print('snr calculation started...')
        snrs = self.LensGalaxyPopulation.get_lensed_snrs(snr_calculator=self.snr_, lensed_param=lensed_param)
        lensed_param.update(snrs)
        print('snr calculation ended...')

        self.lensed_param = lensed_param
        # store all params in json file
        if jsonFile:
            file_name = './lensed_params.json'
            json_dump = json.dumps(lensed_param, cls=NumpyEncoder)
            with open(file_name, "w") as write_file:
                json.dump(json.loads(json_dump), write_file, indent=4)
                
        return(lensed_param)
                
    def lensed_rate(self, size=1000, snr_threshold=8., num_img=2, jsonFile=True, none_as_nan=True):
        '''
        Function to calculate detectable lensed merger rate
        Intput Parameters:
            size (int): number of samples
            snr_threshold (float/array): threshold for detection signal to noise ratio
            num_img (int/array): number of images
                                e.g. For Sub-thershold events, snr_threshold=[8.,6.], num_img=[1,1]
                                The event will contain 1 image with snr>8 and 1 image with snr>6
            jsonFile (bool): if True, store all gravitational waves source parameters in json file
            none_as_nan (bool): if True,  no value is kept as np.nan
                                if False, no value is kept as 0.
        Output Parameters:
            lensed_rate (float): lensed merger rate in yr^-1
        '''
        if self.lensed_param == False:   
            # sample the source parameters
            print('lensed_param not sampled beforehand. Sampling now...')
            self.LensedCBCStatistics(nsamples=size, jsonFile=jsonFile)
        else:
            print('already sampled lensed_param found.')
            size = len(self.lensed_param['zs'])
            print('sample size will be taken as that lensed_param, size=',size)
        
        lensed_param = self.lensed_param.copy()

        snr  = lensed_param['opt_snr_net'] # Dimensions are (nsamples, n_max_images)
        
        snr_threshold, num_img = np.array([snr_threshold]).reshape(-1),np.array([num_img]).reshape(-1) # convert to array
        # sort in descending order of each row
        argTH = (-snr_threshold).argsort() 
        sortedSNR = -np.sort(-snr,axis=1)
        num1 = 0 # tracks the number of images for the current threshold
        num2 = 0 # tracks the column number of the already sorted snr 2D array
        snr_hit = np.full(len(snr),True) # boolean array to store the result of the threshold condition
        pdet_combined = np.full(len(snr),1.) # array to store the result of the pdet condition
        for i in argTH:
            num1 = num_img[i]
            for j in range(num1):
                # snr_hit step function case
                snr_hit = snr_hit&(sortedSNR[:,num2]>snr_threshold[i])
                # pdet for probability of detection
                pdet = 1 - norm.cdf(snr_threshold[i] - np.nan_to_num(sortedSNR[:,num2]))
                pdet_combined = pdet_combined*pdet
                num2 += 1
        lensed_param['pdet_net'] = pdet

        # store all params in json file
        if none_as_nan==True:
            for key, value in lensed_param.items():
                lensed_param[key] = value[snr_hit]
        else:
            for key, value in lensed_param.items():
                lensed_param[key] = np.nan_to_num(value[snr_hit])
            
        # store all params in json file
        if jsonFile:
            file_name = './lensed_params_detectable.json'
            json_dump = json.dumps(lensed_param, cls=NumpyEncoder)
            with open(file_name, "w") as write_file:
                json.dump(json.loads(json_dump), write_file, indent=4)

    
        # montecarlo integration
        # The total rate is Eq. A4 of https://arxiv.org/pdf/2106.06303.pdf
        # R = C0 int Theta(rho-rhoc) p(z) p(theta) dtheta dz_s, where C0 = int R(zs)/(1+zs) dVc/dzs tau(zs) dzs is the normalization constant for p(z)
        # Thus R = C0 <Theta(rho-rhoc)>
        c0 = self.LensGalaxyPopulation.normalization_pdf_z
        total_rate_step = c0 * np.mean(snr_hit)
        print("total unlensed rate with step function: {}".format(total_rate_step))

        # with pdet 
        total_rate_pdet = c0 * np.mean(pdet_combined)
        print("total unlensed rate with pdet function: {}".format(total_rate_pdet))

        return(total_rate_step,total_rate_pdet)

   
class NumpyEncoder(json.JSONEncoder):
    '''
    class for storing a numpy.ndarray or any nested-list composition as JSON file
    '''
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
    