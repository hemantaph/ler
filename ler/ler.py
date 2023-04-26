import numpy as np
from gwsnr import GWSNR
from scipy.stats import norm
import pylab as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad
from astropy.cosmology import Planck18
import json
from ler.lens_galaxy_population import LensGalaxyPopulation
from ler.source_population import CompactBinaryPopulation
import contextlib

# Conversions from SI units to CGS units
C = 299792458.
G = 6.67408*1e-11
Mo = 1.989*1e30
Mpc = 3.086*1e22
# Define the maximum number of images
#------------------------------------------------------------#

#------------------------------------------------------------#
class LeR():
    def __init__(self, nsamples=1000, npool=int(4), z_min=0., z_max=10., m_min=4.59, m_max=86.22, \
                 event_type = "StellarBBH", equal_mass=False, batch_size = 50000, \
                 SNR_finder='gwsnr', **kwargs_snr_function):
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
            SNR_finder: if 'gwsnr', the SNR will be calculated using the gwsnr package
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
        self.batch_size = batch_size
        
        # initializing function for fast SNR calculation (in case of gwsnr)
        if SNR_finder=='gwsnr':
            # default
            self.snr_ = self.gwsnr_intialization(kwargs_snr_function)
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
        # LensGalaxyPopulation already inherits from Lens_Galaxy_Population_Model class form lens_galaxy_population.py
        self.LensGalaxyPopulation = LensGalaxyPopulation(self.CompactBinaryPopulation)

        # dictionary of params for sampler
        # this will be used for kde generation
        # for unlened case
        self.gw_param_sampler_dict = {'nsamples':nsamples, 'm_min':m_min, 'm_max':m_max, 'z_min':z_min, 'z_max':z_max, \
                            'event_type':event_type,'model_pars_gwcosmo': {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, \
                                                                    'mmin': m_min, 'mmax': m_max, 'lambda_peak': 0.08, \
                                                                    'mu_g': 33.07, 'sigma_g': 5.69} }
        # for lened case
        self.lensed_param_sampler_dict = {'nsamples':nsamples, 'min_lensed_images':2, 'max_lensed_images':4, \
                                        'lensModelList':['EPL_NUMBA', 'SHEAR']}
        self.store_ler_params()

        return None
    
    def store_ler_params(self):
        '''
        function to store the parameters of the LER model
        Input parameters:
            None
        Output parameters:
            None
        '''
        # store gw_param_sampler_dict, lensed_param_sampler_dict and snr_calculator_dict
        parameters_dict = {}
        parameters_dict.update({'gw_param_sampler_dict': self.gw_param_sampler_dict})
        parameters_dict.update({'lensed_param_sampler_dict': self.lensed_param_sampler_dict})
        parameters_dict.update({'snr_calculator_dict': self.snr_calculator_dict})
        file_name = './LeR_params.json'
        json_dump = json.dumps(parameters_dict, cls=NumpyEncoder)
        with open(file_name, "w") as write_file:
            json.dump(json.loads(json_dump), write_file, indent=4)

        return None


    def gwsnr_intialization(self, kwargs_snr_function):  
        '''
        function for initializing the gwsnr package
        Input parameters:
            kwargs_snr_function: keyword arguments for the SNR function
        Output Parameters:
            snr_: the gwsnr object
                    gwsnr object is used to calculate the SNR and pdet (probability of detection)
        '''
        gwsnr_param_dict = dict(nsamples_mtot=100, nsamples_mass_ratio=50,\
                                      sampling_frequency=4096.,\
                                      waveform_approximant = "IMRPhenomD", minimum_frequency = 20., \
                                      snr_type = 'interpolation', waveform_inspiral_must_be_above_fmin=True, \
                                      psds=False, psd_file=False)

        for key, value in kwargs_snr_function.items():
            gwsnr_param_dict[key] = value

        snr_ = GWSNR(npool=self.npool, mtot_min=self.m_min*2, mtot_max=self.m_max*2, nsamples_mtot=gwsnr_param_dict['nsamples_mtot'],\
                            nsamples_mass_ratio=gwsnr_param_dict['nsamples_mass_ratio'],\
                            sampling_frequency=gwsnr_param_dict['sampling_frequency'],\
                            waveform_approximant = gwsnr_param_dict['waveform_approximant'], \
                            minimum_frequency = gwsnr_param_dict['minimum_frequency'], \
                            snr_type = gwsnr_param_dict['snr_type'],\
                            waveform_inspiral_must_be_above_fmin=gwsnr_param_dict['waveform_inspiral_must_be_above_fmin'],\
                            psds = gwsnr_param_dict['psds'],\
                            psd_file = gwsnr_param_dict['psd_file'])
        
        self.snr_calculator_dict = gwsnr_param_dict
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


    def UnlensedCBCStatistics(self, nsamples=False, jsonFile=True, **kwargs):
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
        batch_size = self.batch_size
        # gw parameter sampling
        if nsamples==False:
            nsamples = self.gw_param_sampler_dict['nsamples']
        else:
            self.gw_param_sampler_dict['nsamples'] = nsamples
            
        # sampling in batches
        if nsamples%batch_size==0:
            num_batches = int(nsamples/batch_size)
        else:
            num_batches = int(nsamples/batch_size)+1
        
        print("chosen batch size = {}. If you want to change batch size, self.batch_size = new_size".format(batch_size))
        print("There will be {} batche(s)".format(num_batches))
        num_batches = num_batches-1 # consideing couting from 0
        if nsamples>batch_size:
            frac_batches = nsamples-num_batches*batch_size # 
        else:
            frac_batches = nsamples
        track_batches = 0
        #---------------------------------------------------#
        print("Batch no. {}".format(track_batches))
        # get gw params
        print('sampling gw params...')
        gw_param = self.CompactBinaryPopulation.sample_gw_parameters(nsamples=frac_batches)
        # Get all of the signal to noise ratios
        print('calculating snrs...')
        snrs = self.snr_.snr(GWparam_dict = gw_param)
        gw_param.update(snrs)
        #---------------------------------------------------#
        #---------------------------------------------------#
        gw_param_ = None
        for i in range(num_batches):
            track_batches += 1
            print("Batch no. {}".format(track_batches))
            print('sampling gw params...')
            gw_param_ = self.CompactBinaryPopulation.sample_gw_parameters(nsamples=batch_size)
            print('calculating snrs...')
            snrs = self.snr_.snr(GWparam_dict = gw_param_)
            gw_param_.update(snrs)
            for key in gw_param.keys():
                gw_param[key] = np.append(gw_param[key],gw_param_[key], axis=0)
        del gw_param_, snrs
        #---------------------------------------------------#
        
        self.gw_param = gw_param
        # store all params in json file
        if jsonFile:
            file_name = './gw_params.json'
            json_dump = json.dumps(gw_param, cls=NumpyEncoder)
            with open(file_name, "w") as write_file:
                json.dump(json.loads(json_dump), write_file, indent=4)
        
        return gw_param
    
    def unlensed_rate(self, size=False, snr_threshold=8., jsonFile=True):
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
            self.UnlensedCBCStatistics(nsamples=False, snr_threshold=snr_threshold, jsonFile=jsonFile);
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
        self.gw_param_detectable = gw_param
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
    
    def LensedCBCStatistics(self, nsamples=False, jsonFile=True, **kwargs):
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
        
        lens_sampler_dict = self.lensed_param_sampler_dict
        batch_size = self.batch_size
        # if new paramteres are provided
        if nsamples==False:
            nsamples = lens_sampler_dict['nsamples']
        else:
            lens_sampler_dict['nsamples'] = nsamples
        for key, value in kwargs.items():
            if key in lens_sampler_dict:
                lens_sampler_dict[key] = value
        min_img = int(lens_sampler_dict['min_lensed_images'])
        max_img = int(lens_sampler_dict['max_lensed_images'])

        # gw_param will not be kept same as that of unlensed case. So, it is sampled newly
        # sampling in batches
        if nsamples%batch_size==0:
            num_batches = int(nsamples/batch_size)
        else:
            num_batches = int(nsamples/batch_size)+1
        
        print("chosen batch size = {}. If you want to change batch size, self.batch_size = new_size".format(batch_size))
        print("There will be {} batche(s)".format(num_batches))
        num_batches = num_batches-1 # consideing couting from 0
        if nsamples>batch_size:
            frac_batches = nsamples-num_batches*batch_size # 
        else:
            frac_batches = nsamples

        track_batches = 0
        #---------------------------------------------------#
        print("Batch no. {}".format(track_batches))
        lensed_param = self.LensGalaxyPopulation.sample_lens_parameters( size= frac_batches)
        # now get (strongly lensed) image paramters along with lens parameters 
        lensed_param = \
        self.LensGalaxyPopulation.get_image_properties(n_min_images=min_img, n_max_images=max_img,\
                                                       lens_parameters=lensed_param,\
                                                       lensModelList=lens_sampler_dict['lensModelList'], npool=self.npool)
        # Get all of the signal to noise ratios
        print('calculating snrs...')
        snrs = self.LensGalaxyPopulation.get_lensed_snrs(snr_calculator=self.snr_, lensed_param=lensed_param, n_max_images=max_img)
        lensed_param.update(snrs)
        #---------------------------------------------------#
        #---------------------------------------------------#
        lensed_param_ = None
        for i in range(num_batches):
            track_batches += 1
            print("Batch no. {}".format(track_batches))
            #try:
            lensed_param_ = self.LensGalaxyPopulation.sample_lens_parameters( size=batch_size )
            # now get (strongly lensed) image paramters along with lens parameters 
            lensed_param_ = \
            self.LensGalaxyPopulation.get_image_properties(n_min_images=min_img, n_max_images=max_img,\
                                                            lens_parameters=lensed_param_,\
                                                            lensModelList=lens_sampler_dict['lensModelList'], npool=self.npool)
            print('calculating snrs...')
            snrs = self.LensGalaxyPopulation.get_lensed_snrs(snr_calculator=self.snr_, lensed_param=lensed_param_, n_max_images=max_img)
            lensed_param_.update(snrs)
            for key in lensed_param.keys():
                lensed_param[key] = np.append(lensed_param[key],lensed_param_[key], axis=0)
        del lensed_param_, snrs
        #---------------------------------------------------#
        
        self.lensed_param = lensed_param
        # store all params in json file
        if jsonFile:
            file_name = './lensed_params.json'
            json_dump = json.dumps(lensed_param, cls=NumpyEncoder)
            with open(file_name, "w") as write_file:
                json.dump(json.loads(json_dump), write_file, indent=4)
                
        return(lensed_param)
                
    def lensed_rate(self, size=False, snr_threshold=8., num_img=2, jsonFile=True, none_as_nan=True):
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
            self.LensedCBCStatistics(nsamples=size, jsonFile=jsonFile);
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
        weights = lensed_param['weights']
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
        self.lensed_param_detectable = lensed_param
        # montecarlo integration
        # The total rate is Eq. A4 of https://arxiv.org/pdf/2106.06303.pdf
        # R = C0 int Theta(rho-rhoc) p(z) p(theta) dtheta dz_s, where C0 = int R(zs)/(1+zs) dVc/dzs tau(zs) dzs is the normalization constant for p(z)
        # Thus R = C0 <Theta(rho-rhoc)>
        c0 = self.LensGalaxyPopulation.normalization_pdf_z
        total_rate_step = c0 * np.mean(snr_hit*weights)
        print("total unlensed rate with step function: {}".format(total_rate_step))

        # with pdet 
        total_rate_pdet = c0 * np.mean(pdet_combined*weights)
        print("total unlensed rate with pdet function: {}".format(total_rate_pdet))

        return(total_rate_step,total_rate_pdet)
    
    def selecting_n_lensed_detectable_events_from_dict(self, snr_threshold=8., num_img=2, none_as_nan=True, \
                                      lenstype='I'):
        '''
        Function to select n lensed detectable events from self.lensed_param
        Input Parameters:
            snr_threshold (float/array): threshold for detection signal to noise ratio
            num_img (int/array): number of images
                                e.g. For Sub-thershold events, snr_threshold=[8.,6.], num_img=[1,1]
                                The event will contain 1 image with snr>8 and 1 image with snr>6
            none_as_nan (bool): if True,  no value is kept as np.nan
                                if False, no value is kept as 0.
            lenstype (str): lens type, 'I' or 'II'  
        Output Parameters:
            lensed_param (dict): dictionary of lensed parameters
        '''
        try:
            size = len(self.lensed_param['zs'])
        except:
            print('lensed_param not sampled beforehand.')
            return None
        param = self.lensed_param.copy()
        snr  = param['opt_snr_net'] # Dimensions are (nsamples, n_max_images)
        
        # lens type selection dictionary
        lens_type_dict = {'I':[0,1],'II':[2,3]}
        if lenstype=='I':
            snr = snr[:,lens_type_dict[lenstype][0]:lens_type_dict[lenstype][1]+1]
        elif lenstype=='II':
            snr = snr[:,lens_type_dict[lenstype][0]:lens_type_dict[lenstype][1]+1]
        else:
            print('lens type not found. Please choose from',lens_type_dict.keys())
            return None
        
        # dealing with snr_threshold and num_img
        snr_threshold, num_img = np.array([snr_threshold]).reshape(-1),np.array([num_img]).reshape(-1) # convert to array
        # sort in descending order of each row
        argTH = (-snr_threshold).argsort() 
        sortedSNR = -np.sort(-snr,axis=1)
        num1 = 0 # tracks the number of images for the current threshold
        num2 = 0 # tracks the column number of the already sorted snr 2D array
        snr_hit = np.full(len(snr),True) # boolean array to store the result of the threshold condition
        for i in argTH:
            num1 = num_img[i]
            for j in range(num1):
                # snr_hit step function case
                snr_hit = snr_hit&(sortedSNR[:,num2]>snr_threshold[i])
                num2 += 1
        
        if none_as_nan==True:
            for key, value in param.items():
                param[key] = value[snr_hit]
        else:
            for key, value in param.items():
                param[key] = np.nan_to_num(value[snr_hit])
        return param
    
    def selecting_n_lensed_detectable_events_with_sampling(self, snr_threshold=8., num_img=2, none_as_nan=True, \
                                                         size=100, lenstype='I', min_img=2, max_img=4):
        '''
        Function to select n lensed detectable events with sampling
        '''
        n=0
        param_final = {}
        while(n<size):
            # disable print statements
            with contextlib.redirect_stdout(None):
                param = self.LensedCBCStatistics(nsamples=self.batch_size, jsonFile=False, min_lensed_images=min_img, \
                                                 max_lensed_images=max_img)
                snr  = param['opt_snr_net'] # Dimensions are (nsamples, n_max_images)
                
                # lens type selection dictionary
                lens_type_dict = {'I':[0,1],'II':[2,3]}
                if lenstype=='I':
                    snr = snr[:,lens_type_dict[lenstype][0]:lens_type_dict[lenstype][1]+1]
                elif lenstype=='II':
                    snr = snr[:,lens_type_dict[lenstype][0]:lens_type_dict[lenstype][1]+1]
                elif lenstype=='any':
                    pass
                else:
                    print('lens type not found. Please choose from',lens_type_dict.keys())
                    return None
                
                # dealing with snr_threshold and num_img
                snr_threshold, num_img = np.array([snr_threshold]).reshape(-1),np.array([num_img]).reshape(-1) # convert to array
                # sort in descending order of each row
                argTH = (-snr_threshold).argsort() 
                sortedSNR = -np.sort(-snr,axis=1)
                num1 = 0 # tracks the number of images for the current threshold
                num2 = 0 # tracks the column number of the already sorted snr 2D array
                snr_hit = np.full(len(snr),True) # boolean array to store the result of the threshold condition
                for i in argTH:
                    num1 = num_img[i]
                    for j in range(num1):
                        # snr_hit step function case
                        snr_hit = snr_hit&(sortedSNR[:,num2]>snr_threshold[i])
                        num2 += 1
                # update the final param dictionary
                if none_as_nan==True:
                    for key, value in param.items():
                        # if param_final[key] is empty
                        if n==0:
                            param_final[key] = value[snr_hit]
                        else:
                            # concatenate the values of the current param dictionary with the final param dictionary
                            param_final[key] =  np.concatenate((param_final[key],value[snr_hit]), axis=0)
                else:
                    for key, value in param.items():
                        # if param_final[key] is empty
                        if n==0:
                            param_final[key] = value[snr_hit]
                        else:
                            param_final[key] = np.concatenate((param_final[key],np.nan_to_num(value[snr_hit])), axis=0)
                n = len(param_final['opt_snr_net'])
            print('collected number of events = ',n)

        # trim the final param dictionary
        for key, value in param_final.items():
            param_final[key] = value[:size]
        print( 'After trimming, collected number of events = ',len(param_final['opt_snr_net']) )
        return param_final
    
    def selecting_n_unlensed_detectable_events_from_dict(snr_threshold=8.):
        '''
        Function to select n lensed detectable events from self.gw_param
        '''
        try:
            size = len(self.gw_param['zs'])
        except:
            print('unlensed_param not sampled beforehand.')
            return None
        param = self.gw_param.copy()
        snr  = param['opt_snr_net'] # Dimensions are (nsamples)
        
        # selecting only detectable
        idx_detectable = snr > snr_threshold
        # store all detectable params in json file
        for key, value in param.items():
            param[key] = value[idx_detectable]
        return param

    def selecting_n_unlensed_detectable_events_with_sampling(self, snr_threshold=8., size=100):
        '''
        Function to select n lensed detectable events with sampling
        '''
        n=0
        param_final = {}
        while(n<size):
            # disable print statements
            with contextlib.redirect_stdout(None):
                param = self.UnlensedCBCStatistics(nsamples=self.batch_size, jsonFile=False)
                snr  = param['opt_snr_net']

                # selecting only detectable
                idx_detectable = snr > snr_threshold
                # update the final param dictionary
                for key, value in param.items():
                    # param_final is empty
                    if n==0:
                        param_final[key] = value[idx_detectable]
                    else:
                        # concatenate the values of the current param dictionary with the final param dictionary
                        param_final[key] =  np.concatenate((param_final[key],value[idx_detectable]), axis=0)
                n = len(param_final['opt_snr_net'])
            print('collected number of events = ',n)

        # trim the final param dictionary
        for key, value in param_final.items():
            param_final[key] = value[:size]
        print( 'After trimming, collected number of events = ',len(param_final['opt_snr_net']) )
        return param_final

    def selecting_n_detectable_events(self, snr_threshold=8., num_img=2, none_as_nan=True, lensed=True, jsonFile=True, \
                                      lenstype='I', new=False, size=100, min_img=2, max_img=4, batch_size=10000, **kwargs):
        '''
        Function to select n detectable events
        Input parameters:
            snr_threshold (float): the threshold for the SNR
            num_img (int): the number of images
            none_as_nan (bool): if True, then replace None with np.nan
            lensed (bool): if True, then select lensed events
            jsonFile (bool): if True, then save the dictionary as a json file
            lenstype (str): the lens type
                            e.g. 'I' for image type I, 'II' for image type II, 'III' for image type III, 'any' for any type
            new (bool): if True, then sample new events
            size (int): the number of events to be sampled
            min_img (int): the minimum number of images
            max_img (int): the maximum number of images
            batch_size (int): the number of events to be sampled in each iteration
        Output parameters:
            param (dict): the dictionary containing the parameters of the selected events
        '''
        self.batch_size = batch_size
        if lensed==True:
            if new==False:
                param= \
                    self.selecting_n_lensed_detectable_events_from_dict(snr_threshold=snr_threshold, num_img=num_img, \
                                                                                  none_as_nan=none_as_nan, \
                                                                                  lenstype=lenstype)
            else:
                param= \
                    self.selecting_n_lensed_detectable_events_with_sampling(snr_threshold=snr_threshold, num_img=num_img, \
                                                                                  none_as_nan=none_as_nan, \
                                                                                  lenstype=lenstype, \
                                                                                  size=size, min_img=min_img, max_img=max_img)
            self.lensed_param_detectable = param
            # save the dictionary as a json file
            if jsonFile:
                file_name = './lensed_params_detectable.json'
                json_dump = json.dumps(param, cls=NumpyEncoder)
                with open(file_name, "w") as write_file:
                    json.dump(json.loads(json_dump), write_file, indent=4)

        else:
            if new==False:
                param= \
                    self.selecting_n_unlensed_detectable_events_from_dict(snr_threshold=snr_threshold)
            else:
                param= \
                    self.selecting_n_unlensed_detectable_events_with_sampling(snr_threshold=snr_threshold, size=size)
            self.lensed_param_detectable = param
            # save the dictionary as a json file
            if jsonFile:
                file_name = './unlensed_params_detectable.json'
                json_dump = json.dumps(param, cls=NumpyEncoder)
                with open(file_name, "w") as write_file:
                    json.dump(json.loads(json_dump), write_file, indent=4)

        return param

    def rate_comparision(self, size=False, snr_threshold=8., num_img=2, jsonFile=True, none_as_nan=True):
        '''
        Function to compare the detectable lensed merger rate with the unlensed merger rate
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
            unlened_rate (float): unlensed merger rate in yr^-1
            lensed_rate (float): lensed merger rate in yr^-1
            rate_ratio (float): lensed/unlensed merger rate ratio
        '''
        # calculate unlensed rate
        if self.gw_param == False:
            print('unlensed_param not sampled beforehand. Sampling now...')
            self.UnlensedCBCStatistics(nsamples=size, jsonFile=jsonFile)
        unlensed_rate = self.unlensed_rate(size=size, snr_threshold=np.max([snr_threshold]), jsonFile=jsonFile)

        # calculate lensed rate
        if self.lensed_param == False:
            print('lensed_param not sampled beforehand. Sampling now...')
            self.LensedCBCStatistics(nsamples=size, jsonFile=jsonFile)
        lensed_rate = self.lensed_rate(size=size, snr_threshold=snr_threshold, num_img=num_img, jsonFile=jsonFile)

        rate_ratio = (unlensed_rate[0]/lensed_rate[0],unlensed_rate[1]/lensed_rate[1])
        print('unlensed/lensed rate ratio = ', rate_ratio)
        return(unlensed_rate, lensed_rate, rate_ratio)

#--------------------------------------------------#

#--------------------------------------------------#
class NumpyEncoder(json.JSONEncoder):
    '''
    class for storing a numpy.ndarray or any nested-list composition as JSON file
    '''
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
#--------------------------------------------------#