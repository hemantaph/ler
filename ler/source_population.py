import numpy as np
import bilby
from scipy.stats import randint
#from gwcosmo import priors as p
from scipy.interpolate import interp1d
from scipy.integrate import quad
# for redshift to luminosity distance conversion
from astropy.cosmology import Planck18
import astropy.units as u
from astropy import constants as const
# for generating mass distribution
from gwcosmo import priors as p
# for multiprocessing 
# Import helper routines
from helperroutines import rejection_sample

class SourceGalaxyPopulationModel():
    def __init__(self, z_min=0.0001, z_max=10):
        '''
        Functions to create lookup tables for redshifts and distances
        1. Redshift to co-moving distance
        2. Co-moving distance to redshift
        3. Redshift to luminosity distance
        Input parameters:
            z_min (float): minimum redshift of the source population
            z_max (float): maximum redshift of the source population
        Output parameters:
            None
        '''
        self.z_min = z_min
        self.z_max = z_max
        self.create_lookup_table(z_min, z_max)

        # To find the normalization constant of the pdf p(z)
        # Define the merger-rate density function
        merger_rate_density_detector_frame = lambda z: self.merger_rate_density(z)/(1+z) 
        # Define the pdf p(z)
        pdf_unnormalized = lambda z: merger_rate_density_detector_frame(z) * self.differential_comoving_volume(z)
        # Normalize the pdf
        # this normalization factor is common no matter what you choose for z_min and z_max
        self.normalization_pdf_z = quad(pdf_unnormalized, z_min, z_max)[0]
        return None

    def create_lookup_table(self, z_min, z_max):
        '''
        Functions to create lookup tables for redshifts and comoving volume
        Input parameters:
            z_min (float): minimum redshift of the source population
            z_max (float): maximum redshift of the source population
        Output parameters:
            None
        '''
        # initialing cosmological functions for fast calculation through interpolation
        z               = np.linspace(z_min,z_max,500) # redshift
        luminosity_distance              = Planck18.luminosity_distance(z).value # luminosity distance in Mpc
        self.z_to_luminosity_distance    = interp1d( z, luminosity_distance, kind = 'cubic')

        # Create a lookup table for the differential comoving volume
        dVcdz = Planck18.differential_comoving_volume(z).value*4*np.pi
        self.differential_comoving_volume = interp1d(z, dVcdz, kind='linear', fill_value='extrapolate')
        return None

    def sample_source_redshifts(self, size=1000, z_min=0, z_max=12):
        '''
        Function to sample source redshifts from the source galaxy population model
        Input parameters:
            size (int): number of source redshifts to sample
            z_min (float): minimum redshift of the source population
            z_max (float): maximum redshift of the source population
        Output parameters:
            zs (array): array of source redshifts
        '''
        # Define the merger-rate density function
        merger_rate_density_detector_frame = lambda z: self.merger_rate_density(z)/(1+z) 
        # Define the pdf p(z)
        pdf_unnormalized = lambda z: merger_rate_density_detector_frame(z) * self.differential_comoving_volume(z)
        # Normalize the pdf
        normalization = self.normalization_pdf_z
        pdf = lambda z: pdf_unnormalized(z) / normalization
        # Sample the redshifts using rejection sampling
        zs = rejection_sample(pdf, z_min, z_max, size=size)
        return zs
    
    def merger_rate_density(self, zs, R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30):
        '''
        Function to compute the merger rate density
        Input parameters:
            zs (float/array): source redshifts
            R0 (float)      : normalization constant [default: 23.9*1e-9]
            b2 (float)      : fitting paramters [default: 1.6]
            b3 (float)      : fitting paramters [default: 2.0]
            b4 (float)      : fitting paramters [default: 30]
        Output parameters:
            rate_density (float/array): merger rate density
        '''
        rate_density = R0*(b4+1)*np.exp(b2*zs)/(b4+np.exp(b3*zs))
        return rate_density 

class CompactBinaryPopulation(SourceGalaxyPopulationModel):
    def __init__(self, z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "StellarBBH"):
        '''
        Daughter Class to generate a population of compact binaries
        Input parameters:
            m_min (float): minimum mass of the BBHs
            m_max (float): maximum mass of the BBHs
            event_type (string): type of event to generate
                                e.g. "StellarBBH", "BNS"
        Output parameters:
            None
        '''
        # initialized SourceGalaxyPopulationModel mother class
        super().__init__(z_min=z_min, z_max=z_max)

        # self.z_min already defined in SourceGalaxyPopulationModel
        # self.z_max already defined in SourceGalaxyPopulationModel
        self.m_min      = m_min
        self.m_max      = m_max

        if event_type=="BNS":
            # check the mass is for neutron stars
            if m_min<1.4:
                print('WARNING: m_min is too low for neutron stars')
            if m_max>3.0:
                print('WARNING: m_max is too high for neutron stars')
        if event_type=="StellarBBH":
            # check the mass is for neutron stars
            if m_min<4.59:
                print('WARNING: m_min is too low for stellar Black Holes')
            if m_max>86.22:
                print('WARNING: m_max is too high for stellar Black Holes')
        
        # select which function to use according to event type
        self.sample_gw_parameters = getattr(self, 'sample_'+event_type)

    def sample_StellarBBH(self, nsamples=1000, **kwargs):
        '''
        Function to sample BBH parameters from the source galaxy population model
        Input parameters:
            nsamples (int)  : number of BBHs to sample
            **kwargs        : keyword arguments
                            e.g. zs = np.array([0.1,0.2,0.3])
                                model_pars_gwcosmo = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, \
                                   'mmin': m_min, 'mmax': m_max, 'lambda_peak': 0.08, \
                                   'mu_g': 33.07, 'sigma_g': 5.69}
        Output parameters:
            gw_parameters (dict): dictionary of GW parameters
                                e.g. gw_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', \
                                    'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec']
        '''
        z_min = self.z_min
        z_max = self.z_max
        m_min = self.m_min
        m_max = self.m_max
        
        # sampling redshifts and luminosity distances
        try: kwargs['zs']
        except:
            zs = self.sample_source_redshifts(size=nsamples, z_min=z_min, z_max=z_max)
        else:
            zs = kwargs['zs']
            try: len(zs)==nsamples
            except:
                print('WARNING: len(zs) != nsamples. Setting nsamples = len(zs)')
                nsamples = len(zs)
        luminosity_distance = self.z_to_luminosity_distance(zs) # Mpc
        
        # sampling mass_1 from PowerLaw+PEAK model
        # source frame mass samples
        try: kwargs['model_pars_gwcosmo']
        except:
            # default model parameters
            model_pars_gwcosmo = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, \
                                   'mmin': m_min, 'mmax': m_max, 'lambda_peak': 0.08, \
                                   'mu_g': 33.07, 'sigma_g': 5.69}
        else:
            model_pars_gwcosmo = kwargs['model_pars_gwcosmo']

        model_pars_gwcosmo['mmin'] = m_min
        model_pars_gwcosmo['mmax'] = m_max
        self.model_pars_gwcosmo = model_pars_gwcosmo
        model = p.mass_prior('BBH-powerlaw-gaussian', model_pars_gwcosmo)
        mass_1_source, mass_2_source = model.sample(Nsample=nsamples)
        while np.any(mass_2_source>mass_1_source):
            mass_1_source, mass_2_source = model.sample(Nsample=nsamples)
            
        # Sample all other parameters
        # use bilby priors
        bilby.core.utils.logger.disabled = True
        prior_default = bilby.gw.prior.BBHPriorDict()
        # draw associated angles 
        ra = prior_default['ra'].sample(nsamples) 
        dec = prior_default['dec'].sample(nsamples) 
        psi = prior_default['psi'].sample(nsamples) 
        theta_jn = prior_default['theta_jn'].sample(nsamples) 
        phase = prior_default['phase'].sample(nsamples) 
        # compute GPS time
        geocent_time = randint.rvs(20E6, 31.6E6,size=nsamples)
        mass_1, mass_2 = mass_1_source*(1+zs), mass_2_source*(1+zs)
        
        gw_parameters = { 'mass_1': mass_1, 'mass_2': mass_2, 'mass_1_source':mass_1_source, 'mass_2_source':mass_2_source, \
                'zs':zs, 'luminosity_distance':luminosity_distance, 'iota':theta_jn, 'psi':psi, 'phase':phase, \
                'geocent_time':geocent_time, 'ra':ra, 'dec':dec  }
        
        return gw_parameters
    
    def sample_BNS(self, nsamples=1000, **kwargs):
        '''
        Function to sample BNS parameters from the source galaxy population model
        Input parameters:
            nsamples (int)  : number of BBHs to sample
            **kwargs        : keyword arguments
                            e.g. zs = np.array([0.1,0.2,0.3])
        Output parameters:
            gw_parameters (dict): dictionary of GW parameters
                                e.g. gw_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', \
                                    'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec']
        '''
        z_min = self.z_min
        z_max = self.z_max
        m_min = self.m_min
        m_max = self.m_max
            
        # sampling redshifts and luminosity distances
        try: kwargs['zs']
        except:
            zs = self.sample_source_redshifts(size=nsamples, z_min=z_min, z_max=z_max)
        else:
            zs = kwargs['zs']
            try: len(zs)==nsamples
            except:
                print('WARNING: len(zs) != nsamples. Setting nsamples = len(zs)')
                nsamples = len(zs)
        luminosity_distance = self.z_to_luminosity_distance(zs) # Mpc
        
        # sampling mass_1 from PowerLaw+PEAK model
        # source frame mass samples
        mass_1_source = np.random.uniform(m_min,m_max,size=nsamples)
        mass_ratio = np.random.uniform(0.1,1.,size=nsamples)
        mass_2_source = mass_1_source*mass_ratio
            
        # Sample all other parameters
        bilby.core.utils.logger.disabled = True
        prior_default = bilby.gw.prior.BBHPriorDict()
        # draw associated angles 
        ra = prior_default['ra'].sample(nsamples) 
        dec = prior_default['dec'].sample(nsamples) 
        psi = prior_default['psi'].sample(nsamples) 
        theta_jn = prior_default['theta_jn'].sample(nsamples) 
        phase = prior_default['phase'].sample(nsamples) 
        # compute GPS time
        geocent_time = randint.rvs(20E6, 31.6E6,size=nsamples)
        mass_1, mass_2 = mass_1_source*(1+zs), mass_2_source*(1+zs)
        
        return { 'mass_1': mass_1, 'mass_2': mass_2, 'mass_1_source':mass_1_source, 'mass_2_source':mass_2_source, \
                'zs':zs, 'luminosity_distance':luminosity_distance, 'iota':theta_jn, 'psi':psi, 'phase':phase, \
                'geocent_time':geocent_time, 'ra':ra, 'dec':dec  }

