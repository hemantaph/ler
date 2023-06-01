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
from ler.helperroutines import rejection_sample

class SourceGalaxyPopulationModel():
    def __init__(self, z_min=0., z_max=10., event_type='popI_II'):
        '''
        Functions to create lookup tables for redshifts and distances
        1. Redshift to co-moving distance
        2. Co-moving distance to redshift
        3. Redshift to luminosity distance
        Input parameters:
            z_min (float): minimum redshift of the source population
            z_max (float): maximum redshift of the source population
            event_type (string): popI_II,BNS,popIII,primordial,popI_II_Madau_Dickinson
        Output parameters:
            None
        '''
        self.z_min = z_min
        self.z_max = z_max
        self.create_lookup_table(z_min, z_max)

        # To find the normalization constant of the pdf p(z)
        # Define the merger-rate density function
        try:
            self.merger_rate_density = getattr(self, 'merger_rate_density_'+event_type)
        except:
            self.merger_rate_density = getattr(self, 'merger_rate_density_'+'popI_II')

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

    def sample_source_redshifts(self, size=1000, z_min=0., z_max=12.):
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
    
    def merger_rate_density_popI_II(self, zs, R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30):
        '''
        Function to compute the merger rate density (PopI/PopII)
        Input parameters:
            zs (float/array): source redshifts
            R0 (float)      : normalization constant [default: 23.9*1e-9 Mpc^-3 yr^-1]
            b2 (float)      : fitting paramters [default: 1.6]
            b3 (float)      : fitting paramters [default: 2.0]
            b4 (float)      : fitting paramters [default: 30]
        Output parameters:
            rate_density (float/array): merger rate density (Mpc^-3 yr^-1)
        '''
        rate_density = R0*(b4+1)*np.exp(b2*zs)/(b4+np.exp(b3*zs))
        return rate_density 
    
    def merger_rate_density_popI_II_Madau_Dickinson(self, zs, af=2.7,bf=5.6,cf=1.9):
        '''
        Function to compute the merger rate density (PopI/PopII)
        Parameters
        ----------
        zs : `float` 
        
        Returns
        ----------
        rate_density : `float`
        
        '''
        rate_density = (1+zs)**af / (1 + ((1+zs)/cf)**bf)
        return rate_density 
    
    def merger_rate_density_popIII(self, zs, aIII=0.66, bIII=0.3, zIII=11.6):
        '''
        Function to compute the merger rate density (PopIII)
        Parameters
        ----------
        zs : `float` 
        
        Returns
        ----------
        rate_density : `float`
        
        '''
        rate_density = np.exp(aIII*(zs-zIII))/(bIII + aIII*np.exp((aIII+bIII)*(zs-zIII)))
        return rate_density
    
    def merger_rate_density_primordial(self, zs, t0=13.786885302009708):
        '''
        Function to compute the merger rate density
        '''
        rate_density = (Planck18.age(z=zs).value/t0)**(-34/37)
        return rate_density
    
    
class CompactBinaryPopulation(SourceGalaxyPopulationModel):
    def __init__(self, z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II", \
                 model_pars=False):
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
        super().__init__(z_min=z_min, z_max=z_max, event_type = event_type)

        # self.z_min already defined in SourceGalaxyPopulationModel
        # self.z_max already defined in SourceGalaxyPopulationModel
        self.m_min      = m_min
        self.m_max      = m_max
        try:
            self.source_binary_masses = getattr(self, 'binary_masses_'+event_type)
        except:
            pass
        
        if event_type=="popI_II":
            self.model_pars = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, \
                               'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, \
                               'mu_g': 33.07, 'sigma_g': 5.69}
            if m_min<4.59:
                print('WARNING: m_min is too low for popI/II BBHs')
            if m_max>86.22:
                print('WARNING: m_max is too high for popI/II BBHs')
                
        elif event_type=="popIII":
            self.model_pars = None
            if m_min<4.59:
                print('WARNING: m_min is too low for popI/II BBHs')
            if m_max>86.22:
                print('WARNING: m_max is too high for popI/II BBHs')
                
        elif event_type=="BNS":
            self.model_pars = None
            # check the mass is for neutron stars
            if m_min<1.4:
                print('WARNING: m_min is too low for neutron stars')
            if m_max>3.0:
                print('WARNING: m_max is too high for neutron stars')
        
        elif event_type=='popI_II_Madau_Dickinson':
            self.model_pars = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, \
                               'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, \
                               'mu_g': 33.07, 'sigma_g': 5.69}
            self.source_binary_masses = getattr(self, 'binary_masses_'+'popI_II')
          
        elif event_type=='primordial':
            self.model_pars = {'Mc':30.,'sigma':0.3,'beta':1.1}

        else:
            print(f'Event Type: {event_type} is not recognised')
            
        
        # select which function to use according to event type
        #self.sample_gw_parameters = getattr(self, 'sample_'+event_type)

    def sample_gw_parameters(self, nsamples=1000, **kwargs):
        '''
        Function to sample BBH parameters from the source galaxy population model
        Input parameters:
            nsamples (int)  : number of BBHs to sample
            **kwargs        : keyword arguments
                            e.g. zs = np.array([0.1,0.2,0.3])
                                model_pars = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, \
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
        
        def warning(param_):
            try: 
                len(param_)==nsamples
            except:
                print(f'Error: len({param_}) != nsamples')
                return 0
                      
                      
        # sampling redshifts and luminosity distances
        try: kwargs['zs']
        except:
            zs = self.sample_source_redshifts(size=nsamples, z_min=z_min, z_max=z_max)
        else:
            zs = kwargs['zs']
            warning(zs)       
        luminosity_distance = self.z_to_luminosity_distance(zs) # Mpc
        
        # sampling mass1 and mass2
        try: 
            kwargs['mass_1']
            kwargs['mass_2']
        except:
            mass_1_source, mass_2_source = self.source_binary_masses(size=nsamples, model_pars=self.model_pars)
        else:
            mass_1 = kwargs['mass_1']
            warning(mass_1) 
            mass_2 = kwargs['mass_2']
            warning(mass_1)
            mass_1_source, mass_2_source = self.source_binary_masses(size=nsamples, model_pars=kwargs['model_pars'])
            
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
        geocent_time = randint.rvs(19E6, 30.6E6,size=nsamples)
        mass_1, mass_2 = mass_1_source*(1+zs), mass_2_source*(1+zs)
        
        gw_parameters = { 'mass_1': mass_1, 'mass_2': mass_2, 'mass_1_source':mass_1_source, 'mass_2_source':mass_2_source, \
                'zs':zs, 'luminosity_distance':luminosity_distance, 'iota':theta_jn, 'psi':psi, 'phase':phase, \
                'geocent_time':geocent_time, 'ra':ra, 'dec':dec  }
        
        return gw_parameters
    
    def binary_masses_popI_II(self, size, model_pars):
        """
        Function to calculate source mass1 and mass2 with PowerLaw+PEAK model 
        """

        model = p.mass_prior('BBH-powerlaw-gaussian', model_pars)
        mass_1_source, mass_2_source = model.sample(Nsample=size)
        while np.any(mass_2_source>mass_1_source):
            mass_1_source, mass_2_source = model.sample(Nsample=size)

        return(mass_1_source,mass_2_source)
    
    def binary_masses_popIII(self, size, model_pars):
        """
        Function to calculate source mass1 and mass2 with pop III origin
        """
        q = 0.9
        mass_1_source = 20.*np.ones(size)
        mass_2_source = q*mass_1_source

        return(mass_1_source,mass_2_source)
    
    def binary_masses_primordial(self, size, model_pars={'Mc':30.,'sigma':0.3,'beta':1.1}):
        """
        Function to calculate source mass1 and mass2 with PowerLaw+PEAK model 
        """

        Mc = model_pars['Mc']
        sigma = model_pars['sigma']
        beta = model_pars['beta']
        q = self.mass_ratio(size, beta)
        pdf = lambda m: np.exp(-np.log(m/Mc)**2 / (2*sigma**2)) / (np.sqrt(2*np.pi)*sigma*m)
        mass_1_source = rejection_sample(pdf, self.m_min, self.m_max, size=size)
        mass_2_source = q*mass_1_source

        return(mass_1_source,mass_2_source)
    
    def binary_masses_BNS(self, size, model_pars):
        """
        Function to calculate source mass1 and mass2 with pop III origin
        """
        q = self.mass_ratio(size=size)
        mass_1_source = np.random.uniform(1,2.5,size)
        mass_2_source = q*mass_1_source

        return(mass_1_source,mass_2_source)
    
    def mass_ratio(self, size, beta=1.1):
        
        pdf = lambda q: q**beta
        q = rejection_sample(pdf, 0, 1, size=size)
        return(q)
    


