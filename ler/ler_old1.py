import numpy as np
import bilby

# for fast SNR calculatiopn call quintet
#from TaylorF2_with_fLim import quintet as Quintet
from quintet import Quintet
#from quintet import Quintet

from scipy.stats import uniform, randint, gengamma, rayleigh, norm
#from gwcosmo import priors as p
from scipy.interpolate import interp1d
from scipy.integrate import quad

from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from lenstronomy.Util.param_util import phi_q2_ellipticity

# for redshift to luminosity distance conversion
from astropy.cosmology import Planck18
from astropy.cosmology import z_at_value
import astropy.units as u

# for generating mass distribution
from gwcosmo import priors as p

# For outputting results
import pandas

# for multiprocessing 
import multiprocessing 
from multiprocessing import Pool
#import task

# Conversions from SI units to CGS units
C = 299792458.
G = 6.67408*1e-11
Mo = 1.989*1e30
Mpc = 3.086*1e22

# Define the maximum number of images
n_max_images = 5

'''
------------------------------------------------
    Grandmother class containing common methods 
------------------------------------------------
'''
class LeR():
    ####################################################
    #                                                  #
    #             Class initialization                 #
    #                                                  #
    ####################################################
    def __init__(self, npool=int(4), n_for_quintet=200, z_min=0.0001, z_max=10., m_min=4.59, m_max=86.22, event_type = "BBH", det_sensitivity='design',equal_mass=False):
        '''
        -----------------
        input parameters
        -----------------
        n_for_quintet   : number of bilby SNR use by quintet for interpolation
                          higher the number, more is the accuracy
                          no performance gain beyond n_for_quintet=1000
        z_min           : minimum value of redshift for BBH/BNS
                          example -> z_min=0.0001= 0.44312066Mpc =1.44*1e6ly
        z_max           : maximum value of redshift for BBH/BNS
                          example -> z_min=10.= 105999.14Mpc = 3.46*1e11ly
        m_min           : minimum mass of compact object for the chosen type of binaries (BBH or BNS)
        m_max           : maximum mass of compact object for the chosen type of binaries (BBH or BNS)
        event_type      : "BBH" (Binary black hole), "BNS" (Binary neutron star)
        det_sensitivity : sensitivity of the GW detectors
                          examples: "O1","O2","O3"
        '''
        self.z_min      = z_min
        self.z_max      = z_max
        self.m_min      = m_min
        self.m_max      = m_max
        self.event_type = event_type
        self.gw_param   = False          # this needed not repeat sampling 
        self.npool      = npool
        # initializing function for fast SNR calculation
        self.equal_mass = equal_mass
        
        # initializing function for fast SNR calculation
        self.quin_      = Quintet(mtot_min=m_min*2, mtot_max=m_max*2, nsamples=n_for_quintet, sensitivity=det_sensitivity)
        
        # initialing cosmological functions for fast calculation through interpolation
        z               = np.linspace(0,z_max,500) # red-shift
        Dc              = Planck18.comoving_distance(z).value # co-moving distance in Mpc
        luminosity_distance              = Planck18.luminosity_distance(z).value # luminosity distance in Mpc
        self.z_to_Dc    = interp1d( z, Dc, kind = 'cubic')
        self.Dc_to_z    = interp1d( Dc, z, kind = 'cubic')
        self.z_to_luminosity_distance    = interp1d( z, luminosity_distance, kind = 'cubic')
        self.gw_param   = False 
        # for angular diameter distance
        quad_ = [] # refer to quad integradtion from scipy 
        for ii in range(len(z)):
            quad_.append(quad(Planck18._inv_efunc_scalar, 0., z[ii], args=Planck18._inv_efunc_scalar_args)[0])
        quad_ = np.array(quad_)
        self.quad_int = interp1d(z, np.array(quad_), kind = 'cubic')

    ####################################################
    #                                                  #
    #   To sample gravitational wave parameters       #
    #                                                  #
    ####################################################
    def sample_gw_parameters(self, zs=False, nsamples=1000, m_min=4.59, m_max=86.22, z_max=False, event_type = "BBH",
                  model_pars_gwcosmo = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, \
                                       'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, \
                                        'mu_g': 33.07, 'sigma_g': 5.69}, ):
        '''
        -----------------
        Input parameters
        -----------------
        nsamples           : number of samples. nsamples>1e6 recomended for rate estimation.
        m_min              : minimum mass of compact object for the chosen type of binaries (BBH or BNS)
        m_max              : maximum mass of compact object for the chosen type of binaries (BBH or BNS)
        z_max               : maximum value of source redshift (BBH/BNS)
                             example -> z_min=10.= 105999.14Mpc = 3.46*1e11ly
        event_type         : "BBH" (Binary black hole), "BNS" (Binary neutron star)
        model_pars_gwcosmo : contains fitting parameters for PowerLaw+PEAK model for mass_1 and q samples.
        -----------------
        Return values
        -----------------
        gw_params_         : Dict containing array of all parameters.
        
        '''
        #np.random.seed(123) # for reproduceability
        # sampling redshift from Uniform distribution between z_min and z_max
        if z_max==False:
            z_max = self.z_max
        z_min = self.z_min
        
        if np.array(zs).any()==False:
            zs = np.random.uniform(z_min,z_max,size=nsamples)
        else :
            nsamples = len(zs)
            
        luminosity_distance = self.z_to_luminosity_distance(zs) #Mpc    

        # source frame mass samples
        if event_type=="BBH":
            model_pars_gwcosmo['mmin'] = m_min
            model_pars_gwcosmo['mmax'] = m_max
            model=p.mass_prior('BBH-powerlaw-gaussian', model_pars_gwcosmo)
            mass_1, mass_2 = model.sample(Nsample=nsamples)
            while np.any(mass_2>mass_1):
                mass_1, mass_2 = model.sample(Nsample=nsamples)
        elif event_type=="BNS":
            mass_1 = np.random.uniform(m_min,m_max,size=nsamples)
            mass_ratio = np.random.uniform(0.1,1.,size=nsamples)
            mass_2 = mass_1*mass_ratio
        else:
            raise Exception("enter either BBH or BNS for event_type")
        if self.equal_mass:
            mass_2=mass_1
        
        # Sample all other parameters
        prior_default = bilby.gw.prior.BBHPriorDict()
        # draw associated angles 
        ra = prior_default['ra'].sample(nsamples) 
        dec = prior_default['dec'].sample(nsamples) 
        psi = prior_default['psi'].sample(nsamples) 
        theta_jn = prior_default['theta_jn'].sample(nsamples) 
        phase = prior_default['phase'].sample(nsamples) 
        # compute GPS time
        GPStimeValue = randint.rvs(0, 31.6E6,size=nsamples)
        
        # return dict containing array of all parameters
        return({ 'mass_1':mass_1, 'mass_2':mass_2, 'zs':zs, 'luminosity_distance':luminosity_distance, 'iota':theta_jn, 'psi':psi, 'phase':phase, \
                      'geocent_time':GPStimeValue, 'ra':ra, 'dec':dec  })
    
    ########################################################################
    #                                                                      #
    #     merger rate distribution wrt to redshift and comoving volume     #
    #                                                                      #
    ########################################################################
    # source frame merger rate at source redshift z
    def merger_rate(self, zs, R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30):
        '''
        -----------------
        Input parameters
        -----------------
        zs               : redshift of source
        R0               : local merger density rate
        b2, b3, b4       : merger rate density fitting parameters.
        -----------------
        Return values
        -----------------
        merger rate at zs : 
        '''
        # differential_comoving_volume at reshift zs
        dVc_dz = 4*np.pi*Planck18.differential_comoving_volume(zs).value
        
        rate_density = R0*(b4+1)*np.exp(b2*zs)/(b4+np.exp(b3*zs))
        rate = rate_density/(1+zs) * dVc_dz
        
        return(rate)

'''
------------------------------------------------
     father class for unlensed events      
------------------------------------------------
'''
class UnlensedCBCStatistics(LeR):
    def __init__(self, npool=int(4), n_for_quintet=200, z_min=0.0001, z_max=10., m_min=4.59, m_max=86.22, event_type = "BBH", det_sensitivity='design',equal_mass=False):
        '''
        -----------------
        input parameters
        -----------------
        n_for_quintet   : number of bilby SNR use by quintet for interpolation
                          higher the number, more is the accuracy
                          no performance gain beyond n_for_quintet=1000
        z_min           : minimum value of redshift for BBH/BNS
                          example -> z_min=0.0001= 0.44312066Mpc =1.44*1e6ly
        z_max           : maximum value of redshift for BBH/BNS
                          example -> z_min=10.= 105999.14Mpc = 3.46*1e11ly
        m_min           : minimum mass of compact object for the chosen type of binaries (BBH or BNS)
        m_max           : maximum mass of compact object for the chosen type of binaries (BBH or BNS)
        event_type      : "BBH" (Binary black hole), "BNS" (Binary neutron star)
        det_sensitivity : sensitivity of the GW detectors
                          examples: "O1","O2","O3"
        '''  
        super().__init__(npool=int(4), n_for_quintet=200, z_min=0.0001, z_max=10., m_min=4.59, m_max=86.22, event_type = "BBH", det_sensitivity='design',equal_mass=equal_mass)
        
    ########################################################################
    #                                                                      #
    #             detectable unlensed merger rate yr^-1                    #
    #                                                                      #
    ########################################################################
    def rate_detectable(self, zs=False, size=1000, model_pars_gwcosmo = False,):
        '''
        -----------------
        Input parameters
        -----------------
        zs               : redshift of source
        size             : number of samples
        model_pars_gwcosmo : dict containing fitting parameters for PowerLaw+PEAK model for mass_1 and q samples.
        -----------------
        Return values
        -----------------
        rate_pdet    : detectable unlensed merger rate per year with pdet
        rate_step    : detectable unlensed merger rate per year with step function
        '''
        z_max = self.z_max
        z_min = self.z_min
        m_min = self.m_min
        m_max = self.m_max
        if model_pars_gwcosmo==False:
            model_pars_gwcosmo = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, \
                                   'mmin': m_min, 'mmax': m_max, 'lambda_peak': 0.08, \
                                    'mu_g': 33.07, 'sigma_g': 5.69}
            
            
        if np.array(zs).any()==False:
            zs = np.random.uniform(self.z_min, self.z_max, size = size)
            gw_params_ = self.sample_gw_parameters(zs=zs, m_min=m_min, m_max=m_max, z_max=z_max, \
                                    event_type = self.event_type, model_pars_gwcosmo=model_pars_gwcosmo)
            self.gw_param = gw_params_
            #print('a')
        # check whether the gw parameters are already sampled or not
        elif self.gw_param!=False:
            zs   = self.gw_param['zs']
            size = len(zs)
            gw_params_ = self.gw_param
            #print('b')
        else:
            gw_params_ = self.sample_gw_parameters(zs=zs, m_min=m_min, m_max=m_max, z_max=z_max, \
                                    event_type = self.event_type, model_pars_gwcosmo=model_pars_gwcosmo)
            size = len(zs)
            self.gw_param  = gw_params_
            #print('c')

        # source frame
        mass_1 = gw_params_['mass_1']
        mass_2 = gw_params_['mass_2']
        z = gw_params_['zs']
        # detector frame
        mass_1 = mass_1*(1+z)
        mass_2 = mass_2*(1+z)
        luminosity_distance = gw_params_['luminosity_distance'] #Mpc
        theta_jn,ra, dec,psi,phase,geocent_time = gw_params_['iota'], gw_params_['ra'], gw_params_['dec'], \
                                                    gw_params_['psi'], gw_params_['phase'], gw_params_['geocent_time'] 

        # calculate probability of detection for each gw parameters set
        param_ = {'mass_1':mass_1,'mass_2':mass_2,'luminosity_distance':luminosity_distance,'iota':theta_jn,'psi':psi,'phase':phase,'ra':ra,'dec':dec, 'geocent_time':geocent_time}
        pdet_ = self.quin_.pdet(param_)['pdet_net']
        
        # montecarlo integration
        #---------------------------------------------------#
        #                   with pdet
        #---------------------------------------------------#
        merger_rate = self.merger_rate(z)
        rate_pdet = np.sum(merger_rate*pdet_)*(z_max-z_min)/float(size)
        #---------------------------------------------------#
        
        #---------------------------------------------------#
        #                with step function
        #---------------------------------------------------#
        step_            = np.zeros(len(pdet_),dtype=float)
        step_[pdet_>0.5] = 1.
        rate_step        = np.sum(merger_rate*step_)*(z_max-z_min)/float(size)
        #---------------------------------------------------#
        
        return(rate_pdet,rate_step)
    
'''
------------------------------------------------
     mother class for lensed events      
------------------------------------------------
'''
class LensedCBCStatistics(LeR):
    def __init__(self, npool=int(4), n_for_quintet=200, z_min=0.0001, z_max=10., m_min=4.59, m_max=86.22, event_type = "BBH", det_sensitivity='design',equal_mass=False):
        '''
        -----------------
        input parameters
        -----------------
        n_for_quintet   : number of bilby SNR use by quintet for interpolation
                          higher the number, more is the accuracy
                          no performance gain beyond n_for_quintet=1000
        z_min           : minimum value of redshift for BBH/BNS
                          example -> z_min=0.0001= 0.44312066Mpc =1.44*1e6ly
        z_max           : maximum value of redshift for BBH/BNS
                          example -> z_min=10.= 105999.14Mpc = 3.46*1e11ly
        m_min           : minimum mass of compact object for the chosen type of binaries (BBH or BNS)
        m_max           : maximum mass of compact object for the chosen type of binaries (BBH or BNS)
        event_type      : "BBH" (Binary black hole), "BNS" (Binary neutron star)
        det_sensitivity : sensitivity of the GW detectors
                          examples: "O1","O2","O3"
        '''
        super().__init__(npool=int(4), n_for_quintet=200, z_min=0.0001, z_max=10., m_min=4.59, m_max=86.22, event_type = "BBH", det_sensitivity='design',equal_mass=equal_mass)
        
        # create a lookup table for the lens redshift draws
        r = np.linspace(0, 1, num = 100)
        u = 10*r**3 - 15*r**4 + 6*r**5
        self.gausshyperPPF = interp1d(u, r, kind = 'cubic') # Computes r(u)
        
    ########################################################################
    #                                                                      #
    #                optical depth at source redshift                      #
    #                                                                      #
    ########################################################################
    def optical_depth(self, zs):
        ''' Strong lensing optical depth. From Haris et al. 2018
        -----------------
        Input parameters
        -----------------
        zs : redshift of source
        -----------------
        Return values
        -----------------
        optical depth  : measure at redshift z
        '''
        # z to luminosity_distance (luminosity_distance) conversion
        Dc = self.z_to_Dc(zs)*1e-3  # 1e-3 converts Mpc to Gpc
        return( (Dc/62.2)**3 )

    ########################################################################
    #                                                                      #
    #                     Lens parameters sampler                          #
    #                                                                      #
    ########################################################################
    def PEMDGalaxy_param_sampler(self, zs, ndraw=20):
        ''' Sample lens parameters from the PEMD Galaxy catalog.
        -----------------
        Input parameters
        -----------------
        zs     : redshift of source (numpy array random sample, ascending order prefered)
        ndraw  : number of draws within which an event is checked for lensing condition wrt Einstein radius
        
        -----------------
        Return values
        -----------------
        zl     : redshift of lens 
        theta_E : Einstein radius 
        q      : axis ratio of the galaxy
        phi    : axis rotation of the galaxy
        gamma  : spectral index of density profile
        '''
        # making sure zs is numpy array ang its in ascending order
        zs = np.sort(np.array(zs)) # 1D array   
        # Hubble distance:
        H0d = Planck18._hubble_distance.value # unit: Mpc
        
        # defining useful functions
        angular_diameter_distance = lambda zs0 : H0d * (self.quad_int(zs0)) / (zs0 + 1.0)
        angular_diameter_distance_z1z2 = lambda zl0, zs0 : H0d * (self.quad_int(zs0)-self.quad_int(zl0)) / (zs0 + 1.0)
        
        # empty list for appending the results
        zs_final, zl_final, theta_E_final, a_final = [],[],[],[]
        arrLen = len(zs)
        len_ = arrLen
        
        while(True):
            # to make numpy array operation easier
            # create a matrix of shape(len(zs),ndraw)
            # each elements in a row has same value of source redshift
            zs_ = np.transpose(np.array([zs,]))*np.ones((1,ndraw))
        
            # Draw the lens redshift                 
            r = self.gausshyperPPF(np.random.uniform(0, 1, size = (len_,ndraw)))
            DcL = self.z_to_Dc(zs_)*r # corresponding element-wise multiplication between 2 arrays 
            zl = self.Dc_to_z(DcL) # 2D array

            # Draw the velocity dispersion
            a = gengamma.rvs(2.32/2.67, 2.67, size = (len_,ndraw))
            sigma = 161.*a

            # angular diameter distances
            Ds = angular_diameter_distance(zs_) # 2D array with shape=(len(zs),ndraw)
            luminosity_distances = angular_diameter_distance_z1z2(zl, zs_) # shape(len(zs),ndraw)

            # Einstein radius 
            theta_E = 4*np.pi*sigma**2*luminosity_distances/((2.998E5)**2*Ds)*648000./np.pi # unit: arcsec, speed of light: 2.998E5 kms^-1 

            # Rejection sample with array of size (len_,ndraw)
            u = np.random.uniform(0, 3**2, size = (len_,ndraw))
            X = u < theta_E**2 # boolean array
            # first, 'True' value position is found for each row in X (each row correspond to one value of zs)
            # a row can contain no 'True' value or one or more 
            # but if there is no 'True' value in a row, the corresponding zs of that row will be recycle and \ 
            # feed it to next while loop
            idx = np.argmax(X,axis=1) # position (column) of the highest value in the boolean array. 
            idx = np.array([np.arange(0,len_,1),idx]).T # combining row and column position to form a 2D array
            # choose only those row that has satisfy the condition u < theta_E**2
            idx2 = X[idx[:,0],idx[:,1]] == True # 1st column-elements=idx[:,0] , 2nd column-elements=idx[:,1]
            # update idx, to include index for only that satisfy the condition
            idx4 = idx[idx2] 
            # now find the index for those that dont satisfy the condition
            idx5 = np.invert(idx2) # inverting the boolean value
            
            # 1D array of the results
            zs_final.extend(zs[idx4[:,0]].flatten().tolist())
            zl_final.extend(zl[idx4[:,0],idx4[:,1]].flatten().tolist())
            theta_E_final.extend(theta_E[idx4[:,0],idx4[:,1]].flatten().tolist())
            # this parameter 'a' will be used in calculation of axis ratio of the galaxy
            a_final.extend(a[idx4[:,0],idx4[:,1]].flatten().tolist())

            # zs elements that dont satisfy the condition forms new zs 1D array for recycling
            zs = zs[idx5]
            len_ = len(zs)
            if len_==0:
                break # if there is no zs element left to be recycled
        
        # sorting the results according to the input zs 
        matrix = np.array([zs_final,zl_final,theta_E_final,a_final]).T
        matrix = matrix[matrix[:,0].argsort()]
        # delete unncecessary variables to save RAM memory
        # this matters when len(zs)>1e5
        del zs_final,zl_final,theta_E_final,a_final,zs_,r,DcL,zl,a,sigma,Ds,luminosity_distances,theta_E,u,X,idx,idx2,idx4,idx5
        
        #################################################################
        # Sample axis ratio of the galaxy
        
        # 'a' was sampled from gamma distribution
        a_final = matrix[:,3] # 1D array of len = len(zs)
        len_ = arrLen
        idx2 = np.arange(0,arrLen,1) # helps in recycling the rejected samples
        b = np.ones(arrLen)
        # this while loop is simpler as it deals with only 1D array
        while(True):
            s = 0.38 - 0.09177*a_final[idx2]
            # sampled from rayleigh distribution
            b[idx2] = rayleigh.rvs(scale = s)
            # now recycle those index that dont satisfy the condition
            idx2 = np.array(np.where(b >= 0.8)).flatten()
            #print('rejected index array len=',len(idx2))
            if len(idx2) ==0:
                break
        # axis-ratio with size = len(zs)     
        q = 1-b
            
        #################################################################
        # Draw an angle for the axis rotation of the galaxy
        psi = uniform.rvs(size=arrLen, scale = 2*np.pi)
        # Draw a density profile
        gamma = norm.rvs(size=arrLen, loc = 2, scale = 0.2)
        
        return {'zs':matrix[:,0], 'zl':matrix[:,1], 'theta_E':matrix[:,2], 'axis-ratio':q, 'psi':psi, 'gamma':gamma }
        #return np.array(zs_final), np.array(zl_final), np.array(theta_E_final), q, phi, gamma
    
    ########################################################################
    #                                                                      #
    #                      Lens shear sampler                              #
    #                                                                      #
    ########################################################################
    def sample_galaxy_shear(self, size):
        ''' Sample galaxy shears.
        -----------------
        Input parameters
        -----------------
        size  :  size of each sample array
        -----------------
        Return values
        -----------------
        gamma_1, gamma_2 : external shear in two direction
        '''
        # Draw an external shear
        gamma_1 = norm.rvs(size=size, scale = 0.05)
        gamma_2 = norm.rvs(size=size, scale = 0.05)
        return gamma_1, gamma_2
    
    ########################################################################
    #                                                                      #
    #                     lensed rate calculation                          #
    #                                                                      #
    ########################################################################
    def lensed_rate_detectable(self, zs=False,size=1000,ndraws=20):
        ''' Calculate the lensed rate.
        -----------------
        Input parameters
        -----------------
        zs     : redshift of source binary (1d array(float)).  condition: z_min<elements<z_max
                 if no array is given, zs is generated between z_min and z_max with size=size
        npool  : number CPU cores to use for solving lens equations      
        size   : size of zs array
                 if zs array is provided it will be replace with size=len(size)
        ndraws : extra iteration needed for generating lens characteristics
        -----------------
        Return values
        -----------------
        rate_pdet    : detectable lensing event rate with pdet
        rate_step    : detectable lensing event rate with step function
        '''
        # will also be use for unlensed event calculation
        if np.array(zs).any()==False:
            zs = np.random.uniform(self.z_min, self.z_max, size = size)
            self.gw_param = self.sample_gw_parameters(zs=zs)
            #print('a')
        # check whether the gw parameters are already sampled or not
        elif self.gw_param!=False:
            zs   = self.gw_param['zs']
            size = len(zs)
            #print('b')
        else:
            self.gw_param = self.sample_gw_parameters(zs=zs)
            zs            = self.gw_param['zs']
            size          = len(zs)
            #print('c')
        
        # get param corresponding to strongly lensed event
        paramL = self.create_lensed_images(zs=True,npool=self.npool,size=size,ndraws=ndraws)
        self.lens_param = paramL
        zs     = paramL['zs']    # check this line again
        pdet   = paramL['Pdet'] 
        
        # montecarlo integration
        #---------------------------------------------------#
        #                   with pdet
        #---------------------------------------------------#
        rate_pdet = (1/float(size))*(self.z_max-self.z_min)*np.sum(self.merger_rate(zs)*self.optical_depth(zs)*pdet)
        #---------------------------------------------------#
        
        #---------------------------------------------------#
        #                with step function
        #---------------------------------------------------#
        step_           = np.zeros(len(pdet),dtype=float)
        step_[pdet>0.5] = 1.
        rate_step       = (1/float(size))*(self.z_max-self.z_min)*np.sum(self.merger_rate(zs)*self.optical_depth(zs)*step_)
        #---------------------------------------------------#
        
        return(rate_pdet,rate_step)
    
    ########################################################################
    #                                                                      #
    #                     lensed event parameters                          #
    #                                                                      #
    ########################################################################
    def create_lensed_images(self, zs=False, npool=int(4),size=1000,ndraws=20):
        '''
        to find lens parameters 
        -----------------
        Input parameters
        -----------------
        zs     : redshift of source binary (1d array(float)).  condition: z_min<elements<z_max
                 if no array is given, zs is generated between z_min and z_max with size=size
        npool  : number CPU cores to use for solving lens equations      
        size   : size of zs array
                 if zs array is provided it will be replace with size=len(size)
        ndraws : extra iteration needed for generating lens characteristics
        -----------------
        Return values
        -----------------
        dict    : dictionary lensed event parameters. each keys has corresponding 1d array(float) or array(object) with different row len 
                  {'image_type':imageType_, 'magnifications':magnifications, 'time_delays':time_delays,\
                                 'Pdet':pdet2, 'theta_E':einstein_radius, 'zl': zl, 'zs':zs}

        -----------------
        Example
        -----------------
        >>> lensed_event = create_lensed_images()
        >>> print("Image type, magnification, time delay, pdet, theta_E, zl, zs:", lensed_event['image_type'], lensed_event['magnifications'], lensed_event['time_delays'], lensed_event['Pdet'], lensed_event['theta_E'], lensed_event['zl'], lensed_event['zs'])
        '''
        # Get the source redshifts
        # will also be use for unlensed event calculation
        if np.array(zs).any()==False:
            zs = np.random.uniform(self.z_min, self.z_max, size = size)
            self.gw_param = self.sample_gw_parameters(zs=zs)
        # check whether the gw parameters are already sampled or not
        elif self.gw_param!=False:
            zs   = self.gw_param['zs']
            size = len(zs)
        else:
            self.gw_param = self.sample_gw_parameters(zs=zs)
            zs            = self.gw_param['zs']
            size          = len(zs)
        # Get the GW parameters
        gw_param = self.gw_param
        mass_1, mass_2, luminosity_distance, iota, psi, phase, geocent_time, ra, dec, zs = gw_param['mass_1'], gw_param['mass_2'], gw_param['luminosity_distance'], gw_param['iota'], gw_param['psi'], gw_param['phase'], gw_param['geocent_time'], gw_param['ra'], gw_param['dec'], gw_param['zs']

        
        #---------------------------------------------------#
        #          sampling the lens characteristics
        #---------------------------------------------------#
        # get the lens' parameters --> {'zs', 'zl', 'theta_E', 'axis-ratio', 'phi', 'gamma'}
        lens_parameters = self.PEMDGalaxy_param_sampler(zs,ndraw=ndraws) # Note: Includes also the source redshift
        #---------------------------------------------------#
        # external shear params to the 'PEMD' galaxy lens
        gamma1, gamma2 = self.sample_galaxy_shear(size)
        # ellipticity of the galaxy lens
        e1, e2         = phi_q2_ellipticity(lens_parameters['psi'], lens_parameters['axis-ratio'])
        lensModelList  = np.array(['EPL_NUMBA', 'SHEAR'])*np.ones((size,2),dtype=object) # Create the lens model list (note: can be a different lens model for different samples)
        zl          = lens_parameters['zl']
        zs        = lens_parameters['zs']
        gamma = lens_parameters['gamma']
        # free up some memory
        del lens_parameters['psi'],lens_parameters['axis-ratio'],lens_parameters['zl'],lens_parameters['zs']
        
        # Draw random positions relative to the galaxy lens, on the image plane
        etamax = 1.5
        eta      = etamax*uniform.rvs(size=size)**0.5
        phi      = uniform.rvs(size=size, scale = 2*np.pi)            
        x_source = np.array([eta*np.cos(phi), eta*np.sin(phi)])
        
        #---------------------------------------------------#
        #                  Multiprocessing 1
        #---------------------------------------------------#
        # FIXME: I think this bit needs to be clarified. It's a bit difficult to parse here what's going on. 
        # now get it in custic curve in the source palne # FIXME: "It" is not clear here. What is this?
        iterations = np.arange(size)
        input_arguments = np.array([x_source[0],x_source[1],e1,e2,gamma,gamma1,gamma2, zl,zs, iterations], dtype=object).T
        input_arguments = np.concatenate((input_arguments,lensModelList),axis=1)
        # Initialize the image positions and lens argument list here.
        x0_image_positions = np.ones((size, n_max_images))*np.nan; x1_image_positions = np.ones((size, n_max_images))*np.nan; kwargs_lenses = np.ones(size, dtype=object); n_images = np.ones(size, dtype=int)
        # Solve the lens equation in parallel
        with Pool(processes=npool) as pool:
            # call the same function with different data in parallel
            # imap->retain order in the list, while map->doesn't
            for result in pool.map(self.solve_lens_equation,input_arguments): # NOTE: self.solve_lens_equation returns (x0_image_position,x1_image_position,n_images, kwargs_lens, zs)
                # report the value to show progress
                # return(x0_image_position,x1_image_position,len(x0_image_position), kwargs_lens, zs)
                # We assume that there are maximum of n_max_images images. If there are less, we fill the rest with nans. If this choice is confusing, know that this makes it a LOT easier to handle the data later on.
                iter = result[5]
                n_image = result[2]
                n_images[iter] = n_image
                x0_image_position = np.ones(n_max_images)*np.nan
                x1_image_position = np.ones(n_max_images)*np.nan
                x0_image_position[:n_image] = result[0]
                x1_image_position[:n_image] = result[1]
                x0_image_positions[iter] = x0_image_position
                x1_image_positions[iter] = x1_image_position
                kwargs_lenses[iter] = result[3]
        #---------------------------------------------------#
 
        # store only strong lensed values
        # lensed_image_positions[:,2] gives number of image for that lensed event
        pick_strongly_lensed = np.where(n_images>1)[0]
        
        # getting values for only strongly lensed condition 
        # this is done to solve lens equation only for those events that satisfy the condition  
        einstein_radius = lens_parameters['theta_E'][pick_strongly_lensed]
        del lens_parameters
        # mass_1, mass_2, luminosity_distance, iota, psi, phase, geocent_time, ra, dec = gw_param['mass_1'], gw_param['mass_2'], gw_param['luminosity_distance'], gw_param['iota'], gw_param['psi'], gw_param['phase'], gw_param['geocent_time'], gw_param['ra'], gw_param['dec']
        zl                  = zl[pick_strongly_lensed]
        zs                  = zs[pick_strongly_lensed]
        n_images            = n_images[pick_strongly_lensed]
        x0_image_positions  = x0_image_positions[pick_strongly_lensed]
        x1_image_positions  = x1_image_positions[pick_strongly_lensed]
        kwargs_lenses       = kwargs_lenses[pick_strongly_lensed]
        mass_1              = mass_1[pick_strongly_lensed]
        mass_2              = mass_2[pick_strongly_lensed]
        luminosity_distance = luminosity_distance[pick_strongly_lensed]
        iota                = iota[pick_strongly_lensed]
        psi                 = psi[pick_strongly_lensed]
        phi                 = phi[pick_strongly_lensed]
        gamma1              = gamma1[pick_strongly_lensed]
        gamma2              = gamma2[pick_strongly_lensed]
        ra                  = ra[pick_strongly_lensed]
        dec                 = dec[pick_strongly_lensed]
        geocent_time        = geocent_time[pick_strongly_lensed]
        phase               = phase[pick_strongly_lensed]
        gw_param            = {'mass_1':mass_1, 'mass_2':mass_2, 'luminosity_distance':luminosity_distance, 'iota':iota, 'psi':psi, 'phase':phase, 'geocent_time':geocent_time, 'ra':ra, 'dec':dec, 'zs': zs}
        
        # Get the total number of lensed events
        number_of_lensed_events = len(pick_strongly_lensed)       # number of strong lensed events. Each event may contain multiple images
        
        #---------------------------------------------------#
        #                  Multiprocessing 2
        #---------------------------------------------------#
        # lensed_image_properties = (x0_image_position,x1_image_position,len(x0_image_position), kwargs_lens, zs)
        # note: zl and einstein_radius converted to 2d array and transpose before merging to lensed_image_properties
        iterations = np.arange(number_of_lensed_events)
        input_arguments = [(x0_image_positions[i],x1_image_positions[i],number_of_lensed_events,kwargs_lenses[i],zs[i], iterations[i], zl[i], einstein_radius[i]) for i in range(number_of_lensed_events)]  #np.concatenate((x0_image_positions,x1_image_positions,number_of_lensed_events,kwargs_lenses,zs, iterations, zl, einstein_radius))
        magnifications = np.ones((number_of_lensed_events, n_max_images))*np.nan
        time_delays   = np.ones((number_of_lensed_events, n_max_images))*np.nan
        determinants  = np.ones((number_of_lensed_events, n_max_images))*np.nan
        traces       = np.ones((number_of_lensed_events, n_max_images))*np.nan
        # FIXME: This won't work. The lensed images are not in the same order as the lensed_image_positions. They need to be ordered the same.
        # Get the image properties in parallel
        with Pool(processes=npool) as pool:
            # call the same function with different data in parallel
            for result in pool.map(self.get_image_properties, input_arguments):
                # report the value to show progress
                # return(magnifications,time_delays,determinant,trace,zs)
                # We assume that there are maximum of n_max_images images. If there are less, we fill the rest with nans. If this choice is confusing, know that this makes it a LOT easier to handle the data later on.
                n_images_iter = len(result[0])
                magnification = np.ones(n_max_images)*np.nan
                time_delay = np.ones(n_max_images)*np.nan
                determinant = np.ones(n_max_images)*np.nan
                trace = np.ones(n_max_images)*np.nan
                magnification[:n_images_iter] = result[0]
                time_delay[:n_images_iter]    = result[1]
                determinant[:n_images_iter]    = result[2]
                trace[:n_images_iter]          = result[3]
                # Add the magnifications, time delays, determinants, and traces to their respective arrays
                iter                           = result[5]
                magnifications[iter]           = magnification
                time_delays[iter]              = time_delay
                determinants[iter]             = determinant
                traces[iter]                   = trace
        #---------------------------------------------------#
        
        # image type classification
        image_type                  = np.zeros((number_of_lensed_events,n_max_images))
        image_type[traces < 0]       = 3
        image_type[traces > 0]       = 1
        image_type[determinants < 0] = 2
        
        # Get all of the signal to noise ratios
        snrs = self.get_snrs(gw_param, magnifications, time_delays)

        # probability of detection of each event. Each event can contain multiple images
        pdet = self.pdet_lensed(gw_param, magnifications, time_delays)
        
        # FIXME: Todo - add SNR
        return {'image_type':image_type, 'magnifications':magnifications, 'time_delays':time_delays,\
                'Pdet':pdet, 'theta_E':einstein_radius, 'zl': zl, 'zs':zs, 'snrs': snrs, 'x0_image_positions': x0_image_positions, \
                'x1_image_positions': x1_image_positions, 'kwargs_lenses': kwargs_lenses, 'n_images': n_images, \
                'determinants': determinants, 'traces': traces, 'gw_param': gw_param}
 
    ########################################################################
    #                                                                      #
    #                snr calculator for lensed events                    #
    #                                                                      #
    ########################################################################
    def get_snrs(self, gw_param, magnifications, time_delays):
        ''' 
        This function calculates the signal to noise ratio for each image in each event.
        -----------------
        Input parameters
        -----------------
        gw_param        : gravitational wave parameters dict (with 1d arrays(dtype=object)). 
                          { 'mass_1':mass_1, 'mass_2':mass_2, 'zs':zs, 'luminosity_distance':luminosity_distance, 'iota':theta_jn, 'psi':psi, 'phase':phase, \
                            'geocent_time':GPStimeValue, 'ra':ra, 'dec':dec  }
        magnifications  : magnification of all the lensed events. (1d array(dtype=object))
                          Each events can have multiple images, hence multiple magnification.
        time_delays      : time_delays of all the lensed events. (1d array)
        -----------------
        Return values
        -----------------
        snrs            : signal to noise ratio for each image in each event. (dictionary containing 'H1', 'L1', ..., and 'opt_snr_net', which is the network snr, for each image as an array with dimensions (number_of_lensed_events,n_max_images) )
        '''
        # Get the binary parameters
        number_of_lensed_events = len(magnifications)
        mass_1, mass_2, zs, luminosity_distance, iota, psi, phi, ra, dec, geocent_time, phase = gw_param['mass_1'], gw_param['mass_2'], gw_param['zs'], gw_param['luminosity_distance'], gw_param['iota'], gw_param['psi'], gw_param['phase'], gw_param['ra'], gw_param['dec'], gw_param['geocent_time'], gw_param['phase']
        # probability of detection (pdet) calculation for all the images. 
        # array is based on images not on lensed events. len(images)>len(events) 
        # Get the optimal signal-to-noise ratios for each image with n_max_images being the maximum
        detectors = self.quin_.list_of_detectors
        optimal_snrs = dict()
        optimal_snrs['opt_snr_net'] = np.ones((number_of_lensed_events, n_max_images))*np.nan
        for detector in detectors:
            optimal_snrs[detector] = np.ones((number_of_lensed_events, n_max_images))*np.nan
        for i in range(n_max_images):
            # Get the optimal signal to noise ratios for each image
            effective_luminosity_distance = luminosity_distance / np.sqrt(np.abs(magnifications[:,i]))
            effective_geocent_time = geocent_time + time_delays[:,i]
            # Each image has their own effective luminosity distance and effective geocent time
            optimal_snr = self.quin_.snr(mass_1, mass_2, effective_luminosity_distance, iota, psi, phase, effective_geocent_time, ra, dec) # Returns a dictionary
            optimal_snrs['opt_snr_net'][:,i] = optimal_snr['opt_snr_net']
            for detector in detectors:
                optimal_snrs[detector][:,i] = optimal_snr[detector]
        return optimal_snrs
    
    ########################################################################
    #                                                                      #
    #                pdet calculator for leansed events                    #
    #                                                                      #
    ########################################################################
    def pdet_lensed(self, gw_param, magnifications, time_delays):
        '''
        For calculating probability of detection of each image and also for each lensed events 
        -----------------
        Input parameters
        -----------------
        gw_param        : gravitational wave parameters dict (with 1d arrays(dtype=object)). 
                          { 'mass_1':mass_1, 'mass_2':mass_2, 'zs':zs, 'luminosity_distance':luminosity_distance, 'iota':theta_jn, 'psi':psi, 'phase':phase, \
                            'geocent_time':GPStimeValue, 'ra':ra, 'dec':dec  }
        magnifications  : magnification of all the lensed events. (1d array(dtype=object))
                          Each events can have multiple images, hence multiple magnification.
        time_delays      : time_delays of all the lensed events. (1d array)
        shape_          : number of images for each event. (1d array(dtype=int))
        -----------------
        Return values
        -----------------
        pdet_net        : probability of detection of any given image. Dimensions are: (number_of_lensed_events,n_max_images)

        '''
         # Get the binary parameters
        number_of_lensed_events = len(magnifications)
        mass_1, mass_2, zs, luminosity_distance, iota, psi, phi, ra, dec, geocent_time, phase = gw_param['mass_1'], gw_param['mass_2'], gw_param['zs'], gw_param['luminosity_distance'], gw_param['iota'], gw_param['psi'], gw_param['phase'], gw_param['ra'], gw_param['dec'], gw_param['geocent_time'], gw_param['phase']
        # probability of detection (pdet) calculation for all the images. 
        # array is based on images not on lensed events. len(images)>len(events) 
        # Get the optimal signal-to-noise ratios for each image with n_max_images being the maximum
        detectors = self.quin_.list_of_detectors
        pdet_net = np.ones((number_of_lensed_events, n_max_images))*np.nan
        for i in range(n_max_images):
            # Get the optimal signal to noise ratios for each image
            effective_luminosity_distance = luminosity_distance / np.sqrt(np.abs(magnifications[:,i]))
            effective_geocent_time = geocent_time + time_delays[:,i]
            # probability of detection (pdet) calculation for all the images. 
            # array is based on images not on lensed events. len(images)>len(events)
            param_ = {'mass_1':mass_1,'mass_2':mass_2,'luminosity_distance':effective_luminosity_distance,'iota':iota,'psi':psi,'phase':phase,'ra':ra,'dec':dec, 'geocent_time':effective_geocent_time}
            pdet_net[:,i] = self.quin_.pdet(param_)['pdet_net']
        return(pdet_net)
        
    ########################################################################
    #                                                                      #
    #         x-y position on the source plane (multiprocessing)           #
    #                                                                      #
    ########################################################################
    # FIXME: Change the lens parameters so that it's either a dictionary or separate variables (like an array of theta_E, e1, e2, etc). Dictionary might be the best thing to use here because we might want to change the function later on so that it uses a different lens model with a different set of parameters.
    def solve_lens_equation(self, lens_parameters):
        '''
        -----------------
        Input parameters
        -----------------
        lens_parameters       : np.array([x_source[0],x_source[1],e1,e2,lens_parameters['gamma'],gamma1,gamma2,lens_parameters['zl'],lens_parameters['zs']]).T
                       np.concatenate((items,lensModelList),axis=1)
                       # FIXME: It's not very clear here what this means
        -----------------
        Return values
        -----------------
        x0_image_position      : x position of image in the source plane
        x1_image_position      : y position of image in the source plane
        len(x0_image_position) : number of images for that lensed event 
        kwargs_lens  : various arguments associated with that lensed event 
        zs         : redshift of the source binary (use for sorting)
        '''
        zl = lens_parameters[7]
        zs = lens_parameters[8]
        iter = lens_parameters[9]
        lensModel = LensModel(lens_model_list = lens_parameters[10:].tolist(), # FIXME: Why is this 9? What does it mean? 
                          z_lens = zl,
                          z_source = zs )

        class_ = LensEquationSolver(lensModel)

        factor = 1.0
        kwargs_lens = [{'theta_E': 1.0, 'e1': lens_parameters[2], 'e2': lens_parameters[3], 'gamma': lens_parameters[4], \
                      'center_x': 0.0, 'center_y': 0.0}, {'gamma1': lens_parameters[5], 'gamma2': lens_parameters[6]}]
        x0_image_position,x1_image_position = class_.image_position_from_source( \
                                        sourcePos_x = lens_parameters[0]*factor, \
                                        sourcePos_y = lens_parameters[1]*factor, \
                                        kwargs_lens = kwargs_lens, solver='analytical' )
        nsamples = len(x0_image_position)
        return(x0_image_position, x1_image_position, nsamples, kwargs_lens, zs, iter)
    
    ########################################################################
    #                                                                      #
    #               lensed image parameters (multiprocessing)              #
    #                                                                      #
    ########################################################################
    def get_image_properties(self,parameters):
        '''
        -----------------
        Input parameters
        -----------------
        parameters          : np.concatenate((lensed_images, np.array([zl]).T, np.array([zs]).T, np.array([einstein_radius]).T), axis=1)
        -----------------
        Return values
        -----------------
        magnifications  : magnification of all the images
        time_delays      : time_delay of all the images
        determinant     : need for image type classification
        trace           : need for image type classification
        zs         : redshift of the source binary (use for sorting) 
        '''
        zl = parameters[6]
        iter = parameters[5]
        zs = parameters[4]
        einstein_radius = parameters[6]
        lensModel = LensModel(lens_model_list = ['EPL_NUMBA', 'SHEAR'], z_lens = zl, z_source = zs)

        x0_image_position, x1_image_position = parameters[0], parameters[1]
        kwargs_lens = parameters[3]
        imageLen = len(x0_image_position)
        theta_E = einstein_radius*np.ones(imageLen)

        # can have multiple magnification
        magnifications = lensModel.magnification(x0_image_position, x1_image_position, kwargs_lens)

        time_delays = lensModel.arrival_time(x0_image_position, x1_image_position, kwargs_lens)*theta_E**2

        # return: f_xx, f_xy, f_yx, f_yy components
        hessian = lensModel.hessian(x0_image_position, x1_image_position, kwargs_lens)
        determinant = np.array( (1 - hessian[0])*(1 - hessian[3]) - hessian[1]*hessian[2] )
        trace = np.array(2 - hessian[0] - hessian[3])
        
        return(magnifications,time_delays,determinant,trace,zs, iter)
    

'''
------------------------------------------------
  grand daughter class for rate comparision    
------------------------------------------------
'''
class rate_comparison(LensedCBCStatistics,UnlensedCBCStatistics):
    def __init__(self, npool=int(4), n_for_quintet=200, z_min=0.0001, z_max=10., m_min=4.59, m_max=86.22, event_type = "BBH", det_sensitivity='design',equal_mass=False):
        '''
        -----------------
        input parameters
        -----------------
        n_for_quintet   : number of bilby SNR use by quintet for interpolation
                          higher the number, more is the accuracy
                          no performance gain beyond n_for_quintet=1000
        z_min           : minimum value of redshift for BBH/BNS
                          example -> z_min=0.0001= 0.44312066Mpc =1.44*1e6ly
        z_max           : maximum value of redshift for BBH/BNS
                          example -> z_min=10.= 105999.14Mpc = 3.46*1e11ly
        m_min           : minimum mass of compact object for the chosen type of binaries (BBH or BNS)
        m_max           : maximum mass of compact object for the chosen type of binaries (BBH or BNS)
        event_type      : "BBH" (Binary black hole), "BNS" (Binary neutron star)
        det_sensitivity : sensitivity of the GW detectors
                          examples: "O1","O2","O3"
        '''
        super().__init__(npool=int(4),n_for_quintet=200, z_min=0.0001, z_max=10., m_min=4.59, m_max=86.22, event_type = "BBH", det_sensitivity='design',equal_mass=equal_mass)
        
        
    def ratio(self, size=1000):
        '''
        To find ratio of unlensed and lensed event rate
        -----------------
        Input parameters
        -----------------
        size          : size of the sample
        -----------------
        Return values
        -----------------
        rates         : dict containing lensed_rate_detectable, rate_detectable and ratio
        '''
        unlensed = np.array(self.rate_detectable(size=size))
        zs = self.gw_param['zs']
        lensed = np.array(self.lensed_rate_detectable(zs=zs))
        ratio_ = lensed/unlensed
        return({'lensed_rate_detectable':lensed,'rate_detectable':unlensed,'ratio':lensed/unlensed})
        

        
