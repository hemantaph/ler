import numpy as np
from scipy.stats import gengamma, rayleigh, norm
from scipy.interpolate import interp1d
from scipy.integrate import quad
from lenstronomy.Util.param_util import phi_q2_ellipticity
# for redshift to luminosity distance conversion
from astropy.cosmology import Planck18
import astropy.units as u
from astropy import constants as const
# the following .py file will be called if they are not given in the class initialization 
from source_population import CompactBinaryPopulation
from helperroutines import add_dictionaries_together, trim_dictionary
# for multiprocessing 
from multiprocessing import Pool
from tqdm import tqdm

# multiprocessing routines
import solve_lens_equation

min_images = 2
max_images = 5 

class LensGalaxyPopulation():
    def __init__(self, CompactBinaryPopulation_=False):
        '''
        class for lens galaxy population sampling
        Also includes functions to calculate lensed event rate
        Input parameters:
            CompactBinaryPopulation_ (class)    : already initialized CompactBinaryPopulation class (CompactBinaryPopulation for Source parameters sampling)
            z_min (float)                       : minimum redshift of the source population
            z_max (float)                       : maximum redshift of the source population
        Output parameters:
            None
        '''
        
        if CompactBinaryPopulation_==False:
            # initialization of clasess
            # CompactBinaryPopulation already inherits from Source_Galaxy_Population_Model class form source_population.py
            self.CompactBinaryPopulation = CompactBinaryPopulation(z_min=0.0001, z_max=10., m_min=4.59, m_max=86.22, event_type = "StellarBBH")
        else:
            # if the classes are already initialized, then just use them
            self.CompactBinaryPopulation = CompactBinaryPopulation_

        self.z_min = self.CompactBinaryPopulation.z_min
        self.z_max = self.CompactBinaryPopulation.z_max
        self.m_min = self.CompactBinaryPopulation.m_min
        self.m_max = self.CompactBinaryPopulation.m_max
        self.create_lookup_table(z_min=self.z_min, z_max=self.z_max)
        self.gw_param   = False 
        
        # To find the normalization constant of the pdf p(z)
        # this under the assumption that the event is strongly lensed
        # Define the merger-rate density function
        merger_rate_density_detector_frame = lambda z: self.CompactBinaryPopulation.merger_rate_density(z)/(1+z) 
        # Define the pdf p(z)
        pdf_unnormalized = lambda z: merger_rate_density_detector_frame(z) * self.differential_comoving_volume(z)*\
                                    self.strong_lensing_optical_depth(z)
        # Normalize the pdf
        # this normalization factor is common no matter what you choose for z_min and z_max
        self.normalization_pdf_z = quad(pdf_unnormalized, self.z_min, self.z_max)[0]
        return None

    def create_lookup_table(self, z_min, z_max):
        '''
        Functions to create lookup tables
        1. Redshift to co-moving distance
        2. Co-moving distance to redshift
        3. Redshift to luminosity distance
        4. Redshift to angular diameter distance
        5. Lens redshift sampler helper function
        6. Redshift to differential comoving volume
        Input parameters:
            z_min (float): minimum redshift of the source population
            z_max (float): maximum redshift of the source population
        Output parameters:
            None
        '''
        # initialing cosmological functions for fast calculation through interpolation
        z               = np.linspace(z_min,z_max,500) # red-shift
        Dc              = Planck18.comoving_distance(z).value # co-moving distance in Mpc
        self.z_to_Dc    = interp1d( z, Dc, kind = 'cubic')
        self.Dc_to_z    = interp1d( Dc, z, kind = 'cubic')
        self.z_to_luminosity_distance = self.CompactBinaryPopulation.z_to_luminosity_distance
        
        # for angular diameter distance
        quad_ = [] # refers to integration (with quad algorithm) from scipy 
        for ii in range(len(z)):
            quad_.append(quad(Planck18._inv_efunc_scalar, 0., z[ii], args=Planck18._inv_efunc_scalar_args)[0])
        quad_ = np.array(quad_)
        quad_int = interp1d(z, np.array(quad_), kind = 'cubic')

        H0d = Planck18._hubble_distance.value
        self.angular_diameter_distance_z1z2 = lambda zl0, zs0: H0d * (quad_int(zs0)-quad_int(zl0))/(zs0+1.)
        self.angular_diameter_distance = lambda zs0: H0d*quad_int(zs0)/(zs0+1.)

        # create a lookup table for the lens redshift draws
        r = np.linspace(0, 1, num = 100)
        u = 10*r**3 - 15*r**4 + 6*r**5 # See the integral of Eq. A7 of https://arxiv.org/pdf/1807.07062.pdf (cdf)
        self.lens_redshift_sampler_helper_function = interp1d(u, r, kind = 'cubic') # Computes r(u)
        
        # Create a lookup table for the differential comoving volume
        self.differential_comoving_volume = self.CompactBinaryPopulation.differential_comoving_volume
        return None

    def sample_lens_parameters(self, size=1000, lens_parameters_input={}):
        '''
        Function to sample galaxy lens parameters
        Input parameters:
            size : number of lens parameters to sample
            lens_parameters_input : dictionary of lens parameters to sample
        Output parameters:
            lens_parameters : dictionary of lens parameters and source parameters (lens conditions applied)
                            e.g. dictionary keys: 
                            lensing related=>['zl':redshift of lens, 'zs': redshift of source, 'sigma':velocity dispersion, 'q':axis ratios, 'e1':ellipticity, 'e2':ellipticity, 'gamma1':external-shear, 'gamma2':external-shear, 'Dl':angular diameter distance of lens, 'Ds':angular diameter distance of source, 'Dls':angular diameter distance between lens and source, 'theta_E': einstein radius in radian, 'gamma':spectral index of mass density distribution]
                            source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'iota': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination,]
                                    
        '''
        # Sample source redshifts from the source population
        # rejection sampled with optical depth
        source_params_strongly_lensed = self.sample_strongly_lensed_source_parameters(size=size)
        zs = source_params_strongly_lensed['zs']

        # Sample lens redshifts
        zl = self.sample_lens_redshifts(zs)
        # Sample velocity dispersions and axis ratios (note: these should be sampled together because the lensing probability depends on the combination of these two parameters)
        sigma, q = self.sample_velocity_dispersion_axis_ratio(zs)

        # Compute the Einstein radii
        theta_E = self.compute_einstein_radii(sigma, zl, zs)

        # Sample the axis ratio angle
        axis_ratio_angle_phi = self.sample_axis_ratio_angle_phi(size=size)

        # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
        e1, e2 = phi_q2_ellipticity(axis_ratio_angle_phi, q)

        # Sample shears
        gamma1, gamma2 = self.sample_galaxy_shear(size=size)

        # Sample the spectral index of the mass density distribution
        gamma = self.sample_gamma(size=size)

        # Compute the angular diameter distances
        Dl = self.angular_diameter_distance(zl) # for the lens
        Ds = self.angular_diameter_distance(zs) # for the source
        # Compute Dls also
        Dls = self.angular_diameter_distance_z1z2(zl, zs) # for the lens-source pair

        # Create a dictionary of the lens parameters
        lens_parameters = {'zl':zl, 'zs':zs, 'sigma':sigma, 'q':q, 'e1':e1, 'e2':e2, 'gamma1':gamma1, 'gamma2':gamma2, \
                           'Dl':Dl, 'Ds':Ds, 'Dls':Dls, 'theta_E':theta_E, 'gamma':gamma}
        
        # Add source params strongly lensed to the lens params
        lens_parameters.update(source_params_strongly_lensed)

        # Rejection sample based on the lensing probability, that is, rejection sample wrt theta_E
        mask = self.rejection_sample_lensing_probability(theta_E) # proportional to pi theta_E^2
        lens_parameters = {key:val[mask] for key, val in lens_parameters.items()}

        # Add the lensing parameter dictionaries together
        # Note: This is done after rejection sampling, so that we don't have to worry about the size of the dictionary
        lens_parameters = add_dictionaries_together(lens_parameters, lens_parameters_input)
        
        # Check if the lens are larger than requested size
        if len(lens_parameters['zl']) >= size:
            # Trim dicitionary to right size
            lens_parameters = trim_dictionary(lens_parameters, size)
            return lens_parameters
        else:
            # Run iteratively until we have the right number of lensing parmaeters
            return self.sample_lens_parameters(size=size, lens_parameters_input=lens_parameters)
        
    
    def sample_strongly_lensed_source_parameters(self, size=1000):
        '''
        Function to sample source redshifts and other parameters, conditioned on the source being strongly lensed
        PDF of source redshifts : R0(zs)/(1+zs) * dVc/dz * tau(zs)
        where R0(zs) is the source redshift distribution, dVc/dz is the comoving volume element, and tau(zs) is the strong lensing optical depth
        Sampling is done by rejection sampling
        Input parameters:
            size                        : number of source parameters to sample
        Output parameters:
            gw_param_strongly_lensed    : dictionary of source parameters, where zs is conditioned on the source being strongly lensed
                                        e.g. gw_param_strongly_lensed.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'zs']

        '''
        z_max = self.z_max
        z_min = self.z_min
        zs_strongly_lensed = []
        size_ = size
        while size_!=0:
            # get zs
            zs = self.CompactBinaryPopulation.sample_source_redshifts(size=size_, z_min=z_min, z_max=z_max)

            # put strong lensing condition with optical depth
            tau = self.strong_lensing_optical_depth(zs)
            tau_max = np.max(self.strong_lensing_optical_depth(z_max))
            r = np.random.uniform(0, tau_max, size=len(zs))
            pick_strongly_lensed = r<tau # pick strongly lensed sources
            # Add the strongly lensed source redshifts to the list
            zs_strongly_lensed += list(zs[pick_strongly_lensed]) # list concatenation
            
            size_ = abs(size-len(zs_strongly_lensed)) # loop until we have the right number of samples
                
        # gravitional waves source parameter sampling
        gw_param_strongly_lensed = \
        self.CompactBinaryPopulation.sample_gw_parameters(zs=np.array(zs_strongly_lensed), nsamples=size)
        
        return gw_param_strongly_lensed

    def sample_lens_redshifts(self, zs):
        '''
        Function to sample lens redshifts, conditioned on the lens being strongly lensed
        Input parameters:
            zs : source redshifts
        Output parameters:
            zl : lens redshifts
        '''
        # lens redshift distribution
        r = self.lens_redshift_sampler_helper_function(np.random.uniform(0, 1, size =len(zs) ))
        # comoing distance to the lens galaxy
        # on the condition that lens lie between the source and the observer
        lens_galaxy_Dc = self.z_to_Dc(zs)*r # corresponding element-wise multiplication between 2 arrays 
        # lens redshift
        zl = self.Dc_to_z(lens_galaxy_Dc) # 2D array
        return zl
    
    def sample_velocity_dispersion_axis_ratio(self, zs):
        '''
        Function to sample velocity dispersion and axis ratio of the lens galaxy
        Input parameters:
            zs : source redshifts
        Output parameters:
            sigma : velocity dispersion of the lens galaxy
            q : axis ratio of the lens galaxy
        '''
        size = len(zs)
        sigma = []
        q = []
        
        # Draw the velocity dispersions and axis ratios in 'chunks'
        size_ = size
        while size_!=0:
            # Draw the velocity dispersion
            a = gengamma.rvs(2.32/2.67, 2.67, size = size_)
            sigma_ = 161.*a
            
            # Draw the axis ratio see Appendix of https://arxiv.org/pdf/1807.07062.pdf
            s = 0.38 - 0.09177*a
            b = rayleigh.rvs(scale=s, size = size_)
            q_ = 1. - b
            
            # Weed out sigmas and axis ratios that have axis ratio below 0.2
            idx = q_ > 0.2
            sigma_ = sigma_[idx]
            q_ = q_[idx]
            
            # Append the velocity dispersions and axis ratios
            sigma += list(sigma_)
            q += list(q_)
            
            size_ = abs(size-len(sigma))
            
        # Transform to an array
        sigma = np.array(sigma)
        q = np.array(q)
        return sigma, q
    
    def compute_einstein_radii(self, sigma, zl, zs):
        '''
        Function to compute the Einstein radii of the lens galaxies
        Input parameters:
            sigma : velocity dispersion of the lens galaxy
            zl : lens redshifts
            zs : source redshifts
        Output parameters:
            theta_E : Einstein radii of the lens galaxies
        '''
        # Compute the angular diameter distances
        Ds = self.angular_diameter_distance(zs)
        Dls = self.angular_diameter_distance_z1z2(zl, zs)
        # Compute the Einstein radii
        theta_E = 4.*np.pi*(sigma/const.c.to('km/s').value)**2 * Dls/(Ds) # Note: km/s for sigma; Dls, Ds are in Mpc
        return theta_E
        
    def sample_axis_ratio_angle_phi(self, size=1000):
        '''
        Function to sample the axis rotation angle of the elliptical lens galaxy
        Input parameters:
            size : number of samples
        Output parameters:
            q : axis rotation angle
        '''
        # Draw the angles
        phi = np.random.uniform(0, 2*np.pi, size=size)
        return phi
    
    def sample_galaxy_shear(self, size):
        '''
        Function to sample the lens galaxy shear
        Input parameters:
            size : number of samples
        Output parameters:
            gamma_1 : shear component in the x-direction
            gamma_2 : shear component in the y-direction
        '''
        # Draw an external shear
        gamma_1 = norm.rvs(size=size, scale = 0.05)
        gamma_2 = norm.rvs(size=size, scale = 0.05)
        return gamma_1, gamma_2
    
    def sample_gamma(self, size=1000):
        '''
        Function to sample the lens galaxy spectral index of the density profile
        Input parameters:
            size : number of samples
        Output parameters:
            gamma : spectral index of the density profile
        '''
        self.gamma_mean = 2
        self.gamma_std = 0.2
        return np.random.normal(loc = self.gamma_mean, scale = self.gamma_std, size = size) 
    
    def rejection_sample_lensing_probability(self, theta_E):
        '''
        Function to conduct rejection sampling wrt einstein radius
        Input parameters:
            theta_E : einstein radius
        Output parameters:
            idx : boolean array of size len(theta_E) indicating whether the sample is accepted or not
        '''
        size = len(theta_E)
        theta_E_max = np.max(theta_E)
        u = np.random.uniform(0, theta_E_max**2, size=size)
        idx = u < theta_E**2
        return idx


    def strong_lensing_optical_depth(self, zs):
        '''
        Function to compute the strong lensing optical depth
        Input parameters:
            zs : source redshifts
        Output parameters:
            tau : strong lensing optical depth
        '''
        # z to luminosity_distance (luminosity_distance) conversion
        Dc = self.z_to_Dc(zs)*1e-3  # 1e-3 converts Mpc to Gpc
        return( (Dc/62.2)**3 )

    def get_caustics(self, theta_E, gamma, gamma1, gamma2, e1, e2):
        '''
        function to get the caustics of the lens
        Input parameters:
            theta_E : Einstein radius
            gamma : spectral index of the density profile
            gamma1 : shear component in the x-direction
            gamma2 : shear component in the y-direction
            e1 : ellipticity component in the x-direction
            e2 : ellipticity component in the y-direction
        Output parameters:
            caustic_double : double caustic
            caustic_diamond : diamond caustic
        '''
        # set up the lens model
        kwargs_lens = [{'theta_E': theta_E, 'e1': e1, 'e2': e2, 'gamma': gamma, 'center_x': 0.0, 'center_y': 0.0}, 
                       {'gamma1': gamma1, 'gamma2': gamma2, 'ra_0': 0, 'dec_0':0}]
        # Get the lensing diamond and double caustics
        caustic_double_points = caustics_epl_shear(kwargs_lens, return_which='double', maginf=-100)
        caustic_diamond_points = caustics_epl_shear(kwargs_lens, return_which='caustic', maginf=-100)
        caustic_double = Polygon(caustic_double_points.T)
        caustic_diamond = Polygon(caustic_diamond_points.T)
        return caustic_diamond, caustic_double
    
    def get_image_properties(self, lens_parameters, lensModelList=['EPL_NUMBA', 'SHEAR'], npool=4):
        '''
        Function to get the image properties
        Input parameters:
            lens_parameters : dictionary of lens parameters
                            e.g. lens_parameters.keys() = ['zs', 'zl', 'gamma1', 'gamma2', 'e1', 'e2', 'gamma', 'theta_E']
            lensModelList   : list of lens models
                            e.g. lensModelList = ['EPL_NUMBA', 'SHEAR']
            npool           : number of processes to use
                            to check the number of cores available, use os.cpu_count()
        Output parameters:
            lens_parameters : dictionary of lens parameters and image properties
                            e.g. lens_parameters.keys() = ['zs', 'zl', 'gamma1', 'gamma2', 'e1', 'e2', 'gamma', 'theta_E', 'x_image', 'y_image', 'theta_E_image', 'gamma_image', 'e1_image', 'e2_image', 'phi_image', 'caustic_diamond', 'caustic_double']
        '''
        n_max_images     = max_images
        zs               = lens_parameters['zs']
        size             = len(zs)
        zl               = lens_parameters['zl']
        # external shear params to the 'PEMD' galaxy lens
        gamma1, gamma2   = lens_parameters['gamma1'], lens_parameters['gamma2']
        # ellipticity of the galaxy lens
        e1, e2           = lens_parameters['e1'], lens_parameters['e2']
        gamma            = lens_parameters['gamma']
        einstein_radius  = lens_parameters['theta_E']
        # Create the lens model list (note: can be a different lens model for different samples)
        lensModelList    = np.array(lensModelList)*np.ones((size,len(lensModelList)),dtype=object) 
        
        # get image properties (with Multiprocessing)
        iterations = np.arange(size)
        input_arguments = np.array([e1,e2,gamma,gamma1,gamma2,zl,zs,einstein_radius,iterations], dtype=object).T
        input_arguments = np.concatenate((input_arguments,lensModelList),axis=1)
        # Initialize the image positions and lens argument list here.
        x0_image_positions = np.ones((size, n_max_images))*np.nan
        x1_image_positions = np.ones((size, n_max_images))*np.nan
        magnifications = np.ones((size, n_max_images))*np.nan
        time_delays   = np.ones((size, n_max_images))*np.nan
        determinants  = np.ones((size, n_max_images))*np.nan
        traces       = np.ones((size, n_max_images))*np.nan
        n_images = np.ones(size, dtype=int)
        x_source, y_source = np.ones(size)*np.nan, np.ones(size)*np.nan
        eta = np.ones(size)*np.nan
        phi = np.ones(size)*np.nan
        weights = np.ones(size)*np.nan
        
        # Solve the lens equation
        print('solving lens equations, started...')
        with Pool(processes=npool) as pool:
            # call the same function with different data in parallel
            # imap->retain order in the list, while map->doesn't
            for result in tqdm(pool.imap(solve_lens_equation.solve_lens_equation,input_arguments), total=len(input_arguments), ncols= 100, disable=False): 
                #print(result)
                '''
                for i in tqdm(range(size)):
                    result = self.solve_lens_equation(input_arguments[i])
                '''
                x_source_i, y_source_i, eta_i, phi_i, x0_image_position_i, x1_image_position_i, magnifications_i,time_delays_i, \
                n_image_i, determinant_i, trace_i, iter_i, weights_i = result

                n_image_i = min(n_image_i, n_max_images)
                n_images[iter_i] = n_image_i
                x0_image_position = np.ones(n_max_images)*np.nan
                x1_image_position = np.ones(n_max_images)*np.nan
                x0_image_position[:n_image_i] = x0_image_position_i[:n_image_i]
                x1_image_position[:n_image_i] = x1_image_position_i[:n_image_i]
                x0_image_positions[iter_i] = x0_image_position # shape = (size, n_max_images)
                x1_image_positions[iter_i] = x1_image_position # shape = (size, n_max_images)
                magnification = np.ones(n_max_images)*np.nan
                time_delay = np.ones(n_max_images)*np.nan
                determinant = np.ones(n_max_images)*np.nan
                trace = np.ones(n_max_images)*np.nan
                magnification[:n_image_i] = magnifications_i[:n_image_i]
                time_delay[:n_image_i]    = time_delays_i[:n_image_i]
                determinant[:n_image_i]    = determinant_i[:n_image_i]
                trace[:n_image_i]          = trace_i[:n_image_i]
                # Add the magnifications, time delays, determinants, and traces to their respective arrays
                magnifications[iter_i]           = magnification
                time_delays[iter_i]              = time_delay
                determinants[iter_i]             = determinant
                traces[iter_i]                   = trace
                x_source[iter_i]                 = x_source_i
                y_source[iter_i]                 = y_source_i
                eta[iter_i]                     = eta_i
                phi[iter_i]                     = phi_i
                weights[iter_i]                  = weights_i
        print('solving lens equations, ended...')   
        
        # time-delays: convert to positive values
        # time-delays will be relative to the first arrived signal of an lensed event
        time_delays = time_delays - np.array([np.sort(time_delays,axis=1)[:,0]]).T
        
        # select only strongly lensed events are selected
        assert np.all(n_images>=2), "There are events with no images!"
        
        # image type classification (morse phase)
        number_of_lensed_events = size
        image_type                   = np.zeros((number_of_lensed_events,n_max_images))
        image_type[traces < 0]       = 3
        image_type[traces > 0]       = 1
        image_type[determinants < 0] = 2
        
        # Return a dictionary with all of the lens information but also the BBH parameters from gw_param
        image_parameters = {'n_images':n_images, 'x0_image_positions':x0_image_positions, 'x1_image_positions':x1_image_positions, 'magnifications':magnifications, 'time_delays':time_delays, 'traces':traces, 'determinants':determinants, 'image_type':image_type}
        lens_parameters.update(image_parameters)
        
        return lens_parameters

    def get_lensed_snrs(self, snr_calculator, lensed_param):
        ''' 
        Function to calculate the signal to noise ratio for each image in each event.
        Input parameters:
            snr_calculator : snr_calculator class
                            this is an already initialized class that contains a 
                            function (snr_calculator.snr) that actually calculates snr with the given gw_params.
                            snr_calculator.list_of_detectors is a list of detectors that are used in the calculation.
                            e.g. ['H1', 'L1', 'V1']
            lensed_param   : dictionary containing the both already lensed source paramters and image parameters. 
                            e.g. lensed_param.keys() = ['mass_1', 'mass_2', 'zs', 'luminosity_distance', 'iota', 'psi', 'phi', 'ra', 'dec', 'geocent_time', 'phase', 'magnifications', 'time_delays']
        Output parameters:
            snrs            : signal to noise ratio for each image in each event.
                            (dictionary containing 'H1', 'L1', ..., and 'opt_snr_net', which is the network snr, 
                            for each image as an array with dimensions (number_of_lensed_events,n_max_images) )
                            e.g. snrs.keys() = ['H1', 'L1', 'V1', 'opt_snr_net']
        Example:
            snrs = get_lensed_snrs(, lensed_param)
                            
        '''
        # needed to calculate effective luminosity distance and effective time delay
        magnifications = lensed_param['magnifications']
        time_delays    = lensed_param['time_delays']

        n_max_images   = max_images
        # Get the binary parameters
        number_of_lensed_events = len(magnifications)
        mass_1, mass_2, luminosity_distance, iota, psi, ra, dec, geocent_time, phase = \
        lensed_param['mass_1'], lensed_param['mass_2'], lensed_param['luminosity_distance'], lensed_param['iota'], \
        lensed_param['psi'], lensed_param['ra'], lensed_param['dec'], lensed_param['geocent_time'], lensed_param['phase']
        
        # setting up snr dictionary
        detectors = snr_calculator.list_of_detectors
        optimal_snrs = dict()
        optimal_snrs['opt_snr_net'] = np.ones((number_of_lensed_events, n_max_images))*np.nan
        for detector in detectors:
            optimal_snrs[detector] = np.ones((number_of_lensed_events, n_max_images))*np.nan
        
        # LALSimulation cannot handle NaN
        if snr_calculator.snr_type == 'inner_product':
            print('There will be {} progress bar iteration'.format(n_max_images))

        for i in range(n_max_images):
            # Get the optimal signal to noise ratios for each image
            buffer = magnifications[:,i]
            idx = ~np.isnan(buffer) # index of not-nan   
            effective_luminosity_distance = luminosity_distance[idx] / np.sqrt(np.abs(buffer[idx]))
            effective_geocent_time = geocent_time[idx] + time_delays[idx,i]
            # if GPS time is negative, shift it 
            # by a year until it is positive
            effective_geocent_time[effective_geocent_time < 0] += 31556952 # number of seconds in a year
                    
            # Each image has their own effective luminosity distance and effective geocent time
            if len(effective_luminosity_distance)!=0:
                # Returns a dictionary
                optimal_snr = snr_calculator.snr(mass_1[idx], mass_2[idx], effective_luminosity_distance, iota[idx], \
                                            psi[idx], phase[idx], effective_geocent_time, ra[idx], dec[idx], jsonFile=False) 
                
                optimal_snrs['opt_snr_net'][idx,i] = optimal_snr['opt_snr_net']
                for detector in detectors:
                    optimal_snrs[detector][idx,i] = optimal_snr[detector]
                
        self.optimal_snrs = optimal_snrs
        return optimal_snrs
