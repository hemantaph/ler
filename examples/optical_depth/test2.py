# check the folder 'interpolator' exist
# if not, create the folder 'interpolator' and folder 'vel_disp' inside the folder 'interpolator'
# inside the folder 'interpolator/vel_disp' check if the vel_disp_init_dict.pickle exist
# if not, create the vel_disp_init_dict.pickle
# if yes, check inside the vel_disp_init_dict.pickle check if check the input dict is the same as one of the dict inside the pickle file
# If yes, get the dictionary number and check if the interpolator file exist.
# If yes, load the interpolator file
# If no, create the interpolator file and save it

import os
import pickle
import numpy as np
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
from ler.utils import create_inv_cdf
import test

class VelocityDispersionSampler():
    
    def __init__(self,
        nsamples_z=100,
        nsamples_sigma=50,
        param_dict_given={
        'z_min': 0.001,
        'z_max': 10,
        'cosmology': cosmo,
        }, 
        folder="./interpolator_pickle",
        ):
        self.param_dict_given = param_dict_given
        self.folder = folder
        self.nsamples_z = nsamples_z
        self.nsamples_sigma = nsamples_sigma

        path_inv_cdf, it_exist = self.interpolator_pickle_path(
            param_dict_given=self.param_dict_given,
            folder=self.folder,
            )

        if it_exist:
            print(
                f"Inveverse CDF of Velocity dispersion for inverse transform sampling will be loaded from {path_inv_cdf}"
            )
        else:
            print(
                f"Inveverse CDF of Velocity dispersion for inverse transform sampling will be generated at {path_inv_cdf}"
            )
            inv_cdf_list = self.init_inv_cdf_interpolation(nsamples_z, nsamples_sigma, path_inv_cdf)

        # load the interpolator
        with open(path_inv_cdf, "rb") as handle:
            self.inv_cdf_interpolator = pickle.load(handle)

    def sample_vel_disp(self, z, size):
        """
        Function to sample velocity dispersion from the interpolator
        """
        z_max = self.param_dict_given['z_max']
        z_min = self.param_dict_given['z_min']
        zlist = np.geomspace(z_min, z_max, self.nsamples_z)
        # find the index of z in zlist
        idx = np.searchsorted(zlist, z)
        # get the interpolator (inverse cdf) and sample
        u = np.random.uniform(0, 1, size=size)
        sample = self.inv_cdf_interpolator[idx](u)

        return sample


    def init_inv_cdf_interpolation(self, nsamples_z, nsamples_sigma, path_inv_cdf):
        """
        Function to create the interpolator pickle file path
        """
        # create data
        z = np.geomspace(0.001, 10, nsamples_z)
        z = np.ones((nsamples_z, nsamples_sigma))*z[:,None]
        sigma = np.linspace(0, 500, nsamples_sigma)

        inv_cdf_list = []
        for i in range(nsamples_z):
            phi = test.phi(sigma,z[i])
            # creating the interpolator
            inv_cdf_list.append(create_inv_cdf(x=sigma,y=phi)[2])

        # save the interpolator as pickle file
        with open(path_inv_cdf, "wb") as handle:
            pickle.dump(inv_cdf_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return inv_cdf_list
        
    def interpolator_pickle_path(self,
        param_dict_given={
        'z_min': 0.001,
        'z_max': 10,
        'cosmology': cosmo,
        },
        folder="./interpolator_pickle",
        ):
        """
        Function to create the interpolator pickle file path
        """
        
        # check the folder 'interpolator' exist
        if not os.path.exists(folder):
            os.makedirs(folder)
            os.makedirs(folder+'/vel_disp')
        else:
            if not os.path.exists(folder+'/vel_disp'):
                os.makedirs(folder+'/vel_disp')

        # check if param_dict_list.pickle exists
        path1 = folder+'/vel_disp' + '/vel_disp_init_dict.pickle'
        if not os.path.exists(path1):
            dict_list = []
            with open(path1, "wb") as handle:
                pickle.dump(dict_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # check if the input dict is the same as one of the dict inside the pickle file
        param_dict_stored = pickle.load(open(path1, "rb"))
        
        path2 = folder+'/vel_disp'
        len_ = len(param_dict_stored)
        if param_dict_given in param_dict_stored:
            idx = param_dict_stored.index(param_dict_given)
            # load the interpolator
            path_inv_cdf = path2 + "/inv_cdf_" + str(idx) + ".pickle"
            # there will be exception if the file is deleted by mistake
            if os.path.exists(path_inv_cdf):
                it_exist = True
            else:
                it_exist = False
        else:
            it_exist = False
            path_inv_cdf = path2 + "/inv_cdf_" + str(len_) + ".pickle"
            param_dict_stored.append(param_dict_given)
            with open(path1, "wb") as handle:
                pickle.dump(param_dict_stored, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return path_inv_cdf, it_exist

