
from .utils import interpolator_json_path, cubic_spline_interpolator, pdf_cubic_spline_interpolator, cubic_spline_interpolator2d_array, inverse_transform_sampler, pdf_cubic_spline_interpolator2d_array, save_json, load_json, inverse_transform_sampler2d
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.interpolate import  CubicSpline
from scipy.stats import gaussian_kde
from numba import njit


class FunctionConditioning():

    def __init__(self,
        function=None,  # can also be an array of function values 
        x_array=None,
        conditioned_y_array=None,  # if this is not none, 2D interpolation will be used
        y_array=None,
        non_zero_function=False,
        gaussian_kde=False,
        gaussian_kde_kwargs={},
        identifier_dict={},
        directory='./interpolator_json',
        sub_directory='default',
        name='default',
        create_new=False,
        create_function=False,
        create_function_inverse=False,
        create_pdf=False,
        create_rvs=False,
        multiprocessing_function=False,
        callback=None,
    ):
        create = self.create_decision_function(create_function, create_function_inverse, create_pdf, create_rvs)
        self.info = identifier_dict.copy()
        self.callback = callback
        
        if create:
            # create_interpolator input list
            input_list = [function, x_array, conditioned_y_array, create_function_inverse, create_pdf, create_rvs, multiprocessing_function]
            input_list_kde = [x_array, y_array, gaussian_kde_kwargs]

            # check first whether the directory, subdirectory and pickle exist
            path_inv_cdf, it_exist = interpolator_json_path(
                identifier_dict=identifier_dict,
                directory=directory,
                sub_directory=sub_directory,
                interpolator_name=name,
            )

            # if the interpolator exists, load it
            if create_new:
                it_exist = False

            if it_exist:
                print(f"{name} interpolator will be loaded from {path_inv_cdf}")
                # load the interpolator, where the interpolator is a dictionary
                # interpolator = {
                #     'x_array': x_array,
                #     'z_array': z_array,
                #     'conditioned_y_array': conditioned_y_array,
                #     'function_spline': function_spline,
                #     'function_inverse_spline': function_inverse_spline,
                #     'pdf_norm_const': pdf_norm_const,
                #     'cdf_values': cdf_values,
                # }
                interpolator = load_json(path_inv_cdf)
            else:
                print(f"{name} interpolator will be generated at {path_inv_cdf}")
                if not gaussian_kde:
                    interpolator = self.create_interpolator(*input_list)
                    # save the interpolator
                    save_json(path_inv_cdf, interpolator)
                else:
                    # gaussian KDE
                    interpolator = self.create_gaussian_kde(*input_list_kde)

                    #------------
                    # check scipy gaussian kde can be saved as json
                    #------------
                    # save_json(path_inv_cdf, interpolator)

            if not gaussian_kde:
                # Convert loaded JSON lists to NumPy arrays for Numba compatibility
                x_array = np.array(interpolator['x_array'])
                z_array = np.array(interpolator['z_array'])
                conditioned_y_array = np.array(interpolator['conditioned_y_array']) if interpolator['conditioned_y_array'] is not None else None
                y_array = None
                function_spline = np.array(interpolator['function_spline'])
                function_inverse_spline = np.array(interpolator['function_inverse_spline']) if interpolator['function_inverse_spline'] is not None else None
                pdf_norm_const = np.array(interpolator['pdf_norm_const']) if interpolator['pdf_norm_const'] is not None else None
                cdf_values = np.array(interpolator['cdf_values']) if interpolator['cdf_values'] is not None else None

                if conditioned_y_array is None:
                    # function is 1D
                    # njit(lambda x: cubic_spline_interpolator(x, function_spline, x_array))
                    function_any = lambda x: cubic_spline_interpolator(x, function_spline, x_array)

                    def function_non_zero(x):
                        result = cubic_spline_interpolator(x, function_spline, x_array)
                        idx = result<0.0
                        result[idx] = 0.0
                        return result
                    
                    function_final = function_any if non_zero_function else function_non_zero

                    self.function = njit(function_final) if create_function else None

                    # inverse function is 1D
                    function_any = lambda x: cubic_spline_interpolator(x, function_inverse_spline, z_array)

                    def function_non_zero(x):
                        result = cubic_spline_interpolator(x, function_inverse_spline, z_array)
                        idx = result<0.0
                        result[idx] = 0.0
                        return result

                    function_final = function_any if non_zero_function else function_non_zero

                    self.function_inverse = njit(function_final) if create_function_inverse else None

                    # pdf is 1D
                    self.pdf = njit(lambda x: pdf_cubic_spline_interpolator(x, pdf_norm_const, function_spline, x_array)) if create_pdf else None
                    # sampler is 1D
                    self.rvs = njit(lambda size: inverse_transform_sampler(size, cdf_values, x_array)) if create_rvs else None
                    
                else:
                    # function is 2D
                    function_any = lambda x, y: cubic_spline_interpolator2d_array(x, y, function_spline, x_array, conditioned_y_array)

                    def function_non_zero(x, y):
                        result = cubic_spline_interpolator2d_array(x, y, function_spline, x_array, conditioned_y_array)
                        idx = result<0.0
                        result[idx] = 0.0
                        return result
                    
                    function_final = function_any if non_zero_function else function_non_zero
                    self.function = njit(function_final) if create_function else None

                    # inverse function is 2D
                    function_any = lambda x, y: cubic_spline_interpolator2d_array(x, y, function_inverse_spline, z_array, conditioned_y_array)

                    def function_non_zero(x, y):
                        result = cubic_spline_interpolator2d_array(x, y, function_inverse_spline, z_array, conditioned_y_array)
                        idx = result<0.0
                        result[idx] = 0.0
                        return result
                    
                    function_final = function_any if non_zero_function else function_non_zero
                    self.function_inverse = njit(function_final) if create_function_inverse else None

                    self.pdf = njit(lambda x, y: pdf_cubic_spline_interpolator2d_array(x, y, pdf_norm_const, function_spline, x_array, conditioned_y_array)) if create_pdf else None

                    self.rvs = njit(lambda size, y: inverse_transform_sampler2d(size, y, cdf_values, x_array, conditioned_y_array)) if create_rvs else None

                self.x_array = x_array
                self.z_array = z_array
                self.conditioned_y_array = conditioned_y_array
                self.function_spline = function_spline
                self.function_inverse_spline = function_inverse_spline
                self.pdf_norm_const = pdf_norm_const
                self.cdf_values = cdf_values
            else:
                x_array = interpolator['x_array']
                y_array = interpolator['y_array']
                kde_object = interpolator['kde_object']

                
                self.pdf = lambda x: kde_object.pdf(x) if create_pdf else None
                if y_array is None:
                    self.rvs = lambda size: kde_object.resample(size)[0] if create_rvs else None
                else:
                    self.rvs = lambda size: kde_object.resample(size) if create_rvs else None

                self.x_array = x_array
                self.y_array = y_array
                self.kde_object = kde_object
        else:
            self.conditioned_y_array = None
            self.x_array = None
            self.y_array = None
            self.kde_object = None

    def __call__(self, *args):
        args = [np.array(arg).reshape(-1) if isinstance(arg, float) else arg for arg in args]
        return getattr(self, self.callback)(*args)
            

    def create_decision_function(self, create_function, create_function_inverse, create_pdf, create_rvs):

        decision_bool = True
        if not isinstance(create_function, bool) and callable(create_function):
            self.function = create_function
            decision_bool = False
        if not isinstance(create_function_inverse, bool) and callable(create_function_inverse):
            self.function_inverse = create_function_inverse
            decision_bool = False
        if not isinstance(create_pdf, bool) and callable(create_pdf):
            self.pdf = create_pdf
            decision_bool = False
        if not isinstance(create_rvs, bool) and callable(create_rvs):
            self.rvs = create_rvs
            decision_bool = False
        
        return decision_bool


    def create_gaussian_kde(self, x_array, y_array, gaussian_kde_kwargs):

        # 1d KDE
        if y_array is None:
            kde = gaussian_kde(x_array, **gaussian_kde_kwargs)
        else:
            data = np.vstack([x_array, y_array])
            kde = gaussian_kde(data, **gaussian_kde_kwargs)

        return {
            'x_array': x_array,
            'y_array': y_array,
            'kde_object': kde,
        }

    def create_interpolator(self, function, x_array, conditioned_y_array, create_function_inverse, create_pdf, create_rvs, multiprocessing_function):

        # function can be numpy array or callable
        # x_array, z_array are 2D arrays if conditioned_y_array is not None
        x_array, z_array, conditioned_y_array = self.create_z_array(x_array, function, conditioned_y_array, create_pdf, create_rvs, multiprocessing_function)
        del function

        function_spline = self.function_spline_generator(x_array, z_array, conditioned_y_array)

        if create_function_inverse:
            if conditioned_y_array is None:
                idx_sort = np.argsort(z_array)
            else:
                idx_sort = np.argsort(z_array, axis=1)
            x_array_ = x_array[idx_sort]
            z_array_ = z_array[idx_sort]

            # check z_array is strictly increasing
            # if (not np.all(np.diff(z_array) > 0)) or (not np.all(np.diff(z_array) < 0)):
            #     raise ValueError("z_array must be strictly increasing")
            
            function_inverse_spline = self.function_spline_generator(z_array_, x_array_, conditioned_y_array)
        else:
            function_inverse_spline = None

        if create_pdf or create_rvs:
            # cannot have -ve pdf
            pdf_norm_const = self.pdf_norm_const_generator(x_array, function_spline, conditioned_y_array)

            if create_rvs:
                cdf_values = self.cdf_values_generator(x_array, z_array, conditioned_y_array)
            else:
                cdf_values = None
        else:
            pdf_norm_const = None
            cdf_values = None
            
        return {
            'x_array': x_array,
            'z_array': z_array,
            'conditioned_y_array': conditioned_y_array,
            'function_spline': function_spline,
            'function_inverse_spline': function_inverse_spline,
            'pdf_norm_const': pdf_norm_const,
            'cdf_values': cdf_values,
        }

    def create_z_array(self, x_array, function, conditioned_y_array, create_pdf, create_rvs, multiprocessing_function):

        if callable(function):
            # 1D
            if conditioned_y_array is None:
                z_array = function(x_array)
                # remove nan values
                idx = np.argwhere(np.isnan(z_array))
                x_array = np.delete(x_array, idx)
                z_array = np.delete(z_array, idx)
            # 2D
            else:
                # check if x_array is 2D, if not, make it 2D of shape (len(conditioned_y_array), len(x_array))
                if x_array.ndim == 1:
                    x_array = np.array([x_array]*len(conditioned_y_array))

                idx = np.argsort(conditioned_y_array)
                conditioned_y_array = conditioned_y_array[idx]
                # x_array is 2D here, each row corresponds to a different y value
                x_array = x_array[idx]
                # sort each row of x_array
                x_array = np.sort(x_array, axis=1)

                if multiprocessing_function:
                    z_array = function(x_array, conditioned_y_array)
                else:
                    z_list = []
                    for i, y in enumerate(conditioned_y_array):
                        try:
                            z_list.append(function(x_array[i], y*np.ones_like(x_array[i])))
                        except:
                            # print(x_array[i], y)
                            z_list.append(function(x_array[i], y))
                    z_array = np.array(z_list)
        else:
            if conditioned_y_array is None:
                z_array = function
                # remove nan values
                idx = np.argwhere(np.isnan(z_array))
                x_array = np.delete(x_array, idx)
                z_array = np.delete(z_array, idx)
            else:
                if x_array.ndim == 1:
                    x_array = np.array([x_array]*len(conditioned_y_array))
                if function.ndim == 1:
                    raise ValueError('function must be 2D array if conditioned_y_array is not None')
                # row sort
                idx = np.argsort(conditioned_y_array)
                conditioned_y_array = conditioned_y_array[idx]
                # x_array is 2D here, each row corresponds to a different y value
                x_array = x_array[idx]
                z_array = function[idx]

                z_list = []
                x_list = []
                for i in range(len(conditioned_y_array)):
                    # column sort
                    idx = np.argsort(x_array[i])
                    x_list.append(x_array[i][idx])
                    z_list.append(z_array[i][idx])
                x_array = np.array(x_list)
                z_array = np.array(z_list)

        # cannot have -ve pdf
        if create_pdf or create_rvs:
            z_array[z_array < 0.0] = 0.0

        return x_array, z_array, conditioned_y_array

    def cdf_values_generator(self, x_array, z_array, conditioned_y_array):
        # 1D case
        if conditioned_y_array is None:
            # z_array[z_array<0.] = 0. # already done
            cdf_values = cumulative_trapezoid(z_array, x_array, initial=0)
            cdf_values = cdf_values/cdf_values[-1]
        # 2D case
        else:
            cdf_values = []
            for i, y in enumerate(conditioned_y_array):
                z_array_ = z_array[i]
                z_array_[z_array_<0.] = 0.
                cdfs_ = cumulative_trapezoid(z_array_, x_array[i], initial=0)

                cdf_values.append(cdfs_/cdfs_[-1])
                # cdf_values.append(cdfs_)
        
        return np.array(cdf_values)
    
    def pdf_norm_const_generator(self, x_array, function_spline, conditioned_y_array):

        # 1D case
        if conditioned_y_array is None:
            # pdf_unorm = lambda x: cubic_spline_interpolator(np.array([x]), function_spline, x_array)
            def pdf_unorm(x):
                result = cubic_spline_interpolator(np.array([x]), function_spline, x_array)
                idx = result<0.
                result[idx] = 0.0

                return result

            norm = quad(pdf_unorm, min(x_array), max(x_array))[0]
            return norm
        # 2D case
        else:
            norm = []
            for i, y in enumerate(conditioned_y_array):
                # pdf_unorm = lambda x: cubic_spline_interpolator(np.array([x]), function_spline[i], x_array[i])
                def pdf_unorm(x):
                    result = cubic_spline_interpolator(np.array([x]), function_spline[i], x_array[i])
                    idx = result<0.
                    result[idx] = 0.0

                    return result

                norm.append(quad(pdf_unorm, min(x_array[i]), max(x_array[i]))[0])

            return np.array(norm)
    
    def function_spline_generator(self, x_array, z_array, conditioned_y_array):
        # 1D case
        if conditioned_y_array is None:
            function_spline = CubicSpline(x_array, z_array).c
        # 2D case
        else:
            function_spline = []
            for i, y in enumerate(conditioned_y_array):
                
                function_spline.append(CubicSpline(x_array[i], z_array[i]).c)
            
        return np.array(function_spline)