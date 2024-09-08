# -*- coding: utf-8 -*-
""" 
This module contains helper routines for other modules in the ler package.
"""

import os
import pickle
import numpy as np
import json
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy.integrate import quad, cumtrapz
from numba import njit
# import datetime


class NumpyEncoder(json.JSONEncoder):
    """
    Class for storing a numpy.ndarray or any nested-list composition as JSON file. This is required for dealing np.nan and np.inf.

    Parameters
    ----------
    json.JSONEncoder : `class`
        class for encoding JSON file

    Returns
    ----------
    json.JSONEncoder.default : `function`
        function for encoding JSON file

    Example
    ----------
    >>> import numpy as np
    >>> import json
    >>> from ler import helperroutines as hr
    >>> # create a dictionary
    >>> param = {'a': np.array([1,2,3]), 'b': np.array([4,5,6])}
    >>> # save the dictionary as json file
    >>> with open('param.json', 'w') as f:
    >>>     json.dump(param, f, cls=hr.NumpyEncoder)
    >>> # load the dictionary from json file
    >>> with open('param.json', 'r') as f:
    >>>     param = json.load(f)
    >>> # print the dictionary
    >>> print(param)
    {'a': [1, 2, 3], 'b': [4, 5, 6]}
    """

    def default(self, obj):
        """function for encoding JSON file"""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def load_json(file_name):
    """Load a json file.

    Parameters
    ----------
    file_name : `str`
        json file name for storing the parameters.

    Returns
    ----------
    param : `dict`
    """
    with open(file_name, "r", encoding="utf-8") as f:
        param = json.load(f)

    return param

def save_json(file_name, param):
    """Save a dictionary as a json file.

    Parameters
    ----------
    file_name : `str`
        json file name for storing the parameters.
    param : `dict`
        dictionary to be saved as a json file.
    """
    with open(file_name, "w", encoding="utf-8") as write_file:
        json.dump(param, write_file)

def append_json(file_name, new_dictionary, old_dictionary=None, replace=False):
    """
    Append (values with corresponding keys) and update a json file with a dictionary. There are four options:

    1. If old_dictionary is provided, the values of the new dictionary will be appended to the old dictionary and save in the 'file_name' json file.
    2. If replace is True, replace the json file (with the 'file_name') content with the new_dictionary.
    3. If the file does not exist, create a new one with the new_dictionary.
    4. If none of the above, append the new dictionary to the content of the json file.

    Parameters
    ----------
    file_name : `str`
        json file name for storing the parameters. 
    new_dictionary : `dict`
        dictionary to be appended to the json file.
    old_dictionary : `dict`, optional
        If provided the values of the new dictionary will be appended to the old dictionary and save in the 'file_name' json file. 
        Default is None.
    replace : `bool`, optional
        If True, replace the json file with the dictionary. Default is False.

    """

    # check if the file exists
    # time
    # start = datetime.datetime.now()
    if old_dictionary:
        data = old_dictionary
    elif replace:
        data = new_dictionary
    elif not os.path.exists(file_name):
        #print(f" {file_name} file does not exist. Creating a new one...")
        replace = True
        data = new_dictionary
    else:
        #print("getting data from file")
        with open(file_name, "r", encoding="utf-8") as f:
            data = json.load(f)
    # end = datetime.datetime.now()
    # print(f"Time taken to load the json file: {end-start}")

    # start = datetime.datetime.now()
    if not replace:
        data = add_dictionaries_together(data, new_dictionary)
        # data_key = data.keys()
        # for key, value in new_dictionary.items():
        #     if key in data_key:
        #         data[key] = np.concatenate((data[key], value)).tolist()
    # end = datetime.datetime.now()
    # print(f"Time taken to append the dictionary: {end-start}")

    # save the dictionary
    # start = datetime.datetime.now()
    #print(data)
    with open(file_name, "w", encoding="utf-8") as write_file:
        json.dump(data, write_file, indent=4, cls=NumpyEncoder)
    # end = datetime.datetime.now()
    # print(f"Time taken to save the json file: {end-start}")

    return data

# def add_dict_values(dict1, dict2):
#     """Adds the values of two dictionaries together.
    
#     Parameters
#     ----------
#     dict1 : `dict`
#         dictionary to be added.
#     dict2 : `dict`
#         dictionary to be added.

#     Returns
#     ----------
#     dict1 : `dict`
#         dictionary with added values.
#     """
#     data_key = dict1.keys()
#     for key, value in dict2.items():
#         if key in data_key:
#             dict1[key] = np.concatenate((dict1[key], value))

    return dict1

def get_param_from_json(json_file):
    """
    Function to get the parameters from json file.

    Parameters
    ----------
    json_file : `str`
        json file name for storing the parameters.

    Returns
    ----------
    param : `dict`
    """
    with open(json_file, "r", encoding="utf-8") as f:
        param = json.load(f)

    for key, value in param.items():
        param[key] = np.array(value)
    return param

def rejection_sample(pdf, xmin, xmax, size=100, chunk_size=10000):
    """
    Helper function for rejection sampling from a pdf with maximum and minimum arguments.
    
    Parameters
    ----------
    pdf : `function`
        pdf function.
    xmin : `float`
        minimum value of the pdf.
    xmax : `float`
        maximum value of the pdf.
    size : `int`, optional
        number of samples. Default is 100.
    chunk_size : `int`, optional
        chunk size for sampling. Default is 10000.

    Returns
    ----------
    x_sample : `numpy.ndarray`
        samples from the pdf.
    """
    x = np.linspace(xmin, xmax, chunk_size)
    y = pdf(x)
    # Maximum value of the pdf
    ymax = np.max(y)

    # Rejection sample in chunks
    x_sample = []
    while len(x_sample) < size:
        x_try = np.random.uniform(xmin, xmax, size=chunk_size)
        pdf_x_try = pdf(x_try) # Calculate the pdf at the random x values
        # this is for comparing with the pdf value at x_try, will be used to accept or reject the sample
        y_try = np.random.uniform(0, ymax, size=chunk_size)
        
        # Update the maximum value of the pdf
        ymax = max(ymax, np.max(pdf_x_try))  
        
        # applying condition to accept the sample
        # Add while retaining 1D shape of the list
        x_sample += list(x_try[y_try < pdf_x_try])

    # Transform the samples to a 1D numpy array
    x_sample = np.array(x_sample).flatten()
    # Return the correct number of samples
    return x_sample[:size]


def rejection_sample2d(pdf, xmin, xmax, ymin, ymax, size=100, chunk_size=10000):
    """
    Helper function for rejection sampling from a 2D pdf with maximum and minimum arguments.

    Parameters
    ----------
    pdf : `function`
        2D pdf function.
    xmin : `float`
        minimum value of the pdf in the x-axis.
    xmax : `float`
        maximum value of the pdf in the x-axis.
    ymin : `float`
        minimum value of the pdf in the y-axis.
    ymax : `float`
        maximum value of the pdf in the y-axis.
    size : `int`, optional
        number of samples. Default is 100.
    chunk_size : `int`, optional
        chunk size for sampling. Default is 10000.

    Returns
    ----------
    x_sample : `numpy.ndarray`
        samples from the pdf in the x-axis.
    """

    x = np.random.uniform(xmin, xmax, chunk_size)
    y = np.random.uniform(ymin, ymax, chunk_size)
    z = pdf(x, y)
    # Maximum value of the pdf
    zmax = np.max(z)

    # Rejection sample in chunks
    x_sample = []
    y_sample = []
    while len(x_sample) < size:
        x_try = np.random.uniform(xmin, xmax, size=chunk_size)
        y_try = np.random.uniform(ymin, ymax, size=chunk_size)
        pdf_xy_try = pdf(x_try, y_try)
        # this is for comparing with the pdf value at x_try, will be used to accept or reject the sample
        z_try = np.random.uniform(0, zmax, size=chunk_size)
        
        # Update the maximum value of the pdf
        zmax = max(zmax, np.max(pdf_xy_try))

        x_sample += list(x_try[z_try < pdf_xy_try])
        y_sample += list(y_try[z_try < pdf_xy_try])

    # Transform the samples to a 1D numpy array
    x_sample = np.array(x_sample).flatten()
    y_sample = np.array(y_sample).flatten()
    # Return the correct number of samples
    return x_sample[:size], y_sample[:size]


def add_dictionaries_together(dictionary1, dictionary2):
    """
    Adds two dictionaries with the same keys together.
    
    Parameters
    ----------
    dictionary1 : `dict`
        dictionary to be added.
    dictionary2 : `dict`
        dictionary to be added.

    Returns
    ----------
    dictionary : `dict`
        dictionary with added values.
    """
    dictionary = {}
    # Check if either dictionary empty, in which case only return the dictionary with values
    if len(dictionary1) == 0:
        return dictionary2
    elif len(dictionary2) == 0:
        return dictionary1
    # Check if the keys are the same
    if dictionary1.keys() != dictionary2.keys():
        raise ValueError("The dictionaries have different keys.")
    for key in dictionary1.keys():
        value1 = dictionary1[key]
        value2 = dictionary2[key]

        # check if the value is empty
        bool0 = len(value1) == 0 or len(value2) == 0
        # check if the value is an ndarray or a list
        bool1 = isinstance(value1, np.ndarray) and isinstance(value2, np.ndarray)
        bool2 = isinstance(value1, list) and isinstance(value2, list)
        bool3 = isinstance(value1, np.ndarray) and isinstance(value2, list)
        bool4 = isinstance(value1, list) and isinstance(value2, np.ndarray)
        bool4 = bool4 or bool3
        bool5 = isinstance(value1, dict) and isinstance(value2, dict)

        if bool0:
            if len(value1) == 0 and len(value2) == 0:
                dictionary[key] = np.array([])
            elif len(value1) != 0 and len(value2) == 0:
                dictionary[key] = np.array(value1)
            elif len(value1) == 0 and len(value2) != 0:
                dictionary[key] = np.array(value2)
        elif bool1:
            dictionary[key] = np.concatenate((value1, value2))
        elif bool2:
            dictionary[key] = value1 + value2
        elif bool4:
            dictionary[key] = np.concatenate((np.array(value1), np.array(value2)))
        elif bool5:
            dictionary[key] = add_dictionaries_together(
                dictionary1[key], dictionary2[key]
            )
        else:
            raise ValueError(
                "The dictionary contains an item which is neither an ndarray nor a dictionary."
            )
    return dictionary


def trim_dictionary(dictionary, size):
    """
    Filters an event dictionary to only contain the size.
    
    Parameters
    ----------
    dictionary : `dict`
        dictionary to be trimmed.
    size : `int`
        size to trim the dictionary to.

    Returns
    ----------
    dictionary : `dict`
        trimmed dictionary.
    """
    for key in dictionary.keys():
        # Check if the item is an ndarray
        if isinstance(dictionary[key], np.ndarray):
            dictionary[key] = dictionary[key][:size]  # Trim the array
        # Check if the item is a nested dictionary
        elif isinstance(dictionary[key], dict):
            dictionary[key] = trim_dictionary(
                dictionary[key], size
            )  # Trim the nested dictionary
        else:
            raise ValueError(
                "The dictionary contains an item which is neither an ndarray nor a dictionary."
            )
    return dictionary

def create_func_pdf_invcdf(x, y, category="function"):
    """
    Function to create a interpolated function, inverse function or inverse cdf from the input x and y.

    Parameters
    ----------
    x : `numpy.ndarray`
        x values. This has to sorted in ascending order.
    y : `numpy.ndarray`
        y values. Corresponding to the x values.
    category : `str`, optional
        category of the function. Default is "function". Other options are "function_inverse", "pdf" and "inv_cdf".

    Returns
    ----------
    pdf : `pdf function`
        interpolated pdf function.
    inv_pdf : `function inverse`
        interpolated inverse pdf function.
    inv_cdf : `function`
        interpolated inverse cdf.
    """

    idx = np.argwhere(np.isnan(y))
    x = np.delete(x, idx)
    y = np.delete(y, idx)

    # create pdf with interpolation
    pdf_unorm = interp1d(x, y, kind="cubic", fill_value="extrapolate")
    if category == "function":
        return pdf_unorm
    if category == "function_inverse":
        # create inverse function
        return interp1d(y, x, kind="cubic", fill_value="extrapolate")

    min_, max_ = min(x), max(x)
    norm = quad(pdf_unorm, min_, max_)[0]
    y = y / norm
    if category == "pdf" or category is None:
        # normalize the pdf
        pdf = interp1d(x, y, kind="cubic", fill_value="extrapolate")
        return pdf
    # cdf
    cdf_values = cumtrapz(y, x, initial=0)
    idx = np.argwhere(cdf_values > 0)[0][0]
    cdf_values = cdf_values[idx:]
    x = x[idx:]
    inv_cdf = interp1d(cdf_values, x, kind="cubic", fill_value="extrapolate")
    if category == "inv_cdf":
        return inv_cdf
    if category == "all":
        return([pdf, inv_cdf])
    
def create_conditioned_pdf_invcdf(x, conditioned_y, pdf_func, category):
    """
    pdf_func is the function to calculate the pdf of x given y
    x is an array and the output of pdf_func is an array
    y is the condition
    we consider parameter plane of x and y

    Parameters
    ----------
    x : `numpy.ndarray`
        x values.
    conditioned_y : `numpy.ndarray`
        conditioned y values.
    pdf_func : `function`
        function to calculate the pdf of x given y.
    category : `str`, optional
        category of the function. Default is "function". Other options are "function_inverse", "pdf" and "inv_cdf".
    """

    list_ = []
    for y in conditioned_y:
        phi = pdf_func(x,y)
        # append pdf for each y along the x-axis
        list_.append(create_func_pdf_invcdf(x, phi, category=category))
        
    return list_

def create_func(x, y):
    """
    Function to create a spline interpolated function from the input x and y.

    Parameters
    ----------
    x : `numpy.ndarray`
        x values.
    y : `numpy.ndarray`
        y values.

    Returns
    ----------
    c : `numpy.ndarray`
        spline coefficients.
    """

    idx = np.argwhere(np.isnan(y))
    x = np.delete(x, idx)
    y = np.delete(y, idx)
    return CubicSpline(x, y).c, x

def create_func_inv(x, y):
    """
    Function to create a spline interpolated inverse function from the input x and y.

    Parameters
    ----------
    x : `numpy.ndarray`
        x values.
    y : `numpy.ndarray`
        y values.

    Returns
    ----------
    c : `numpy.ndarray`
        spline coefficients.
    """

    idx = np.argwhere(np.isnan(y))
    x = np.delete(x, idx)
    y = np.delete(y, idx)
    return CubicSpline(y, x).c, y

def create_pdf(x, y):
    """
    Function to create a spline interpolated normalized pdf from the input x and y.

    Parameters
    ----------
    x : `numpy.ndarray`
        x values.
    y : `numpy.ndarray`
        y values.   
    
    Returns
    ----------
    c : `numpy.ndarray`
        spline coefficients.
    """
    idx = np.argwhere(np.isnan(y))
    x = np.delete(x, idx)
    y = np.delete(y, idx)
    pdf_unorm = interp1d(x, y, kind="cubic", fill_value="extrapolate")
    min_, max_ = min(x), max(x)
    norm = quad(pdf_unorm, min_, max_)[0]
    y = y / norm
    return CubicSpline(x, y).c, x

def create_inv_cdf_array(x, y):
    """
    Function to create a spline interpolated inverse cdf from the input x and y.

    Parameters
    ----------
    x : `numpy.ndarray`
        x values.
    y : `numpy.ndarray`
        y values.

    Returns
    ----------
    c : `numpy.ndarray`
        spline coefficients.
    """

    idx = np.argwhere(np.isnan(y))
    x = np.delete(x, idx)
    y = np.delete(y, idx)
    cdf_values = cumtrapz(y, x, initial=0)
    cdf_values = cdf_values / cdf_values[-1]
    # to remove duplicate values on x-axis before interpolation
    # idx = np.argwhere(cdf_values > 0)[0][0]
    # cdf_values = cdf_values[idx:]
    # x = x[idx:]
    # cdf_values = np.insert(cdf_values, 0, 0)
    # x = np.insert(x, 0, x[idx-1])
    return np.array([cdf_values, x])

def create_conditioned_pdf(x, conditioned_y, pdf_func):
    """
    Function to create a conditioned pdf from the input x and y.

    Parameters
    ----------
    x : `numpy.ndarray`
        x values.
    conditioned_y : `numpy.ndarray`
        conditioned y values.
    pdf_func : `function`
        function to calculate the pdf of x given y.

    Returns
    ----------
    list_ : `list`
        list of pdfs.
    """
    list_ = []
    for y in conditioned_y:
        phi = pdf_func(x,y)
        list_.append(create_pdf(x, phi))
        
    return np.array(list_)

def create_conditioned_inv_cdf_array(x, conditioned_y, pdf_func):
    """
    Function to create a conditioned inv_cdf from the input x and y.

    Parameters
    ----------
    x : `numpy.ndarray`
        x values.
    conditioned_y : `numpy.ndarray`
        conditioned y values.
    pdf_func : `function`
        function to calculate the pdf of x given y.

    Returns
    ----------
    list_ : `list`
        list of inv_cdfs.
    """

    list_ = []
    for y in conditioned_y:
        phi = pdf_func(x,y)
        list_.append(create_inv_cdf_array(x, phi))
        
    return np.array(list_)

def interpolator_from_pickle(
    param_dict_given, directory, sub_directory, name, x, pdf_func=None, y=None, conditioned_y=None, dimension=1,category="pdf", create_new=False
):
    """
    Function to decide which interpolator to use.

    Parameters
    ----------
    param_dict_given : `dict`
        dictionary of parameters.
    directory : `str`
        directory to store the interpolator.
    sub_directory : `str`
        sub-directory to store the interpolator.
    name : `str`
        name of the interpolator.
    x : `numpy.ndarray`
        x values.
    pdf_func : `function`
        function to calculate the pdf of x given y.
    y : `numpy.ndarray`
        y values.
    conditioned_y : `numpy.ndarray`
        conditioned y values.
    dimension : `int`
        dimension of the interpolator. Default is 1.
    category : `str`
        category of the function. Default is "pdf".
    create_new : `bool`
        if True, create a new interpolator. Default is False.

    Returns
    ----------
    interpolator : `function`
        interpolator function.
    """

    # check first whether the directory, subdirectory and pickle exist
    path_inv_cdf, it_exist = interpolator_pickle_path(
        param_dict_given=param_dict_given,
        directory=directory,
        sub_directory=sub_directory,
        interpolator_name=name,
    )
    if create_new:
        it_exist = False
    if it_exist:
        print(f"{name} interpolator will be loaded from {path_inv_cdf}")
        # load the interpolator
        with open(path_inv_cdf, "rb") as handle:
            interpolator = pickle.load(handle)
        return interpolator
    else:
        print(f"{name} interpolator will be generated at {path_inv_cdf}")

        # create the interpolator
        if dimension==1:
            if y is None:
                y = pdf_func(x)
            if category=="function":
                interpolator = create_func(x, y)
            elif category=="function_inverse":
                interpolator = create_func_inv(x, y)
            elif category=="pdf":
                interpolator = create_pdf(x, y)
            elif category=="inv_cdf":
                interpolator = create_inv_cdf_array(x, y)
            else:
                raise ValueError("The category given should be function, function_inverse, pdf or inv_cdf.")
        elif dimension==2:
            if category=="pdf":
                interpolator = create_conditioned_pdf(x, conditioned_y, pdf_func)
            elif category=="inv_cdf":
                interpolator = create_conditioned_inv_cdf_array(x, conditioned_y, pdf_func)
        else:
            raise ValueError("The dimension is not supported.")
        # save the interpolator
        with open(path_inv_cdf, "wb") as handle:
            pickle.dump(interpolator, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return interpolator

def interpolator_pickle_path(
    param_dict_given,
    directory,
    sub_directory,
    interpolator_name,
):
    """
    Function to create the interpolator pickle file path.

    Parameters
    ----------
    param_dict_given : `dict`
        dictionary of parameters.
    directory : `str`
        directory to store the interpolator.
    sub_directory : `str`
        sub-directory to store the interpolator.
    interpolator_name : `str`
        name of the interpolator.

    Returns
    ----------
    path_inv_cdf : `str`
        path of the interpolator pickle file.
    it_exist : `bool`
        if True, the interpolator exists.
    """

    # check the dir 'interpolator' exist
    full_dir = directory + "/" + sub_directory
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.makedirs(full_dir)
    else:
        if not os.path.exists(full_dir):
            os.makedirs(full_dir)

    # check if param_dict_list.pickle exists
    path1 = full_dir + "/init_dict.pickle"
    if not os.path.exists(path1):
        dict_list = []
        with open(path1, "wb") as handle:
            pickle.dump(dict_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # check if the input dict is the same as one of the dict inside the pickle file
    param_dict_stored = pickle.load(open(path1, "rb"))

    path2 = full_dir
    len_ = len(param_dict_stored)
    if param_dict_given in param_dict_stored:
        idx = param_dict_stored.index(param_dict_given)
        # load the interpolator
        path_inv_cdf = path2 + "/" + interpolator_name + "_" + str(idx) + ".pickle"
        # there will be exception if the file is deleted by mistake
        if os.path.exists(path_inv_cdf):
            it_exist = True
        else:
            it_exist = False
    else:
        it_exist = False
        path_inv_cdf = path2 + "/" + interpolator_name + "_" + str(len_) + ".pickle"
        param_dict_stored.append(param_dict_given)
        with open(path1, "wb") as handle:
            pickle.dump(param_dict_stored, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return path_inv_cdf, it_exist

def interpolator_pdf_conditioned(x, conditioned_y, y_array, interpolator_list):
    """
    Function to find the pdf interpolator coefficients from the conditioned y.

    Parameters
    ----------
    x : `numpy.ndarray`
        x values.
    conditioned_y : `float`
        conditioned y value.
    y_array : `numpy.ndarray`
        y values.
    interpolator_list : `list`
        list of interpolators.

    Returns
    ----------
    interpolator_list[idx](x) : `numpy.ndarray`
        samples from the interpolator.
    """
    # find the index of z in zlist
    idx = np.searchsorted(y_array, conditioned_y)

    return interpolator_list[idx](x)

def interpolator_sampler_conditioned(conditioned_y, y_array, interpolator_list, size=1000):
    """
    Function to find sampler interpolator coefficients from the conditioned y.

    Parameters
    ----------
    conditioned_y : `float`
        conditioned y value.
    y_array : `numpy.ndarray`
        y values.
    interpolator_list : `list`
        list of interpolators.
    size : `int`
        number of samples.

    Returns
    ----------
    """

    # find the index of z in zlist
    idx = np.searchsorted(y_array, conditioned_y)
    u = np.random.uniform(0, 1, size=size)
    return interpolator_list[idx](u)

@njit
def cubic_spline_interpolator(xnew, coefficients, x):
    """
    Function to interpolate using cubic spline.

    Parameters
    ----------
    xnew : `numpy.ndarray`
        new x values.
    coefficients : `numpy.ndarray`
        coefficients of the cubic spline.
    x : `numpy.ndarray`
        x values.

    Returns
    ----------
    result : `numpy.ndarray`
        interpolated values.
    """

    # Handling extrapolation
    i = np.searchsorted(x, xnew) - 1
    idx1 = xnew <= x[0]
    idx2 = xnew > x[-1]
    i[idx1] = 0
    i[idx2] = len(x) - 2

    # Calculate the relative position within the interval
    dx = xnew - x[i]

    # Calculate the interpolated value
    # Cubic polynomial: a + b*dx + c*dx^2 + d*dx^3
    a, b, c, d = coefficients[:, i]
    #result = a + b*dx + c*dx**2 + d*dx**3
    result = d + c*dx + b*dx**2 + a*dx**3
    return result

@njit
def inverse_transform_sampler(size, cdf, x):
    """
    Function to sample from the inverse transform method.

    Parameters
    ----------
    size : `int`
        number of samples.
    cdf : `numpy.ndarray`
        cdf values.
    x : `numpy.ndarray`
        x values.

    Returns
    ----------
    samples : `numpy.ndarray`
        samples from the cdf.
    """

    u = np.random.uniform(0, 1, size)
    idx = np.searchsorted(cdf, u)
    x1, x0, y1, y0 = cdf[idx], cdf[idx-1], x[idx], x[idx-1]
    samples = y0 + (y1 - y0) * (u - x0) / (x1 - x0)
    return samples

def batch_handler(size, batch_size, sampling_routine, output_jsonfile, save_batch=True, resume=False, param_name='parameters'):
    """
    Function to run the sampling in batches.

    Parameters
    ----------
    size : `int`
        number of samples.
    batch_size : `int`
        batch size.
    sampling_routine : `function`
        sampling function. It should have 'size' as input and return a dictionary.
    output_jsonfile : `str`
        json file name for storing the parameters.
    save_batch : `bool`, optional
        if True, save sampled parameters in each iteration. Default is True.
    resume : `bool`, optional
        if True, resume sampling from the last batch. Default is False.
    param_name : `str`, optional
        name of the parameter. Default is 'parameters'.

    Returns
    ----------
    dict_buffer : `dict`
        dictionary of parameters.
    """

    # sampling in batches
    if resume and os.path.exists(output_jsonfile):
        # get sample from json file
        dict_buffer = get_param_from_json(output_jsonfile)
    else:
        dict_buffer = None
    
    # if size is multiple of batch_size
    if size % batch_size == 0:
        num_batches = size // batch_size
    # if size is not multiple of batch_size
    else:
        num_batches = size // batch_size + 1

    print(
        f"chosen batch size = {batch_size} with total size = {size}"
    )
    print(f"There will be {num_batches} batche(s)")

    # note frac_batches+(num_batches-1)*batch_size = size
    if size > batch_size:
        frac_batches = size - (num_batches - 1) * batch_size
    # if size is less than batch_size
    else:
        frac_batches = size
    track_batches = 0  # to track the number of batches

    if not resume:
        # create new first batch with the frac_batches
        track_batches, dict_buffer = create_batch_params(sampling_routine, frac_batches, dict_buffer, save_batch, output_jsonfile, track_batches=track_batches)
    else:
        # check where to resume from
        # identify the last batch and assign current batch number
        # try-except is added to avoid the error when the file does not exist or if the file is empty or corrupted or does not have the required key.
        try:
            print(f"resuming from {output_jsonfile}")
            len_ = len(list(dict_buffer.values())[0])
            track_batches = (len_ - frac_batches) // batch_size + 1
        except:
            # create new first batch with the frac_batches
            track_batches, dict_buffer = create_batch_params(sampling_routine, frac_batches, dict_buffer, save_batch, output_jsonfile, track_batches=track_batches)

    # loop over the remaining batches
    min_, max_ = track_batches, num_batches
    # print(f"min_ = {min_}, max_ = {max_}")
    save_param = False
    if min_ == max_:
        print(f"{param_name} already sampled.")
    elif min_ > max_:
        len_ = len(list(dict_buffer.values())[0])
        print(f"existing {param_name} size is {len_} is more than the required size={size}. It will be trimmed.")
        dict_buffer = trim_dictionary(dict_buffer, size)
        save_param = True
    else:
        for i in range(min_, max_):
            _, dict_buffer = create_batch_params(sampling_routine, batch_size, dict_buffer, save_batch, output_jsonfile, track_batches=i, resume=True)

        if save_batch:
            # if save_batch=True, then dict_buffer is only the last batch
            dict_buffer = get_param_from_json(output_jsonfile)
        else:  # dont save in batches
            # this if condition is required if there is nothing to save
            save_param = True
    
    if save_param:
        # store all params in json file
        print(f"saving all {param_name} in {output_jsonfile} ")
        append_json(output_jsonfile, dict_buffer, replace=True)

    return dict_buffer

def create_batch_params(sampling_routine, frac_batches, dict_buffer, save_batch, output_jsonfile, track_batches, resume=False):
    """
    Helper function to batch_handler. It create batch parameters and store in a dictionary.

    Parameters
    ----------
    sampling_routine : `function`
        sampling function. It should have 'size' as input and return a dictionary.
    frac_batches : `int`
        batch size.
    dict_buffer : `dict`
        dictionary of parameters.
    save_batch : `bool`
        if True, save sampled parameters in each iteration.
    output_jsonfile : `str`
        json file name for storing the parameters.
    track_batches : `int`
        track the number of batches.
    resume : `bool`, optional
        if True, resume sampling from the last batch. Default is False.

    Returns
    ----------
    track_batches : `int`
        track the number of batches.
    """

    track_batches = track_batches + 1
    print(f"Batch no. {track_batches}")
    param = sampling_routine(size=frac_batches, save_batch=save_batch, output_jsonfile=output_jsonfile, resume=resume)

    # adding batches and hold it in the buffer attribute.
    if not save_batch:
        # in the new batch (new sampling run), dict_buffer will be None
        if dict_buffer is None:
            dict_buffer = param
        else:
            for key, value in param.items():
                dict_buffer[key] = np.concatenate((dict_buffer[key], value))
    else:
        # store all params in json file
        dict_buffer = append_json(file_name=output_jsonfile, new_dictionary=param,  old_dictionary=dict_buffer, replace=not (resume))

    return track_batches, dict_buffer