# -*- coding: utf-8 -*-
""" 
This module contains helper routines for other modules in the ler package.
"""

import os
import pickle
import h5py
import numpy as np
import json
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from numba import njit
from importlib import resources
# import datetime
from numba.core.registry import CPUDispatcher


def is_njitted(func):
    return isinstance(func, CPUDispatcher)

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

def load_pickle(file_name):
    """Load a pickle file.

    Parameters
    ----------
    file_name : `str`
        pickle file name for storing the parameters.

    Returns
    ----------
    param : `dict`
    """
    with open(file_name, "rb") as handle:
        param = pickle.load(handle)

    return param

def save_pickle(file_name, param):
    """Save a dictionary as a pickle file.

    Parameters
    ----------
    file_name : `str`
        pickle file name for storing the parameters.
    param : `dict`
        dictionary to be saved as a pickle file.
    """
    with open(file_name, "wb") as handle:
        pickle.dump(param, handle, protocol=pickle.HIGHEST_PROTOCOL)

# hdf5
def load_hdf5(file_name):
    """Load a hdf5 file.

    Parameters
    ----------
    file_name : `str`
        hdf5 file name for storing the parameters.

    Returns
    ----------
    param : `dict`
    """

    return h5py.File(file_name, 'r')

def save_hdf5(file_name, param):
    """Save a dictionary as a hdf5 file.

    Parameters
    ----------
    file_name : `str`
        hdf5 file name for storing the parameters.
    param : `dict`
        dictionary to be saved as a hdf5 file.
    """
    with h5py.File(file_name, 'w') as f:
        for key, value in param.items():
            f.create_dataset(key, data=value)

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

    Example
    ----------
    >>> from ler.utils import save_json, load_json
    >>> param = dict_list = [{"name": "lens_mass", "type": "cubic_spline", "x": np.array([0, 1, 2, 3]), "y": [0, 1, 0, 1]},]
    >>> save_json("param.json", param)
    >>> param = load_json("param.json")
    """
    try:
        with open(file_name, "w", encoding="utf-8") as write_file:
            json.dump(param, write_file)
    except:
        with open(file_name, "w", encoding="utf-8") as write_file:
            json.dump(param, write_file, indent=4, cls=NumpyEncoder)

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

def concatenate_dict_values(dict1, dict2):
    """Adds the values of two dictionaries together.
    
    Parameters
    ----------
    dict1 : `dict`
        dictionary to be added.
    dict2 : `dict`
        dictionary to be added.

    Returns
    ----------
    dict1 : `dict`
        dictionary with added values.
    """
    data_key = dict1.keys()
    for key, value in dict2.items():
        if key in data_key:
            dict1[key] = np.concatenate((dict1[key], value))

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

def load_txt_from_module(package, directory, filename):
    """
    Function to load a specific dataset from a .txt file within the package

    Parameters
    ----------
    package : str
        name of the package
    directory : str
        name of the directory within the package
    filename : str
        name of the .txt file

    Returns
    ----------
    data : `dict`
        Dictionary loaded from the .txt file
    """

    with resources.path(package + '.' + directory, filename) as txt_path:
        return np.loadtxt(txt_path)
        
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
    cdf_values = cumulative_trapezoid(y, x, initial=0)
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
    cdf_values = cumulative_trapezoid(y, x, initial=0)
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
    size = len(x)
    for y in conditioned_y:
        phi = pdf_func(x,y*np.ones(size))
        list_.append(create_inv_cdf_array(x, phi))
        
    return np.array(list_)

def interpolator_from_json(
    identifier_dict, directory, sub_directory, name, x, pdf_func=None, y=None, conditioned_y=None, dimension=1,category="pdf", create_new=False
):
    """
    Function to decide which interpolator to use.

    Parameters
    ----------
    identifier_dict : `dict`
        dictionary of identifiers.
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

    # check first whether the directory, subdirectory and json exist
    path_inv_cdf, it_exist = interpolator_json_path(
        identifier_dict=identifier_dict,
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
            interpolator = json.load(handle)
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
        save_json(path_inv_cdf, interpolator)
        return interpolator

def interpolator_json_path(
    identifier_dict,
    directory,
    sub_directory,
    interpolator_name,
):
    """
    Function to create the interpolator json file path.

    Parameters
    ----------
    identifier_dict : `dict`
        dictionary of identifiers.
    directory : `str`
        directory to store the interpolator.
    sub_directory : `str`
        sub-directory to store the interpolator.
    interpolator_name : `str`
        name of the interpolator.

    Returns
    ----------
    path_inv_cdf : `str`
        path of the interpolator json file.
    it_exist : `bool`
        if True, the interpolator exists.
    """

    # convert values of identifier_dict to string
    identifier_dict = {k: str(v) for k, v in identifier_dict.items()}

    # check the dir 'interpolator' exist
    full_dir = directory + "/" + sub_directory
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.makedirs(full_dir)
    else:
        if not os.path.exists(full_dir):
            os.makedirs(full_dir)

    # check if identifier_dict_list.json exists
    path1 = full_dir + "/init_dict.json"
    if not os.path.exists(path1):
        dict_list = []
        save_json(path1, dict_list)

    # check if the input dict is the same as one of the dict inside the json file
    identifier_dict_stored = load_json(path1)

    path2 = full_dir
    len_ = len(identifier_dict_stored)
    if identifier_dict in identifier_dict_stored:
        idx = identifier_dict_stored.index(identifier_dict)
        # load the interpolator
        path_inv_cdf = path2 + "/" + interpolator_name + "_" + str(idx) + ".json"
        # there will be exception if the file is deleted by mistake
        if os.path.exists(path_inv_cdf):
            it_exist = True
        else:
            it_exist = False
    else:
        it_exist = False
        path_inv_cdf = path2 + "/" + interpolator_name + "_" + str(len_) + ".json"
        identifier_dict_stored.append(identifier_dict)
        save_json(path1, identifier_dict_stored)

    return path_inv_cdf, it_exist

# def interpolator_pdf_conditioned(x, conditioned_y, y_array, interpolator_list):
#     """
#     Function to find the pdf interpolator coefficients from the conditioned y.

#     Parameters
#     ----------
#     x : `numpy.ndarray`
#         x values.
#     conditioned_y : `float`
#         conditioned y value.
#     y_array : `numpy.ndarray`
#         y values.
#     interpolator_list : `list`
#         list of interpolators.

#     Returns
#     ----------
#     interpolator_list[idx](x) : `numpy.ndarray`
#         samples from the interpolator.
#     """
#     # find the index of z in zlist
#     idx = np.searchsorted(y_array, conditioned_y)

#     return interpolator_list[idx](x)

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

def monte_carlo_integration(function, uniform_prior, size=10000):
    """
    Function to perform Monte Carlo integration.

    Parameters
    ----------
    function : `function`
        function to be integrated.
    prior : `function`
        prior function.

    Returns
    ----------
    integral : `float`
        integral value.
    """

    # sample from the prior
    x = uniform_prior(size=size)
    # calculate the function
    y = np.array([function(x[i]) for i in range(size)])
    # calculate the integral
    integral = np.mean(y)*(np.max(x)-np.min(x))
    return integral

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
    #i = np.array(i, dtype=np.int32)
    # x = np.array(x, dtype=np.float64)

    # Calculate the relative position within the interval
    #print(f"i = {i}\n xnew = {xnew}\n x = {x}")
    #print(x[i])
    dx = xnew - x[i]

    # Calculate the interpolated value
    # Cubic polynomial: a + b*dx + c*dx^2 + d*dx^3
    # coefficients = np.array(coefficients, dtype=np.float64)
    #print(f"coefficients = {coefficients}")
    a, b, c, d = coefficients[:, i]
    #result = a + b*dx + c*dx**2 + d*dx**3
    result = d + c*dx + b*dx**2 + a*dx**3
    return result

@njit
def pdf_cubic_spline_interpolator(xnew, norm_const, coefficients, x):
    """
    Function to interpolate pdf using cubic spline.

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
    #i = np.array(i, dtype=np.int32)
    # x = np.array(x, dtype=np.float64)

    # Calculate the relative position within the interval
    #print(f"i = {i}\n xnew = {xnew}\n x = {x}")
    #print(x[i])
    dx = xnew - x[i]

    # Calculate the interpolated value
    # Cubic polynomial: a + b*dx + c*dx^2 + d*dx^3
    # coefficients = np.array(coefficients, dtype=np.float64)
    #print(f"coefficients = {coefficients}")
    a, b, c, d = coefficients[:, i]
    #result = a + b*dx + c*dx**2 + d*dx**3
    result = d + c*dx + b*dx**2 + a*dx**3
    return result/norm_const

@njit
def pdf_cubic_spline_interpolator2d_array(xnew_array, ynew_array, norm_array, coefficients, x, y):
    """
    Function to calculate the interpolated value of snr_partialscaled given the mass ratio (ynew) and total mass (xnew). This is based off 2D bicubic spline interpolation.

    Parameters
    ----------
    xnew_array : `numpy.ndarray`
        Total mass of the binary.
    ynew_array : `numpy.ndarray`
        Mass ratio of the binary.
    coefficients : `numpy.ndarray`
        Array of coefficients for the cubic spline interpolation.
    x : `numpy.ndarray`
        Array of total mass values for the coefficients.
    y : `numpy.ndarray`
        Array of mass ratio values for the coefficients.

    Returns
    -------
    result : `float`
        Interpolated value of snr_partialscaled.
    """
    result_array = []
    for i in range(len(xnew_array)):
        xnew = xnew_array[i]
        ynew = ynew_array[i]

        len_y = len(y)
        # find the index nearest to the ynew in y
        y_idx = np.searchsorted(y, ynew) - 1 if ynew > y[0] else 0

        if (ynew>y[0]) and (ynew<y[1]):
            if ynew > y[y_idx] + (y[y_idx+1] - y[y_idx]) / 2:
                y_idx = y_idx + 1
            result=pdf_cubic_spline_interpolator(np.array([xnew]), norm_array[y_idx], coefficients[y_idx], x[y_idx])[0]
        elif y_idx == 0:  # lower end point
            result=pdf_cubic_spline_interpolator(np.array([xnew]), norm_array[0], coefficients[0], x[0])[0]
            # print("a")
        elif y_idx+1 == len_y:  # upper end point
            result=pdf_cubic_spline_interpolator(np.array([xnew]), norm_array[-1], coefficients[-1], x[-1])[0]
            # print("b")
        elif y_idx+2 == len_y:  # upper end point
            result=pdf_cubic_spline_interpolator(np.array([xnew]), norm_array[-1], coefficients[-1], x[-1])[0]
            # print("b")
        else:
            y_idx1 = y_idx - 1
            y_idx2 = y_idx
            y_idx3 = y_idx + 1
            y_idx4 = y_idx + 2
            coeff_low, coeff_high = 4, 8
            # print("c")
            y1, y2, y3, y4 = y[y_idx1], y[y_idx2], y[y_idx3], y[y_idx4]
            z1 = pdf_cubic_spline_interpolator(np.array([xnew]), norm_array[y_idx1], coefficients[y_idx1], x[y_idx1])[0]
            z2 = pdf_cubic_spline_interpolator(np.array([xnew]), norm_array[y_idx2], coefficients[y_idx2], x[y_idx2])[0]
            z3 = pdf_cubic_spline_interpolator(np.array([xnew]), norm_array[y_idx3], coefficients[y_idx3], x[y_idx3])[0]
            z4 = pdf_cubic_spline_interpolator(np.array([xnew]), norm_array[y_idx4], coefficients[y_idx4], x[y_idx4])[0]

            coeff = coefficients_generator_ler(y1, y2, y3, y4, z1, z2, z3, z4)
            matrixD = coeff[coeff_low:coeff_high]
            matrixB = np.array([ynew**3, ynew**2, ynew, 1])
            result = np.dot(matrixB, matrixD) 

        result_array.append(result)

    result_array = np.array(result_array)
    idx = result_array < 0.0
    result_array[idx] = 0.0

    return result_array

@njit
def cubic_spline_interpolator2d_array(xnew_array, ynew_array, coefficients, x, y):
    """
    Function to calculate the interpolated value of snr_partialscaled given the mass ratio (ynew) and total mass (xnew). This is based off 2D bicubic spline interpolation.

    Parameters
    ----------
    xnew_array : `numpy.ndarray`
        Total mass of the binary.
    ynew_array : `numpy.ndarray`
        Mass ratio of the binary.
    coefficients : `numpy.ndarray`
        Array of coefficients for the cubic spline interpolation.
    x : `numpy.ndarray`
        Array of total mass values for the coefficients.
    y : `numpy.ndarray`
        Array of mass ratio values for the coefficients.

    Returns
    -------
    result : `float`
        Interpolated value of snr_partialscaled.
    """
    result_array = []
    for i in range(len(xnew_array)):
        xnew = xnew_array[i]
        ynew = ynew_array[i]

        len_y = len(y)
        # find the index nearest to the ynew in y
        y_idx = np.searchsorted(y, ynew) - 1 if ynew > y[0] else 0

        if (ynew>y[0]) and (ynew<y[1]):
            if ynew > y[y_idx] + (y[y_idx+1] - y[y_idx]) / 2:
                y_idx = y_idx + 1
            result_array.append(cubic_spline_interpolator(np.array([xnew]), coefficients[y_idx], x[y_idx])[0])
            # print(f"a) idx = {y_idx}")
        elif y_idx == 0:  # lower end point
            result_array.append(cubic_spline_interpolator(np.array([xnew]), coefficients[0], x[0])[0])
            # print(f"a) idx = {y_idx}")
        elif y_idx+1 == len_y:  # upper end point
            result_array.append(cubic_spline_interpolator(np.array([xnew]), coefficients[-1], x[-1])[0])
            # print(f"b) idx = {y_idx}")
        elif y_idx+2 == len_y:  # upper end point
            result_array.append(cubic_spline_interpolator(np.array([xnew]), coefficients[-1], x[-1])[0])
            # print(f"c) idx = {y_idx}")
        else:
            y_idx1 = y_idx - 1
            y_idx2 = y_idx
            y_idx3 = y_idx + 1
            y_idx4 = y_idx + 2
            coeff_low, coeff_high = 4, 8
            # print(f"d) idx1, idx2, idx3, idx4 = {y_idx1}, {y_idx2}, {y_idx3}, {y_idx4}")
            y1, y2, y3, y4 = y[y_idx1], y[y_idx2], y[y_idx3], y[y_idx4]
            z1 = cubic_spline_interpolator(np.array([xnew]), coefficients[y_idx1], x[y_idx1])[0]
            z2 = cubic_spline_interpolator(np.array([xnew]), coefficients[y_idx2], x[y_idx2])[0]
            z3 = cubic_spline_interpolator(np.array([xnew]), coefficients[y_idx3], x[y_idx3])[0]
            z4 = cubic_spline_interpolator(np.array([xnew]), coefficients[y_idx4], x[y_idx4])[0]

            coeff = coefficients_generator_ler(y1, y2, y3, y4, z1, z2, z3, z4)
            matrixD = coeff[coeff_low:coeff_high]
            matrixB = np.array([ynew**3, ynew**2, ynew, 1])
            result_array.append(np.dot(matrixB, matrixD))

    return np.array(result_array)

@njit
def coefficients_generator_ler(y1, y2, y3, y4, z1, z2, z3, z4):
    """
    Function to generate the coefficients for the cubic spline interpolation of fn(y)=z.

    Parameters
    ----------
    y1, y2, y3, y4, z1, z2, z3, z4: `float`
        Values of y and z for the cubic spline interpolation.

    Returns
    -------
    coefficients: `numpy.ndarray`
        Coefficients for the cubic spline interpolation.
    """
    matrixA = np.array([
        [y1**3, y1**2, y1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [y2**3, y2**2, y2, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, y2**3, y2**2, y2, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, y3**3, y3**2, y3, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, y3**3, y3**2, y3, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, y4**3, y4**2, y4, 1],
        [3*y2**2, 2*y2, 1, 0, -3*y2**2, -2*y2, -1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 3*y3**2, 2*y3, 1, 0, -3*y3**2, -2*y3, -1, 0],
        [6*y2, 2, 0, 0, -6*y2, -2, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 6*y3, 2, 0, 0, -6*y3, -2, 0, 0],
        [6*y1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 6*y4, 2, 0, 0],
    ])
    matrixC = np.array([z1, z2, z2, z3, z3, z4, 0, 0, 0, 0, 0, 0])
    return np.dot(np.linalg.inv(matrixA), matrixC)

@njit
def inverse_transform_sampler2d(size, conditioned_y, cdf2d, x2d, y1d):
    """
    Function to find sampler interpolator coefficients from the conditioned y.

    Parameters
    ----------
    size: `int`
        Size of the sample.
    conditioned_y: `float`
        Conditioned y value.
    cdf2d: `numpy.ndarray`
        2D array of cdf values.
    x2d: `numpy.ndarray`
        2D array of x values.
    y1d: `numpy.ndarray`
        1D array of y values.

    Returns
    -------
    samples: `numpy.ndarray`
        Samples of the conditioned y.
    """

    samples = np.zeros(size)
    for i, y in enumerate(conditioned_y):
        # find the index nearest to the conditioned_y in y1d
        y_idx = np.searchsorted(y1d, y) - 1 if y > y1d[0] else 0
        if y > y1d[y_idx] + (y1d[y_idx+1] - y1d[y_idx]) / 2:
            y_idx = y_idx + 1

        # linear scaling
        # new cdf values 
        x1, x0, y1, y0 = y1d[y_idx], y1d[y_idx-1], cdf2d[y_idx], cdf2d[y_idx-1]
        cdf = y0 + (y1 - y0) * (y - x0) / (x1 - x0)
        # new x values
        # x=x2d[0], if all rows are same 
        x1, x0, y1, y0 = y1d[y_idx], y1d[y_idx-1], x2d[y_idx], x2d[y_idx-1]
        x = y0 + (y1 - y0) * (y - x0) / (x1 - x0)

        while True:
            u = np.random.uniform(0, 1)
            idx = np.searchsorted(cdf, u)
            x1, x0, y1, y0 = cdf[idx], cdf[idx-1], x[idx], x[idx-1]
            samples[i] = y0 + (y1 - y0) * (u - x0) / (x1 - x0)
            if not np.isnan(samples[i]):
                break
    return samples


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
    # Clip idx to valid range to prevent boundary issues:
    # - When idx=0 (u <= cdf[0]), idx-1 would wrap to -1 (last element)
    # - When idx=len(cdf) (u > cdf[-1]), idx would be out of bounds
    idx = np.clip(idx, 1, len(cdf) - 1)
    x1, x0, y1, y0 = cdf[idx], cdf[idx-1], x[idx], x[idx-1]
    samples = y0 + (y1 - y0) * (u - x0) / (x1 - x0)
    return samples

@njit
def normal_pdf(x, mean=0., std=0.05):
    """
    Calculate the value of a normal probability density function.

    Parameters
    ----------
    x : `float` or `numpy.ndarray`
        The value(s) at which to evaluate the PDF.
    mean : `float`, optional
        The mean of the normal distribution. Default is 0.
    std : `float`, optional
        The standard deviation of the normal distribution. Default is 0.05.

    Returns
    -------
    pdf : `float` or `numpy.ndarray`
        The probability density function value(s) at x.
    """

    return (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std)**2)


@njit
def normal_pdf_2d(x, y, mean_x=0., std_x=0.05, mean_y=0., std_y=0.05):
    """
    Calculate the value of a 2D normal probability density function.

    Parameters
    ----------
    x : `float`
        The x-coordinate for which the PDF is evaluated.
    y : `float`
        The y-coordinate for which the PDF is evaluated.
    mean_x : `float`, optional
        The mean of the normal distribution along the x-axis. Default is 0.
    std_x : `float`, optional
        The standard deviation of the normal distribution along the x-axis. Default is 0.05.
    mean_y : `float`, optional
        The mean of the normal distribution along the y-axis. Default is 0.
    std_y : `float`, optional
        The standard deviation of the normal distribution along the y-axis. Default is 0.05.

    Returns
    -------
    `float`
        The probability density function value at the given x and y coordinates.
    """

    factor_x = (1 / (std_x * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean_x) / std_x)**2)
    factor_y = (1 / (std_y * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((y - mean_y) / std_y)**2)
    return factor_x * factor_y


def load_txt_from_module(package, directory, filename):
    """
    """

    with resources.path(package + '.' + directory, filename) as txt_path:
        return np.loadtxt(txt_path)
    
@njit
def cumulative_trapezoid(y, x=None, dx=1.0, initial=0.0):
    """
    Compute the cumulative integral of a function using the trapezoidal rule.
    """
    if x is None:
        x = np.arange(len(y)) * dx

    # Calculate the cumulative integral using trapezoidal rule
    cumsum = np.zeros_like(y)
    cumsum[0] = initial
    for i in range(1, len(y)):
        cumsum[i] = cumsum[i - 1] + (y[i - 1] + y[i]) * (x[i] - x[i - 1]) / 2.0

    return cumsum
    
@njit
def sample_from_powerlaw_distribution(size, alphans, mminns, mmaxns):
    """
    Inverse transform sampling for a power-law mass distribution:
    p(m) ∝ m^{-alphans}, m in [mminns, mmaxns]
    
    Parameters
    ----------
    size : int
        Number of samples to generate.
    alphans : float
        Power-law index (alpha).
    mminns : float
        Minimum neutron star mass (lower bound).
    mmaxns : float
        Maximum neutron star mass (upper bound).
    random_state : int, np.random.Generator, or None
        Seed or random generator for reproducibility.
    
    Returns
    -------
    m : ndarray
        Array of sampled neutron star masses.
    """

    u = np.random.uniform(0, 1, size)

    if alphans == 1.0:
        # Special case α=1
        m = mminns * (mmaxns / mminns) ** u
    elif alphans == 0.0:
        # Special case α=0 (uniform distribution)
        m = mminns + (mmaxns - mminns) * u
    else:
        pow1 = 1.0 - alphans
        mmin_pow = mminns ** pow1
        mmax_pow = mmaxns ** pow1
        m = (u * (mmax_pow - mmin_pow) + mmin_pow) ** (1.0 / pow1)

    return m
    
