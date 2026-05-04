# -*- coding: utf-8 -*-
"""
This module contains helper routines for other modules in the ler package.
"""

import os
import h5py
import numpy as np
import json
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from numba import njit, prange
from importlib import resources
# import datetime
from numba.core.registry import CPUDispatcher


def is_njitted(func):
    return isinstance(func, CPUDispatcher)

# ---------------------
# IO bound functions
# ---------------------
def remove_file(file_name):
    """Remove a file."""
    print(f"removing {file_name} if it exists")
    if os.path.exists(file_name):
        os.remove(file_name)

# pickle
def load_pickle(file_name):
    """Load a pickle file.

    Parameters
    ----------
    file_name : `str`
        pickle file name for storing the parameters.

    Returns
    ----------
    param : `dict`

    Note
    ----
    Uses dill instead of pickle to support serialization of local functions,
    lambdas, and closures.
    """
    import dill

    with open(file_name, "rb") as handle:
        param = dill.load(handle)

    return param

def save_pickle(file_name, param):
    """Save a dictionary as a pickle file.

    Parameters
    ----------
    file_name : `str`
        pickle file name for storing the parameters.
    param : `dict`
        dictionary to be saved as a pickle file.

    Note
    ----
    Uses dill instead of pickle to support serialization of local functions,
    lambdas, and closures.
    """
    import dill

    with open(file_name, "wb") as handle:
        dill.dump(param, handle, protocol=dill.HIGHEST_PROTOCOL)

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

    return h5py.File(file_name, "r")

def save_hdf5(file_name, param):
    """Save a dictionary as a hdf5 file.

    Parameters
    ----------
    file_name : `str`
        hdf5 file name for storing the parameters.
    param : `dict`
        dictionary to be saved as a hdf5 file.
    """
    with h5py.File(file_name, "w") as f:
        for key, value in param.items():
            f.create_dataset(key, data=value)

# json
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
    >>> from ler.utils import NumpyEncoder
    >>> # create a dictionary
    >>> param = {'a': np.array([1,2,3]), 'b': np.array([4,5,6])}
    >>> # save the dictionary as json file
    >>> with open('param.json', 'w') as f:
    >>>     json.dump(param, f, cls=NumpyEncoder)
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
    except TypeError:
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
        # print(f" {file_name} file does not exist. Creating a new one...")
        replace = True
        data = new_dictionary
    else:
        # print("getting data from file")
        with open(file_name, "r", encoding="utf-8") as f:
            data = json.load(f)
        # json.load returns lists; convert to numpy arrays for compatibility
        for key, value in data.items():
            if isinstance(value, list):
                data[key] = np.array(value)
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
    # print(data)
    with open(file_name, "w", encoding="utf-8") as write_file:
        json.dump(data, write_file, indent=4, cls=NumpyEncoder)
    # end = datetime.datetime.now()
    # print(f"Time taken to save the json file: {end-start}")

    return data

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

# # dictionary handling
# def concatenate_dict_values(dict1, dict2):
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

#     return dict1

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
        len1 = value1.shape[0]
        len2 = value2.shape[0]
        bool0 = len1 == 0 or len2 == 0
        # check if the value is an ndarray or a list
        bool1 = isinstance(value1, np.ndarray) and isinstance(value2, np.ndarray)
        bool2 = isinstance(value1, list) and isinstance(value2, list)
        bool3 = isinstance(value1, np.ndarray) and isinstance(value2, list)
        bool4 = isinstance(value1, list) and isinstance(value2, np.ndarray)
        bool4 = bool4 or bool3
        bool5 = isinstance(value1, dict) and isinstance(value2, dict)

        if bool0:
            if len1 == 0 and len2 == 0:
                dictionary[key] = np.array([])
            elif len1 != 0 and len2 == 0:
                dictionary[key] = np.array(value1)
            elif len1 == 0 and len2 != 0:
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

# module loading
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

    with resources.path(package + "." + directory, filename) as txt_path:
        return np.loadtxt(txt_path)

# ----------------------
# Rejection sampling 
# ----------------------
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
        pdf_x_try = pdf(x_try)  # Calculate the pdf at the random x values
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

# def create_func_pdf_invcdf(x, y, category="function"):
#     """
#     Function to create a interpolated function, inverse function or inverse cdf from the input x and y.

#     Parameters
#     ----------
#     x : `numpy.ndarray`
#         x values. This has to sorted in ascending order.
#     y : `numpy.ndarray`
#         y values. Corresponding to the x values.
#     category : `str`, optional
#         category of the function. Default is "function". Other options are "function_inverse", "pdf" and "inv_cdf".

#     Returns
#     ----------
#     pdf : `pdf function`
#         interpolated pdf function.
#     inv_pdf : `function inverse`
#         interpolated inverse pdf function.
#     inv_cdf : `function`
#         interpolated inverse cdf.
#     """

#     idx = np.argwhere(np.isnan(y))
#     x = np.delete(x, idx)
#     y = np.delete(y, idx)

#     # create pdf with interpolation
#     pdf_unorm = interp1d(x, y, kind="cubic", fill_value="extrapolate")
#     if category == "function":
#         return pdf_unorm
#     if category == "function_inverse":
#         # create inverse function
#         return interp1d(y, x, kind="cubic", fill_value="extrapolate")

#     min_, max_ = min(x), max(x)
#     norm = quad(pdf_unorm, min_, max_)[0]
#     y = y / norm
#     if category == "pdf" or category is None:
#         # normalize the pdf
#         pdf = interp1d(x, y, kind="cubic", fill_value="extrapolate")
#         return pdf
#     # cdf
#     cdf_values, x, _ = cumulative_trapezoid(y, x, initial=0)

#     inv_cdf = interp1d(cdf_values, x, kind="cubic", fill_value="extrapolate")
#     if category == "inv_cdf":
#         return inv_cdf
#     if category == "all":
#         return [pdf, inv_cdf]

# def create_conditioned_pdf_invcdf(x, conditioned_y, pdf_func, category):
#     """
#     pdf_func is the function to calculate the pdf of x given y
#     x is an array and the output of pdf_func is an array
#     y is the condition
#     we consider parameter plane of x and y

#     Parameters
#     ----------
#     x : `numpy.ndarray`
#         x values.
#     conditioned_y : `numpy.ndarray`
#         conditioned y values.
#     pdf_func : `function`
#         function to calculate the pdf of x given y.
#     category : `str`, optional
#         category of the function. Default is "function". Other options are "function_inverse", "pdf" and "inv_cdf".
#     """

#     list_ = []
#     for y in conditioned_y:
#         phi = pdf_func(x, y)
#         # append pdf for each y along the x-axis
#         list_.append(create_func_pdf_invcdf(x, phi, category=category))

#     return list_

# ------------------------
# interpolator creation
# ------------------------
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
    cdf_values, x, _ = cumulative_spline(y=y, x=x, initial=0)
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
        phi = pdf_func(x, y)
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
    size = x.shape[0]
    for y in conditioned_y:
        phi = pdf_func(x, y * np.ones(size))
        list_.append(create_inv_cdf_array(x, phi))

    return np.array(list_)

def interpolator_from_json(
    identifier_dict,
    directory,
    sub_directory,
    name,
    x,
    pdf_func=None,
    y=None,
    conditioned_y=None,
    dimension=1,
    category="pdf",
    create_new=False,
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
        if dimension == 1:
            if y is None:
                y = pdf_func(x)
            if category == "function":
                interpolator = create_func(x, y)
            elif category == "function_inverse":
                interpolator = create_func_inv(x, y)
            elif category == "pdf":
                interpolator = create_pdf(x, y)
            elif category == "inv_cdf":
                interpolator = create_inv_cdf_array(x, y)
            else:
                raise ValueError(
                    "The category given should be function, function_inverse, pdf or inv_cdf."
                )
        elif dimension == 2:
            if category == "pdf":
                interpolator = create_conditioned_pdf(x, conditioned_y, pdf_func)
            elif category == "inv_cdf":
                interpolator = create_conditioned_inv_cdf_array(
                    x, conditioned_y, pdf_func
                )
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

# ------------------------
# Integration
# ------------------------
# NOTE: cache=True intentionally NOT set here. ``function`` and
# ``uniform_prior`` are first-class Numba dispatchers, whose types vary
# per call site, so Numba cannot persist a stable cache key for this
# function across processes.
@njit(fastmath=True)
def monte_carlo_integration(function, uniform_prior, size=10000):
    """
    Function to perform Monte Carlo integration.

    Parameters
    ----------
    function : `function`
        function to be integrated.
    uniform_prior : `function`
        uniform prior function.

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
    integral = np.mean(y) * (np.max(x) - np.min(x))
    return integral

# ------------------------
# Interpolation
# ------------------------
@njit(cache=True, fastmath=True)
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
    i[idx2] = x.shape[0] - 2

    # Calculate the relative position within the interval
    dx = xnew - x[i]

    # Calculate the interpolated value
    # Cubic polynomial: a + b*dx + c*dx^2 + d*dx^3
    a, b, c, d = coefficients[:, i]
    # result = a + b*dx + c*dx**2 + d*dx**3
    result = d + c * dx + b * dx**2 + a * dx**3

    return result

@njit(cache=True, fastmath=True)
def cubic_spline_interpolator_scalar(xnew, coefficients, x):
    """
    Function to interpolate using cubic spline.

    Parameters
    ----------
    xnew : `float`
        new x value.
    coefficients : `numpy.ndarray`
        coefficients of the cubic spline.
    x : `numpy.ndarray`
        x values.

    Returns
    ----------
    result : `float`
        interpolated value.
    """

    # Handling extrapolation
    i = np.searchsorted(x, xnew) - 1

    if xnew <= x[0]:
        i = 0
    elif xnew > x[-1]:
        i = x.shape[0] - 2

    # Calculate the relative position within the interval
    dx = xnew - x[i]

    # Calculate the interpolated value
    a, b, c, d = coefficients[:, i]
    # result = a + b*dx + c*dx**2 + d*dx**3
    result = d + c * dx + b * dx**2 + a * dx**3
    return result

@njit(cache=True, fastmath=True)
def pdf_cubic_spline_interpolator(xnew, norm_const, coefficients, x):
    """
    Function to interpolate pdf using cubic spline.

    Parameters
    ----------
    xnew : `numpy.ndarray`
        new x values.
    norm_const : `float`
        normalization constant.
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
    i[idx2] = x.shape[0] - 2

    # Calculate the relative position within the interval
    dx = xnew - x[i]

    # Calculate the interpolated value
    # Cubic polynomial: a + b*dx + c*dx^2 + d*dx^3
    a, b, c, d = coefficients[:, i]
    result = d + c * dx + b * dx**2 + a * dx**3
    
    return result / norm_const

@njit(cache=True, fastmath=True)
def cubic_hermite_interpolation(x, x_array, y_array):
    """
    Function to perform cubic hermite interpolation to find y.
    This is a scalar function.

    Parameters
    ----------
    x : `float`
        x value to find the corresponding y value.
    x_array : `numpy.ndarray`
        x-axis values of the known 4 points with x_array[1] <= x <= x_array[2].
    y_array : `numpy.ndarray`
        y-axis values of the known 4 points.

    Returns
    -------
    y : `float`
        Interpolated y value.
    """
    x0, x1, x2, x3 = x_array
    y0, y1, y2, y3 = y_array

    # Compute tangents (derivatives) at x1 and x2 using finite differences
    m1 = ((y2 - y1) / (x2 - x1)) * ((x1 - x0) / (x2 - x0)) + ((y1 - y0) / (x1 - x0)) * ((x2 - x1) / (x2 - x0))
    m2 = ((y3 - y2) / (x3 - x2)) * ((x2 - x1) / (x3 - x1)) + ((y2 - y1) / (x2 - x1)) * ((x3 - x2) / (x3 - x1))  

    # Compute the relative position within the interval
    denom = x2 - x1
    t = (x - x1) / denom

    # Hermite basis polynomials
    h00 = 2.0 * t**3 - 3.0 * t**2 + 1.0
    h10 = t**3 - 2.0 * t**2 + t
    h01 = -2.0 * t**3 + 3.0 * t**2
    h11 = t**3 - t**2

    # Final interpolated value
    y = h00 * y1 + h10 * m1 * denom + h01 * y2 + h11 * m2 * denom
    return y

@njit(cache=True, fastmath=True)
def cubic_spline_interpolator2d(xnew_array, ynew_array, coefficients, x, y):
    """
    Function to calculate the interpolated value given the conditioned variable (ynew) and primary variable (xnew). This is based off 2D bicubic spline interpolation.

    Parameters
    ----------
    xnew_array : `numpy.ndarray`
        New x values at which to interpolate.
    ynew_array : `numpy.ndarray`
        New y values (conditioned variable) at which to interpolate.
    coefficients : `numpy.ndarray`
        Array of coefficients for the cubic spline interpolation.
    x : `numpy.ndarray`
        Array of x values for the coefficients.
    y : `numpy.ndarray`
        Array of y values (conditioned variable) for the coefficients.

    Returns
    -------
    result : `numpy.ndarray`
        Interpolated values.
    """
    n_samples = xnew_array.shape[0]
    result_array = np.empty(n_samples, dtype=np.float64)
    y_idx_arr = np.searchsorted(y, ynew_array) - 1 # left side index
    len_y = y.shape[0]
    y_idx_arr = np.clip(y_idx_arr, 0, len_y - 2) # clip to ensure valid index range
    
    for i in range(n_samples):
        xnew = xnew_array[i]
        ynew = ynew_array[i]
        xnew_single = np.empty(1, dtype=np.float64)  # local per-thread to avoid race condition
        xnew_single[0] = xnew

        y_idx = y_idx_arr[i]
        condi = (y_idx == 0) or (y_idx == len_y - 2)

        if not condi:
            y_idx1 = y_idx - 1
            y_idx2 = y_idx
            y_idx3 = y_idx + 1
            y_idx4 = y_idx + 2
            # print(f"d) idx1, idx2, idx3, idx4 = {y_idx1}, {y_idx2}, {y_idx3}, {y_idx4}")
            y1, y2, y3, y4 = y[y_idx1], y[y_idx2], y[y_idx3], y[y_idx4]
            z1 = cubic_spline_interpolator(
                xnew_single, coefficients[y_idx1], x[y_idx1]
            )[0]
            z2 = cubic_spline_interpolator(
                xnew_single, coefficients[y_idx2], x[y_idx2]
            )[0]
            z3 = cubic_spline_interpolator(
                xnew_single, coefficients[y_idx3], x[y_idx3]
            )[0]
            z4 = cubic_spline_interpolator(
                xnew_single, coefficients[y_idx4], x[y_idx4]
            )[0]

            # # ---------------------------------------------------------
            # # Cubic Spline Interpolation with Coefficients (LER method)
            # # ---------------------------------------------------------
            # coeff_low, coeff_high = 4, 8
            # coeff = coefficients_generator_ler(y1, y2, y3, y4, z1, z2, z3, z4)
            # matrixD = coeff[coeff_low:coeff_high]
            # matrixB = np.array([ynew**3, ynew**2, ynew, 1.0])
            # result_array[i] = np.dot(matrixB, matrixD)

            # ---------------------------------------------------------
            # Cubic Hermite Interpolation (no coefficients needed, but requires tangents)
            # ---------------------------------------------------------
            result_array[i] = cubic_hermite_interpolation(ynew, np.array([y1, y2, y3, y4]), np.array([z1, z2, z3, z4]))
        # lower or upper point
        else:
            if y_idx == 0:  # lower end point
                y_idx1 = y_idx 
                y_idx2 = y_idx + 1
            else:  # upper end point
                y_idx1 = y_idx - 1
                y_idx2 = y_idx

            y1, y2 = y[y_idx1], y[y_idx2]
            z1 = cubic_spline_interpolator(
                xnew_single, coefficients[y_idx1], x[y_idx1]
            )[0]
            z2 = cubic_spline_interpolator(
                xnew_single, coefficients[y_idx2], x[y_idx2]
            )[0]
            # use linear interpolation for the lower end point
            result_array[i] = z1 + (z2 - z1) * (ynew - y1) / (y2 - y1)

    return result_array

@njit(cache=True, fastmath=True)
def pdf_cubic_spline_interpolator2d(
    xnew_array, ynew_array, norm_array, coefficients, x, y
):
    """
    Function to calculate the interpolated PDF value given the conditioned variable (ynew) and primary variable (xnew). This is based off 2D bicubic spline interpolation with per-slice normalization.

    Parameters
    ----------
    xnew_array : `numpy.ndarray`
        New x values at which to interpolate.
    ynew_array : `numpy.ndarray`
        New y values (conditioned variable) at which to interpolate.
    norm_array : `numpy.ndarray`
        Array of normalization constants for each y slice.
    coefficients : `numpy.ndarray`
        Array of coefficients for the cubic spline interpolation.
    x : `numpy.ndarray`
        Array of x values for the coefficients.
    y : `numpy.ndarray`
        Array of y values (conditioned variable) for the coefficients.

    Returns
    -------
    result : `numpy.ndarray`
        Interpolated PDF values (clipped to non-negative).
    """
    n_samples = xnew_array.shape[0]
    result_array = np.empty(n_samples, dtype=np.float64)
    y_idx_arr = np.searchsorted(y, ynew_array) - 1 # left side index
    len_y = y.shape[0]
    y_idx_arr = np.clip(y_idx_arr, 0, len_y - 2) # clip to ensure valid index range
    
    for i in range(n_samples):
        xnew = xnew_array[i]
        ynew = ynew_array[i]
        xnew_single = np.empty(1, dtype=np.float64)  # local per-thread to avoid race condition
        xnew_single[0] = xnew

        # find the index nearest to the ynew in y
        # y_idx = np.searchsorted(y, ynew) - 1 if ynew > y[0] else 0
        y_idx = y_idx_arr[i]

        condi = (y_idx == 0) or (y_idx == len_y - 2)

        if not condi:
            y_idx1 = y_idx - 1
            y_idx2 = y_idx
            y_idx3 = y_idx + 1
            y_idx4 = y_idx + 2
            # print(f"d) idx1, idx2, idx3, idx4 = {y_idx1}, {y_idx2}, {y_idx3}, {y_idx4}")
            y1, y2, y3, y4 = y[y_idx1], y[y_idx2], y[y_idx3], y[y_idx4]
            z1 = pdf_cubic_spline_interpolator(
                xnew_single, norm_array[y_idx1], coefficients[y_idx1], x[y_idx1]
            )[0]
            z2 = pdf_cubic_spline_interpolator(
                xnew_single, norm_array[y_idx2], coefficients[y_idx2], x[y_idx2]
            )[0]
            z3 = pdf_cubic_spline_interpolator(
                xnew_single, norm_array[y_idx3], coefficients[y_idx3], x[y_idx3]
            )[0]
            z4 = pdf_cubic_spline_interpolator(
                xnew_single, norm_array[y_idx4], coefficients[y_idx4], x[y_idx4]
            )[0]

            # # ---------------------------------------------------------
            # # Cubic Spline Interpolation with Coefficients (LER method)
            # # ---------------------------------------------------------
            # coeff_low, coeff_high = 4, 8
            # coeff = coefficients_generator_ler(y1, y2, y3, y4, z1, z2, z3, z4)
            # matrixD = coeff[coeff_low:coeff_high]
            # matrixB = np.array([ynew**3, ynew**2, ynew, 1.0])
            # result_array[i] = np.dot(matrixB, matrixD)

            # ---------------------------------------------------------
            # Cubic Hermite Interpolation (no coefficients needed, but requires tangents)
            # ---------------------------------------------------------
            result_array[i] = cubic_hermite_interpolation(ynew, np.array([y1, y2, y3, y4]), np.array([z1, z2, z3, z4]))
        # lower or upper point
        else:
            if y_idx == 0:  # lower end point
                y_idx1 = y_idx 
                y_idx2 = y_idx + 1
            else:  # upper end point
                y_idx1 = y_idx - 1
                y_idx2 = y_idx

            y1, y2 = y[y_idx1], y[y_idx2]
            z1 = pdf_cubic_spline_interpolator(
                xnew_single, norm_array[y_idx1], coefficients[y_idx1], x[y_idx1]
            )[0]
            z2 = pdf_cubic_spline_interpolator(
                xnew_single, norm_array[y_idx2], coefficients[y_idx2], x[y_idx2]
            )[0]
            # use linear interpolation for the lower end point
            result_array[i] = z1 + (z2 - z1) * (ynew - y1) / (y2 - y1)

    # Clip negative values to zero
    result_array = np.clip(result_array, a_min=0, a_max=None)

    return result_array

# @njit
# def coefficients_generator_ler(y1, y2, y3, y4, z1, z2, z3, z4):
#     """
#     Function to generate the coefficients for the cubic spline interpolation of fn(y)=z.

#     Parameters
#     ----------
#     y1, y2, y3, y4, z1, z2, z3, z4: `float`
#         Values of y and z for the cubic spline interpolation.

#     Returns
#     -------
#     coefficients: `numpy.ndarray`
#         Coefficients for the cubic spline interpolation.
#     """
#     matrixA = np.array(
#         [
#             [y1**3, y1**2, y1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
#             [y2**3, y2**2, y2, 1, 0, 0, 0, 0, 0, 0, 0, 0],
#             [0, 0, 0, 0, y2**3, y2**2, y2, 1, 0, 0, 0, 0],
#             [0, 0, 0, 0, y3**3, y3**2, y3, 1, 0, 0, 0, 0],
#             [0, 0, 0, 0, 0, 0, 0, 0, y3**3, y3**2, y3, 1],
#             [0, 0, 0, 0, 0, 0, 0, 0, y4**3, y4**2, y4, 1],
#             [3 * y2**2, 2 * y2, 1, 0, -3 * y2**2, -2 * y2, -1, 0, 0, 0, 0, 0],
#             [0, 0, 0, 0, 3 * y3**2, 2 * y3, 1, 0, -3 * y3**2, -2 * y3, -1, 0],
#             [6 * y2, 2, 0, 0, -6 * y2, -2, 0, 0, 0, 0, 0, 0],
#             [0, 0, 0, 0, 6 * y3, 2, 0, 0, -6 * y3, -2, 0, 0],
#             [6 * y1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#             [0, 0, 0, 0, 0, 0, 0, 0, 6 * y4, 2, 0, 0],
#         ]
#     )
#     matrixC = np.array([z1, z2, z2, z3, z3, z4, 0, 0, 0, 0, 0, 0])
#     return np.dot(np.linalg.inv(matrixA), matrixC)

@njit(parallel=True)
def inverse_transform_sampler2d(size, conditioned_y, cdf2d, x2d, y1d, trim=True):
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
    trim: `bool`
        Whether to trim the CDF and x values. This trims leading zeros and trailing ones from the cdf.

    Returns
    -------
    samples: `numpy.ndarray`
        Samples of the conditioned y.
    """

    samples = np.zeros(size)

    n_y = y1d.shape[0]
    n_x = x2d.shape[1]
    y_idx_arr = np.searchsorted(y1d, conditioned_y)
    y_idx_arr = np.clip(y_idx_arr, 1, n_y - 1)

    # for i, y in enumerate(conditioned_y):
    for i in prange(size):
        
        y = conditioned_y[i]
        # find the index nearest to the conditioned_y in y1d
        y_idx = y_idx_arr[i]
        

        # linear scaling
        # new cdf values
        y1, y0, z1, z0 = y1d[y_idx], y1d[y_idx - 1], cdf2d[y_idx], cdf2d[y_idx - 1]
        cdf1d = z0 + (z1 - z0) * (y - y0) / (y1 - y0)
        # new x values
        # x=x2d[0], if all rows are same
        x1, x0 = x2d[y_idx], x2d[y_idx - 1]
        x1d = x0 + (x1 - x0) * (y - y0) / (y1 - y0)

        # sampling part
        if trim:
            cdf1d, x1d = trim_cdf(cdf1d, x1d)

        while True:
            u = np.random.uniform(0, 1)
            idx = np.searchsorted(cdf1d, u)
            len_cdf = cdf1d.shape[0]
            # Clamp idx to valid range (np.clip on scalar int is not supported by Numba)
            if idx < 1:
                idx = 1
            elif idx > len_cdf - 1:
                idx = len_cdf - 1

            cdf_idx_hi, cdf_idx_lo = cdf1d[idx], cdf1d[idx - 1]
            x_idx_hi, x_idx_lo = x1d[idx], x1d[idx - 1]

            samples[i] = x_idx_lo + (x_idx_hi - x_idx_lo) * (u - cdf_idx_lo) / (cdf_idx_hi - cdf_idx_lo)

            if not np.isnan(samples[i]):
                break

    return samples

@njit(parallel=True)
def inverse_transform_sampler2d_spline(size, conditioned_y, cdf2d, x2d, y1d, trim=True):
    """
    Sample from a 2D inverse transform with spline interpolation in both dimensions.

    Performs 1D inverse transform sampling for a set of conditioned y values by \n
    interpolating the 2D CDF surface and using cubic spline-based inverse sampling.

    Parameters
    ----------
    size : `int`
        Number of samples to generate.
    conditioned_y : `numpy.ndarray`
        Conditioned y values (1D array of length size).
    cdf2d : `numpy.ndarray`
        2D array of CDF values (shape: len(y1d), len(x2d[0])).
    x2d : `numpy.ndarray`
        2D array of x-axis points (shape: len(y1d), len(x2d[0])).
    y1d : `numpy.ndarray`
        1D array of y-axis conditioning points.
    trim : `bool`, optional
        Whether to trim CDF leading zeros and trailing ones. Default is True.

    Returns
    -------
    samples : `numpy.ndarray`
        Generated samples (1D array of length size).
    """

    samples = np.zeros(size)
    n_y = y1d.shape[0]
    n_x = x2d.shape[1]
    
    y_idx_arr = np.searchsorted(y1d, conditioned_y) - 1  # left side index
    y_idx_arr = np.clip(y_idx_arr, 0, n_y - 2)  # clip to ensure valid index range
    
    for i in prange(size):
        y_i = conditioned_y[i]
        y_idx = y_idx_arr[i]
        

        condi = (y_idx == 0) or (y_idx == n_y - 2)

        if not condi:
            # Hermite interpolation for y direction (using 4 points)
            y_idx_list = np.array([y_idx - 1, y_idx, y_idx + 1, y_idx + 2], dtype=np.int64)
            y_arr = y1d[y_idx_list]

            cdf_new = np.zeros(n_x)
            for j in range(n_x):
                cdf_arr = cdf2d[y_idx_list, j]
                cdf_new[j] = cubic_hermite_interpolation(y_i, y_arr, cdf_arr)

            # new x values
            x1d = x2d[y_idx] + (x2d[y_idx + 1] - x2d[y_idx]) * (y_i - y_arr[1]) / (y_arr[2] - y_arr[1])
        else:
            # linear interpolation
            if y_idx == 0:  # lower end point
                y1 = y1d[y_idx]
                y2 = y1d[y_idx + 1]
                z1 = cdf2d[y_idx]
                z2 = cdf2d[y_idx + 1]
            else:  # upper end point
                y1 = y1d[y_idx - 1]
                y2 = y1d[y_idx]
                z1 = cdf2d[y_idx - 1]
                z2 = cdf2d[y_idx]

            # z1,z2 are 1d array, while y1,y2 are floats
            cdf_new = z1 + (z2 - z1) * (y_i - y1) / (y2 - y1)

            # new x values; keep interpolation consistent with the selected y bracket
            if y_idx == 0:
                x1d = x2d[y_idx] + (x2d[y_idx + 1] - x2d[y_idx]) * (y_i - y1) / (y2 - y1)
            else:
                x1d = x2d[y_idx - 1] + (x2d[y_idx] - x2d[y_idx - 1]) * (y_i - y1) / (y2 - y1)

        # Sample using spline-based inverse transform
        if trim:
            cdf_new, x = trim_cdf(cdf_new, x1d)
        else:
            x = x1d
        
        # inverse cdf spline (cdf → x)
        coeff = get_pchip_spline_coeffs(x=cdf_new, y=x)
        
        # Sample one point using the spline
        u = np.random.uniform(0, 1)
        x_sample = cubic_spline_interpolator_scalar(u, coeff, cdf_new)
        samples[i] = x_sample
    
    return samples

@njit(cache=True, fastmath=True)
def inverse_transform_sampler(size, cdf, x, trim=True):
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

    if trim:
        cdf, x = trim_cdf(cdf, x)

    u = np.random.uniform(0, 1, size)
    idx = np.searchsorted(cdf, u)

    # Clip idx to valid range to prevent boundary issues:
    # - When idx=0 (u <= cdf[0]), idx-1 would wrap to -1 (last element)
    # - When idx=len(cdf) (u > cdf[-1]), idx would be out of bounds
    idx = np.clip(idx, 1, len(cdf) - 1)
    x1, x0, y1, y0 = cdf[idx], cdf[idx - 1], x[idx], x[idx - 1]
    samples = y0 + (y1 - y0) * (u - x0) / (x1 - x0)

    return samples

@njit(cache=True, fastmath=True)
def inverse_transform_sampler_spline(size, cdf, x, trim=True):
    """
    Function to sample from the inverse transform method using spline interpolation.

    Parameters
    ----------
    size : `int`
        number of samples.
    cdf : `numpy.ndarray`
        cdf values.
    x : `numpy.ndarray`
        x values corresponding to the cdf.
    trim : `bool`
        Whether to trim leading zeros and trailing ones from the cdf.

    Returns
    ----------
    samples : `numpy.ndarray`
        samples from the cdf.
    """

    if trim:
        cdf, x = trim_cdf(cdf, x)

    cdf_new = np.random.uniform(0, 1, size)

    inv_coeff = get_pchip_spline_coeffs(x=cdf, y=x)
    
    x_new = cubic_spline_interpolator(cdf_new, inv_coeff, cdf)

    return x_new

@njit(cache=True, fastmath=True)
def trim_cdf(cdf, x):
    """
    Function to trim leading zeros and trailing ones from the cdf, but guarantee at least 2 points for downstream spline interpolation.

    Parameters
    ----------
    cdf : `numpy.ndarray`
        CDF values.
    x : `numpy.ndarray`
        x values corresponding to the CDF values.

    Returns
    -------
    trimmed_cdf : `numpy.ndarray`
        Trimmed CDF values.
    trimmed_x : `numpy.ndarray`
        x values corresponding to the trimmed CDF values.
    """

    zeros = np.where(cdf <= 0)[0]
    ones = np.where(cdf >= 1)[0]
    idx_lower = zeros[-1] if zeros.shape[0] > 0 else 0
    idx_upper = (ones[0] + 1) if ones.shape[0] > 0 else len(cdf)
    
    # Safety: ensure at least 2 output points for downstream spline interpolation
    if idx_upper - idx_lower < 2:
        idx_lower = 0
        idx_upper = cdf.shape[0]

    return cdf[idx_lower:idx_upper], x[idx_lower:idx_upper]

# ----------------------
# Normal PDF
# ----------------------
# 1D normal PDF
@njit(cache=True, fastmath=True)
def normal_pdf(x, mean=0.0, std=0.05):
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

    return (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)

# 2D normal PDF
@njit(cache=True, fastmath=True)
def normal_pdf_2d(x, y, mean_x=0.0, std_x=0.05, mean_y=0.0, std_y=0.05):
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

    factor_x = (1 / (std_x * np.sqrt(2 * np.pi))) * np.exp(
        -0.5 * ((x - mean_x) / std_x) ** 2
    )
    factor_y = (1 / (std_y * np.sqrt(2 * np.pi))) * np.exp(
        -0.5 * ((y - mean_y) / std_y) ** 2
    )
    return factor_x * factor_y

# def load_txt_from_module(package, directory, filename):
#     """ """

#     with resources.path(package + "." + directory, filename) as txt_path:
#         return np.loadtxt(txt_path)

# ----------------------
# CDF from PDF
# ----------------------
@njit(parallel=True, fastmath=True)
def cumulative_spline(x, y, initial=0.0):
    """
    Compute CDF from a PDF array using cubic spline + analytical integration.
    """
    n = x.shape[0]
    m = n - 1
    cs = get_pchip_spline_coeffs(x, y)  # Get coefficients for each interval

    # 1. Array to store the independent area of each interval
    integrals = np.zeros(n, dtype=np.float64)
    
    # Set the starting value
    integrals[0] = initial
    
    # 2. Parallel loop: Compute each interval's integral independently
    for i in prange(m):
        h = x[i + 1] - x[i]
        # Integrate S_i(u) = d[i] + c[i]*u + b[i]*u^2 + a[i]*u^3 from 0 to h
        integrals[i+1] = cs[3, i] * h + cs[2, i] * h**2 / 2.0 + cs[1, i] * h**3 / 3.0 + cs[0, i] * h**4 / 4.0

    # 3. Sequential cumulative sum
    cdf = np.cumsum(integrals)
    norm = cdf[m]
    cdf /= norm  # Normalize to ensure the last value is 1.0

    return cdf, x, norm

@njit(cache=True, fastmath=True)
def cumulative_trapezoid(y, x, initial=0.0):
    """
    Compute the cumulative integral using the trapezoidal rule. The output is conditioned to be between 0 and 1, and the leading zeros and trailing ones are trimmed. This is used for computing the inverse CDF coefficients. 

    Parameters
    ----------
    y : ``numpy.ndarray``
        Function values to integrate.
    x : ``numpy.ndarray``
        x-coordinates corresponding to y values.
    initial : ``float``
        Initial value for the cumulative sum. \n
        default: 0.0

    Returns
    -------
    cumsum : ``numpy.ndarray``
        Cumulative integral values.
    """

    n = x.shape[0]
    m = int(n-1)
    cdf = np.zeros(n, dtype=np.float64)
    cdf[0] = initial
    # for i in range(1, len(y)):
    #   cdf[i] = cdf[i - 1] + (y[i - 1] + y[i]) * (x[i] - x[i - 1]) / 2.0
    dx = np.diff(x)
    cdf[1:] = np.cumsum((y[:m] + y[1:]) * dx / 2.0) + initial

    norm = cdf[m]
    cdf /= norm  # Normalize to ensure the last value is 1.0

    return cdf, x, norm

@njit(parallel=True, fastmath=True)
def get_pchip_spline_coeffs(x, y):
    """
    Computes the PCHIP cubic spline coefficients for 1D arrays x and y.
    Returns a (4, N-1) array matching SciPy's PchipInterpolator(x, y).c
    """
    n = len(x)
    c = np.empty((4, n - 1), dtype=np.float64)
    
    if n < 2:
        raise ValueError("At least 2 points are required for interpolation.")
        
    h = np.empty(n - 1, dtype=np.float64)
    s = np.empty(n - 1, dtype=np.float64)
    d = np.empty(n, dtype=np.float64)

    # 1. Compute secant slopes
    for i in range(n - 1):
        h[i] = x[i + 1] - x[i]
        s[i] = (y[i + 1] - y[i]) / h[i]

    # 2. Compute PCHIP derivatives at each point
    # Interior points (weighted harmonic mean)
    for i in range(1, n - 1):
        # If secant slopes change sign or either is zero, derivative is zero
        if s[i - 1] * s[i] <= 0.0:
            d[i] = 0.0
        else:
            w1 = 2.0 * h[i] + h[i - 1]
            w2 = h[i] + 2.0 * h[i - 1]
            d[i] = (w1 + w2) / (w1 / s[i - 1] + w2 / s[i])

    # Endpoints (Shape-preserving 3-point extrapolation)
    if n == 2:
        d[0] = s[0]
        d[1] = s[0]
    else:
        # Left endpoint
        d0 = ((2.0 * h[0] + h[1]) * s[0] - h[0] * s[1]) / (h[0] + h[1])
        if np.sign(d0) != np.sign(s[0]):
            d[0] = 0.0
        elif np.sign(s[0]) != np.sign(s[1]) and abs(d0) > 3.0 * abs(s[0]):
            d[0] = 3.0 * s[0]
        else:
            d[0] = d0

        # Right endpoint
        dN = ((2.0 * h[-1] + h[-2]) * s[-1] - h[-1] * s[-2]) / (h[-1] + h[-2])
        if np.sign(dN) != np.sign(s[-1]):
            d[-1] = 0.0
        elif np.sign(s[-1]) != np.sign(s[-2]) and abs(dN) > 3.0 * abs(s[-1]):
            d[-1] = 3.0 * s[-1]
        else:
            d[-1] = dN

    # 3. Compute the coefficients for each interval
    # P(x) = c0*(x-xi)^3 + c1*(x-xi)^2 + c2*(x-xi) + c3
    for i in prange(n - 1):
        c[0, i] = (d[i] + d[i + 1] - 2.0 * s[i]) / (h[i] * h[i])
        c[1, i] = (3.0 * s[i] - 2.0 * d[i] - d[i + 1]) / h[i]
        c[2, i] = d[i]
        c[3, i] = y[i]
        
    return c

# ----------------------
# batch sampling
# ----------------------
def batch_handler(
    size,
    batch_size,
    sampling_routine,
    output_jsonfile,
    save_batch=True,
    resume=False,
    param_name="parameters",
):
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
    if not output_jsonfile:
        # no output file requested — skip all file I/O
        save_batch = False
        resume = False
        dict_buffer = None
    elif resume and os.path.exists(output_jsonfile):
        # get sample from json file
        dict_buffer = get_param_from_json(output_jsonfile)
    else:
        remove_file(output_jsonfile)
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
    freshly_sampled = False  # track whether initial batch was just sampled

    if not resume:
        # create new first batch with the frac_batches
        track_batches, dict_buffer = create_batch_params(
            sampling_routine,
            frac_batches,
            dict_buffer,
            save_batch,
            output_jsonfile,
            track_batches=track_batches,
        )
        freshly_sampled = True
    else:
        # check where to resume from
        # identify the last batch and assign current batch number
        # try-except is added to avoid the error when the file does not exist or if the file is empty or corrupted or does not have the required key.
        try:
            print(f"resuming from {output_jsonfile}")
            len_ = len(list(dict_buffer.values())[0])
            track_batches = (len_ - frac_batches) // batch_size + 1
        except (IndexError, TypeError, KeyError, AttributeError):
            # create new first batch with the frac_batches
            track_batches, dict_buffer = create_batch_params(
                sampling_routine,
                frac_batches,
                dict_buffer,
                save_batch,
                output_jsonfile,
                track_batches=track_batches,
            )
            freshly_sampled = True

    # loop over the remaining batches
    min_, max_ = track_batches, num_batches
    # print(f"min_ = {min_}, max_ = {max_}")
    save_param = False
    if min_ == max_:
        if freshly_sampled and not save_batch:
            # Data was just sampled in the initial batch but not yet written to file
            save_param = True
        else:
            print(f"{param_name} already sampled.")
            save_param = False
    elif min_ > max_:
        len_ = len(list(dict_buffer.values())[0])
        print(
            f"existing {param_name} size is {len_} is more than the required size={size}. It will be trimmed."
        )
        dict_buffer = trim_dictionary(dict_buffer, size)
        save_param = True
    else:
        for i in range(min_, max_):
            _, dict_buffer = create_batch_params(
                sampling_routine,
                batch_size,
                dict_buffer,
                save_batch,
                output_jsonfile,
                track_batches=i,
                resume=True,
            )

        if save_batch:
            # if save_batch=True, then dict_buffer is only the last batch
            # batch saving is already done in create_batch_params function
            dict_buffer = get_param_from_json(output_jsonfile)
        else:  # dont save in batches
            # this if condition is required if there is nothing to save
            save_param = True

    if save_param and output_jsonfile:
        # store all params in json file
        print(f"saving all {param_name} in {output_jsonfile} ")
        append_json(output_jsonfile, dict_buffer, replace=True)

    return dict_buffer

def create_batch_params(
    sampling_routine,
    frac_batches,
    dict_buffer,
    save_batch,
    output_jsonfile,
    track_batches,
    resume=False,
):
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
    param = sampling_routine(
        size=frac_batches,
        save_batch=save_batch,
        output_jsonfile=output_jsonfile,
        resume=resume,
    )

    # adding batches and hold it in the buffer attribute.
    if not save_batch:
        # in the new batch (new sampling run), dict_buffer will be None
        if dict_buffer is None:
            dict_buffer = param
        else:
            for key, value in param.items():
                try:
                    dict_buffer[key] = np.concatenate((dict_buffer[key], value))
                except Exception as e:
                    raise ValueError(
                        f"For key {key}, concatenate failed for dictionary 1 {dict_buffer[key]} and dictionary 2 {value}. Error: {e}"
                    ) from e
    else:
        # store all params in json file
        print(f"saving batch {track_batches} in {output_jsonfile} ")
        dict_buffer = append_json(
            file_name=output_jsonfile,
            new_dictionary=param,
            old_dictionary=dict_buffer,
            replace=not (resume),
        )

    return track_batches, dict_buffer


def KStest(
    lens_param1: dict,
    lens_param2: dict,
    *,
    keys=None,
    alternative: str = "two-sided",
    mode: str = "auto",
    nan_policy: str = "omit",
    return_pvalue: bool = True,
):
    """Two-sample Kolmogorov-Smirnov tests for lens-parameter dictionaries.

    For each key in common (or the keys you pass), the samples are compared
    with ``scipy.stats.ks_2samp``. The KS test is non-parametric: it does not
    assume a parametric form for the parent distribution(s).

    The reported statistic ``D`` (SciPy's ``statistic``) is the maximum absolute
    gap between the two empirical cumulative distribution functions. It ranges
    from 0 to 1. Smaller ``D`` means the two empirical CDFs track each other
    more closely; under the null hypothesis that both samples are i.i.d. draws
    from the same *continuous* distribution, large ``D`` implies a small
    p-value (evidence against that null). Interpretation depends on sample size
    and on whether the two-sided or one-sided ``alternative`` you choose is
    appropriate for your science question.

    Parameters
    ----------
    lens_param1, lens_param2 : dict
        Dicts like the output of `ler.sample_lens_parameters`, mapping parameter
        names -> array-like samples.
    keys : iterable[str] | None
        Which keys to test. If None, uses the intersection of keys.
    alternative : str
        Passed to `scipy.stats.ks_2samp`.
    mode : str
        Passed to `scipy.stats.ks_2samp`.
    nan_policy : {'omit','propagate'}
        If 'omit', drops non-finite values before testing.
    return_pvalue : bool
        If True, return both KS statistic and p-value.

    Returns
    -------
    out : dict
        out[key] = {'D': <ks statistic>, 'pvalue': <pvalue>, 'n1': <int>, 'n2': <int>}
        If return_pvalue=False, 'pvalue' is omitted.
    """
    try:
        from scipy.stats import ks_2samp
    except Exception as e:
        raise ImportError(
            "KStest requires scipy. Install it (e.g. `pip install scipy`)."
        ) from e

    if keys is None:
        keys = sorted(set(lens_param1.keys()).intersection(lens_param2.keys()))

    out = {}
    for k in keys:
        x = np.asarray(lens_param1[k])
        y = np.asarray(lens_param2[k])

        if nan_policy == "omit":
            x = x[np.isfinite(x)]
            y = y[np.isfinite(y)]

        n1 = int(x.size)
        n2 = int(y.size)
        if n1 == 0 or n2 == 0:
            out[k] = {"D": np.nan, "pvalue": np.nan, "n1": n1, "n2": n2} if return_pvalue else {"D": np.nan, "n1": n1, "n2": n2}
            continue

        res = ks_2samp(x, y, alternative=alternative, mode=mode)
        if return_pvalue:
            out[k] = {"D": float(res.statistic), "pvalue": float(res.pvalue), "n1": n1, "n2": n2}
        else:
            out[k] = {"D": float(res.statistic), "n1": n1, "n2": n2}

    return out