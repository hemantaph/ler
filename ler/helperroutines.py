# -*- coding: utf-8 -*-
""" 
This module contains helper routines for other modules in the ler package.
"""

import numpy as np
import json

chunk_size = 10000


class NumpyEncoder(json.JSONEncoder):
    """
    class for storing a numpy.ndarray or any nested-list composition as JSON file

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
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def append_json(file_name, dictionary, replace=False):
    """Append and update a json file with a dictionary.

    Parameters
    ----------
    file_name : `str`
        json file name for storing the parameters.
    dictionary : `dict`
        dictionary to be appended to the json file.
    replace : `bool`, optional
        If True, replace the json file with the dictionary. Default is False.

    """

    # check if the file exists
    try:
        with open(file_name, "r", encoding="utf-8") as f:
            data = json.load(f)
    except:
        "File does not exist. Creating a new one..."
        replace = True

    if replace:
        data = dictionary
    else:
        data_key = data.keys()
        for key, value in dictionary.items():
            if key in data_key:
                data[key] = np.concatenate((data[key], value))

    json_dump = json.dumps(data, cls=NumpyEncoder)
    with open(file_name, "w", encoding="utf-8") as write_file:
        json.dump(json.loads(json_dump), write_file, indent=4)


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


def rejection_sample(pdf, xmin, xmax, size=100):
    """
    Helper function for rejection sampling from a pdf with maximum and minimum arguments.
    Input parameters:
        pdf: the pdf to sample from
        xmin: the minimum argument of the pdf
        xmax: the maximum argument of the pdf
        size: the number of samples to draw
    Output:
        samples: the samples drawn from the pdf
    """
    x = np.linspace(xmin, xmax, 1000)
    y = pdf(x)
    ymax = np.max(y)
    # Rejection sample in chunks
    x_sample = []
    while len(x_sample) < size:
        x_try = np.random.uniform(xmin, xmax, size=chunk_size)
        y_try = np.random.uniform(0, ymax, size=chunk_size)
        ymax = max(ymax, np.max(y_try))
        # Add while retaining 1D shape of the list
        x_sample += list(x_try[y_try < pdf(x_try)])
    # Transform the samples to a 1D numpy array
    x_sample = np.array(x_sample).flatten()
    # Return the correct number of samples
    return x_sample[:size]

def rejection_sample2d(pdf, xmin, xmax, ymin, ymax, size=100):
    
    chunk_size = 10000
    
    x = np.linspace(xmin, xmax, 1000)
    y = np.linspace(ymin, ymax, 1000)
    z = pdf(x,y)
    zmax = np.max(z)
    
    
    # Rejection sample in chunks
    x_sample = []
    y_sample = []
    while len(x_sample) < size:
        x_try = np.random.uniform(xmin, xmax, size=chunk_size)
        y_try = np.random.uniform(ymin, ymax, size=chunk_size)
        
        z_try = np.random.uniform(0, zmax, size=chunk_size)
        zmax = max(zmax, np.max(z_try))

        x_sample += list(x_try[z_try < pdf(x_try, y_try)])
        y_sample += list(y_try[z_try < pdf(x_try, y_try)])
        
    # Transform the samples to a 1D numpy array
    x_sample = np.array(x_sample).flatten()
    y_sample = np.array(y_sample).flatten()
    # Return the correct number of samples
    return x_sample[:size], y_sample[:size]

def add_dictionaries_together(dictionary1, dictionary2):
    """Adds two dictionaries with the same keys together."""
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
        # Check if the item is an ndarray
        if isinstance(dictionary1[key], np.ndarray):
            dictionary[key] = np.concatenate((dictionary1[key], dictionary2[key]))
        # Check if the item is a nested dictionary
        elif isinstance(dictionary1[key], dict):
            dictionary[key] = add_dictionaries_together(
                dictionary1[key], dictionary2[key]
            )
        else:
            raise ValueError(
                "The dictionary contains an item which is neither an ndarray nor a dictionary."
            )
    return dictionary


def trim_dictionary(dictionary, size):
    """Filters an event dictionary to only contain the size."""
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
