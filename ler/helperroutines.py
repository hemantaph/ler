''' Helper functions for various tasks used in LER, for example combining dictionaries together. Really this is a place for routines which don't seem to fit into anywhere else. '''

import numpy as np

chunk_size = 10000

def trim_dictionary(dictionary, size):
    ''' Filters an event dictionary to only contain the size. '''
    for key in dictionary.keys():
        # Check if the item is an ndarray
        if isinstance(dictionary[key], np.ndarray):
            dictionary[key] = dictionary[key][:size] # Trim the array
        # Check if the item is a nested dictionary
        elif isinstance(dictionary[key], dict):
            dictionary[key] = trim_dictionary(dictionary[key], size) # Trim the nested dictionary
        else:
            raise ValueError("The dictionary contains an item which is neither an ndarray nor a dictionary.")
    return dictionary

def trim_dictionary_by_indices(dictionary, indices):
    ''' Filters an event dictionary to only contain the indices. '''
    for key in dictionary.keys():
        # Check if the item is an ndarray
        if isinstance(dictionary[key], np.ndarray):
            dictionary[key] = dictionary[key][indices] # Trim the array
        # Check if the item is a nested dictionary
        elif isinstance(dictionary[key], dict):
            dictionary[key] = trim_dictionary_by_indices(dictionary[key], indices) # Trim the nested dictionary
        else:
            raise ValueError("The dictionary contains an item which is neither an ndarray nor a dictionary.")
    return dictionary

def add_dictionaries_together(dictionary1, dictionary2):
    ''' Adds two dictionaries with the same keys together. '''
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
            dictionary[key] = add_dictionaries_together(dictionary1[key], dictionary2[key])
        else:
            raise ValueError("The dictionary contains an item which is neither an ndarray nor a dictionary.")
    return dictionary

# Helper function for rejection sampling
def rejection_sample(pdf, xmin, xmax, size=100):
    ''' Helper function for rejection sampling from a pdf with maximum and minimum arguments.
    Input parameters:
        pdf: the pdf to sample from
        xmin: the minimum argument of the pdf
        xmax: the maximum argument of the pdf
        size: the number of samples to draw
    Output:
        samples: the samples drawn from the pdf
    '''
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


# Helpelr function
def combine_lens_parameter_dictionaries(lensed_parameters, lensed_parameters_draw, idx, n_images):
    ''' Adds lensed_parameters_draw to lensed_parameters dictionary for selected events idx and for n_images.

        Input parameters:
        lensed_parameters (dict): Dictionary of lensed parameters
        lensed_parameters_draw (dict): Dictionary of lensed parameters to be added to lensed_parameters
        idx (int): Index of the events to be added to lensed_parameters
        n_images (int): Number of images

        Output parameters:
        lensed_parameters (dict): Dictionary of lensed parameters
    '''
    for key in lensed_parameters_draw.keys():
        item = lensed_parameters_draw[key]
        # Check if the key exists, if not, create it
        if key not in lensed_parameters.keys():
            # Check if  a nested dictionary. Add the new events to the dictionary for each key and create new one if it does not exist
            if isinstance(item, dict):
                lensed_parameters[key] = {}
                lensed_parameters[key] = combine_lens_parameter_dictionaries(lensed_parameters[key], item, idx, n_images)
            else:
                lensed_parameters[key] = item[idx]
        # Check if ndarray, and if so, add the value to the list
        elif isinstance(lensed_parameters[key], np.ndarray):
            lensed_parameters[key] = np.concatenate((lensed_parameters[key], item[idx]))
        # Check if a nested dictionary. Add the new events to the dictionary for each key and create new one if it does not exist
        else:
            # Nested dictionary. Add the new events to the dictionary for each key and create new one if it does not exist
            lensed_parameters[key] = combine_lens_parameter_dictionaries(lensed_parameters[key], item, idx, n_images)
    return lensed_parameters

def impose_snr_cut(lensed_parameters, snr_threshold, n_images):
    # Get rid of the events that do not have an snr threshold for the n_images images
    snrs = lensed_parameters['snr_opt_snr_net'] # Shape (nevents, n_max_images)
    # Take out the snrs that are not used
    snrs = snrs[:,:n_images]
    # Get the indices of the events that have an snr threshold for each image
    idx = np.where(np.all(snrs > snr_threshold, axis=1))[0]
    # Get rid of the events that do not have an snr threshold for each image
    lensed_parameters = trim_dictionary_by_indices(lensed_parameters, idx)
    return lensed_parameters

def save_dictionary_to_numpy_txt_file(detectable_lensed_event_parameters, fname= 'detectable_lensed_event_parameters.txt' ):
    ''' Saves a dictionary to a numpy txt file.

    Parameters:
        detectable_lensed_event_parameters (dict): Dictionary to be saved
        fname (str): Name of the file to be saved

    Example:
    from ler import helperroutines as hr
    # Save the detectable lensed event parameters dictionary as numpy txt file
    hr.save_dictionary_to_numpy_txt_file(detectable_lensed_event_parameters, fname= 'detectable_lensed_event_parameters.txt' )
    # Load the detectable lensed event parameters dictionary
    data = np.genfromtxt('detectable_lensed_4_image_event_parameters.txt', names=True)
    names = data.dtype.names
    '''
    # Save the detectable lensed event parameters dictionary as numpy txt file
    keys = detectable_lensed_event_parameters.keys()
    data = []
    names = []
    for key in keys:
        # Check if 1D or 2D numpy array
        if len(detectable_lensed_event_parameters[key].shape) == 1:
            data.append(detectable_lensed_event_parameters[key])
            names.append(key)
        else:
            for i in range(detectable_lensed_event_parameters[key].shape[1]):
                data.append(detectable_lensed_event_parameters[key][:,i])
                names.append(key + '_' + str(i))
    np.savetxt(fname, np.transpose(data), header=' '.join(names))

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)
