"""
All functions that are required in the simulations and independent of input data and detection method are stored here.
"""

import numpy as np
import random
import os
import pickle


def sample_wr(data, n_samples=1):
    """
    Create a list of data by random sampling from a data set with replacement
    Inputs:
        data            List of data
        n_samples       number of samples to draw
    Return:
         sample         list of samples drawn
    """
    sample = []
    for ind in range(0, n_samples):
        sample.append(random.choice(data))
    return sample


def new_leak_count(time, gas_field):
    """
    Calculate the number of new leaks to generate
    Inputs:
        time        a Time object 
        gas_field   a GasField object
    Return:
         Number of new leaks
    """
    return np.random.poisson(gas_field.leak_production_rate * time.delta_t * gas_field.component_count)


def save_results(dir_out, results):
    """
    Save results to a file
    Inputs:
        dir_out             Name of output file to save
        results             A results object
    """
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    n_realization = len(os.listdir(dir_out))
    file_out = dir_out + '/realization' + str(n_realization) + '.p'
    pickle.dump(results, open(file_out, 'wb'))


def set_kwargs_attrs(obj_in, kwargs, only_existing=True):
    """
    Function for overwriting parameters with key word arguments
    Inputs:
        obj_in          Object with parameters to be updated
        kwargs          Dict containing new parameter values
        only_existing   Determines with new parameters can be created with kwargs
    """
    for key in kwargs.keys():
        # If only_existing is true, only set attributes that already exist
        if only_existing:
            if not hasattr(obj_in, key):
                raise ValueError("Tried to set invalid attribute. Class: ", type(obj_in), 'attempted attribute:', key)
        setattr(obj_in, key, kwargs[key])


# Calculate the concentration of methane at arbitrary points downwind of a leak source.
def gauss_leak_model(x, y, z, leak, atm, current_time, index=None):
    """
    Returns the concentration of methane at a position downwind of the leak
    Inputs:
        x, y, z     array of points at which to calculate the concentration [m]
                  x is measured parallel to the wind, z is vertical distance above the ground and y completes a right
                  hand coordinate system
        leak      Leak object containing numpy arrays of leak data
        atm       Atmosphere object containing atmosphere data
        index     index of leak to be analyzed from within 'leak.'
    Outputs:
        phi    array storing the concentration of methane at the desired
               coordinates [g/m^3]
    """
    # Definitions:
    # sigmay is the standard deviation of the concentration in the y direction.
    # sigmaz is the standard deviation of the concentration in the z direction.
    # sigmay and sigmaz are calculated as a function of x based on the Pasquill
    # stability category. They are linear fits to the 100 meter data on the
    # Pasquill Gifford curves.
    n_loc = len(x)
    # If an index (or list of indexes) is passed in, only consider those leaks. Otherwise, compute phi for all leaks
    if index is None:
        rng = np.array(range(0, leak.n_leaks))
        phi = np.zeros([n_loc, leak.n_leaks])
    elif type(index) is not np.ndarray:
        try:
            len(index)
            rng = np.array(index)
        except TypeError:
            rng = np.array([index])
        phi = np.zeros([n_loc, len(rng)])
    else:
        rng = index
        phi = np.zeros([n_loc, len(rng)])
    # Iterate through each location
    for loc in range(0, n_loc):
        if z[loc] < 0:
            continue
        normx = x[loc]-leak.x[rng]
        # poslocations limits the calculation to points downwind of the leak and above ground.
        poslocations = np.where(normx > 0)
        normx = normx[poslocations]
        temp_rng = rng[poslocations]
        # sigmay and sigmaz are arrays of dimension normx
        sigmaz = atm.l[current_time]*normx/(1+normx/atm.a[current_time])**atm.q[current_time]
        sigmay = atm.k[current_time]*normx/(1+normx/atm.a[current_time])**atm.p[current_time]
        #  a and b are variables used in later calculations to account for
        #  the buoyancy of the plume.
        b = -leak.z[temp_rng]-1.6 * leak.f_one_third[temp_rng] * normx**(2/3)/atm.wind_speed[current_time]
        a = z[loc]+b
        # Concentration Calculation
        fy = np.exp(-(y[loc]-leak.y[temp_rng])**2 / (2*sigmay**2))
        fz = np.exp(-a**2./(2*sigmaz**2))
        fg = np.exp(-(a-2*b)**2./(2*sigmaz**2))
        phi[loc, poslocations] = leak.flux[temp_rng]/(atm.wind_speed[current_time] * 2*np.pi*sigmaz*sigmay)*fy*(fz+fg)
        # g/m^3
    return phi
