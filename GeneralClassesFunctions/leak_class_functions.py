"""
Leak data, leak distribution properties, and leak objects are created in this module
"""
import pickle
import numpy as np
import random

# Constants:
g = 9.8  # g is the strength of gravity [m/s^2]
RHO_AIR = 1225  # density of air [g/m^3]
RHO_METHANE = 681  # density of methane at atmospheric pressure [g/m^3]


class Leak:
    """
        Stores a list of leaks
    """
    def __init__(self, x_pos=(), y_pos=(), z_pos=(), flux=(), f_one_third=(), leaks_detected=(), capacity=0):
        """
        Inputs:
        x_pos               East-west position of leak with respect to the center of the well (m)
        y_pos               North-south position of leak with respect to the center o the well (m)
        z_pos               Altitude of leak with respect to the ground (m)
        flux                leak size (g/s)
        f_one_third         A plume dispersion factor
        leaks_detected      Binary value to save whether the leak has been detected or not (1 if detected, 0 otherwise)
        capacity            Expected total number of leaks to be stored in this instance of Leak (allows for faster
                            extend method)
        """
        if f_one_third is () and flux is not ():
            f_one_third = f_one_third_calc(flux)
        if leaks_detected is () and flux is not ():
            leaks_detected = np.zeros(len(flux))
        try:
            length_in = len(x_pos)
            if not len(x_pos) == len(y_pos) == len(z_pos) == len(flux) == len(f_one_third):
                raise ValueError("x_pos, y_pos, z_pos, flux, f_one_third and leaks_detected must be equal length")
        except TypeError:
            length_in = 1

        if capacity == 0:
            self.x = np.array(x_pos)
            self.y = np.array(y_pos)
            self.z = np.array(z_pos)
            self.flux = np.array(flux)
            self.f_one_third = np.array(f_one_third)
            self.leaks_detected = np.array(leaks_detected) if leaks_detected != () else np.zeros(length_in)
        else:
            self.x = np.zeros(capacity)
            self.y = np.zeros(capacity)
            self.z = np.zeros(capacity)
            self.flux = np.zeros(capacity)
            self.f_one_third = np.zeros(capacity)
            self.leaks_detected = np.zeros(capacity)
            self.x[0:length_in] = x_pos
            self.y[0:length_in] = y_pos
            self.z[0:length_in] = z_pos
            self.flux[0:length_in] = flux
            self.f_one_third[0:length_in] = f_one_third
            self.leaks_detected[0:length_in] = np.array(leaks_detected) if leaks_detected != () else np.zeros(length_in)
        self.n_leaks = length_in

    def extend(self, leak_obj_in):
        """
        Add a new leak
        Inputs:
            leak_obj_in     a Leak object
        """
        if len(self.x) - leak_obj_in.n_leaks - self.n_leaks >= 0:
            self.x[self.n_leaks: self.n_leaks + leak_obj_in.n_leaks] = leak_obj_in.x
            self.y[self.n_leaks: self.n_leaks + leak_obj_in.n_leaks] = leak_obj_in.y
            self.z[self.n_leaks: self.n_leaks + leak_obj_in.n_leaks] = leak_obj_in.z
            self.flux[self.n_leaks: self.n_leaks + leak_obj_in.n_leaks] = leak_obj_in.flux
            self.f_one_third[self.n_leaks: self.n_leaks + leak_obj_in.n_leaks] = leak_obj_in.f_one_third
            self.leaks_detected[self.n_leaks: self.n_leaks + leak_obj_in.n_leaks] = leak_obj_in.leaks_detected
        else:
            self.x = np.append(self.x, leak_obj_in.x)
            self.y = np.append(self.y, leak_obj_in.y)
            self.z = np.append(self.z, leak_obj_in.z)
            self.flux = np.append(self.flux, leak_obj_in.flux)
            self.f_one_third = np.append(self.f_one_third, leak_obj_in.f_one_third)
            self.leaks_detected = np.append(self.leaks_detected, leak_obj_in.leaks_detected)
        self.n_leaks += leak_obj_in.n_leaks

    def delete_leaks(self, indexes_to_delete):
        """
        Delete all parameters associated with leaks at indexes 'indexes_to_delete'
        indexes_to_delete       A list of leak indexes to delete, or the string 'all'
        """
        if type(indexes_to_delete) is str:
            if indexes_to_delete == 'all':
                indexes_to_delete = list(range(0, self.n_leaks))
            else:
                raise ValueError('indexes_to_delete must be a scalar, an array or the str "all"')
        self.x = np.delete(self.x, indexes_to_delete)
        self.y = np.delete(self.y, indexes_to_delete)
        self.z = np.delete(self.z, indexes_to_delete)
        self.flux = np.delete(self.flux, indexes_to_delete)
        self.f_one_third = np.delete(self.f_one_third, indexes_to_delete)
        self.leaks_detected = np.delete(self.leaks_detected, indexes_to_delete)
        try:
            self.n_leaks -= len(indexes_to_delete)
        except TypeError:
            self.n_leaks -= 1


def leak_objects_generator(dist_type, leak_data_path):
    """
    leak_objects is a parent function that will be called to initialize gas fields
    Inputs:
        dist_type           Type of leak distribution to be used
        leak_data_path      Path to a leak data file
    """
    leak_data_path = 'InputData/DataObjectInstances/' + leak_data_path
    leak_data = pickle.load(open(leak_data_path, 'rb'))

    # Define leak params and leak_size_maker based on the leak distribution type
    if dist_type == 'bootstrap':
        leak_params = leak_data
        leak_size_maker = bootstrap_leak_maker
    else:
        raise NameError('Leak distribution type unsupported in GasField')

    # Number of leaking components at each well (Poisson distribution)
    detection_types, leaks_per_well = leak_data.leak_sizes.keys(), 0
    # This sums the leaks found per well by every method in leak_data
    for key in detection_types:
        leaks_per_well += len(leak_data.leak_sizes[key]) / leak_data.well_counts[key]
    return leak_size_maker, leak_params, leaks_per_well


def bootstrap_leak_maker(n_leaks_in, gas_field):
    """
    Create leaks using a bootstrap method
    n_leaks_in              number of leaks to generate
    gas_field               a GasField object
    """
    leak_params = gas_field.leak_params
    if gas_field.dist_type.lower() == 'bootstrap':
        detection_methods = list(leak_params.leak_sizes.keys())
        flux = []
        round_err, leaks_per_well, counter = [], [], -1
        # Calculate leaks per well identified with each detection method stored in leak_params
        for method in detection_methods:
            n_leaks = len(leak_params.leak_sizes[method])
            n_wells = leak_params.well_counts[method]
            leaks_per_well.append(n_leaks/n_wells)
        # Generate the appropriate number of leaks from the distribution associated with each detection method
        for method in detection_methods:
            counter += 1
            n_leaks_key = leaks_per_well[counter] / sum(leaks_per_well) * n_leaks_in
            for ind in range(0, int(n_leaks_key)):
                flux.append(random.choice(leak_params.leak_sizes[method]))
            round_err.append(n_leaks_key % 1)
        # Add leaks omitted due to inability to add fractional leaks
        # The "round" function in the following line is intended to eliminate floating point errors. 
        chooser = np.random.uniform(0, sum(round_err), round(sum(round_err)))
        error_intervals = np.cumsum(round_err)
        for choose in chooser:
            ind = 0
            # Add a leak from the appropriate detection method
            while choose > error_intervals[ind]:
                ind += 1
            flux.append(random.choice(leak_params.leak_sizes[detection_methods[ind]]))
        random.shuffle(flux)
    else:
        raise NameError('Leak size distribution type is unsupported')

    x = list(np.random.uniform(-gas_field.well_length/2, gas_field.well_length/2, n_leaks_in))
    y = list(np.random.uniform(-gas_field.well_length/2, gas_field.well_length/2, n_leaks_in))
    z = list(np.random.uniform(0, gas_field.h0_max, n_leaks_in))
    f_one_third = f_one_third_calc(flux)
    return Leak(x_pos=x, y_pos=y, z_pos=z, flux=flux, f_one_third=f_one_third, leaks_detected=[0]*len(x))


def f_one_third_calc(flux):
    """
    Computes the f_one_third leak dispersion parameter given a leak flux
    Inputs:
        flux        Leakage rate. Must be an iterable list or array [g/s]
    Return:
        f_one_third     Leak dispersion parameter
    """
    f_factor = g / np.pi*(1/RHO_METHANE-1/RHO_AIR)  # Buoyancy flux [m^4/s^3]
    f_one_third = []
    for item in flux:
        f_one_third.append((f_factor * item)**(1/3))
    return f_one_third
