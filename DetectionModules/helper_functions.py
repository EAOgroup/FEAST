"""
These functions can be called by DetectionModules to aid in calculations
"""
import copy
import numpy as np


def coord_transform(theta, x, y):
    """
    Transforms coordinates so that x is measured in the direction of the wind.
    NOTE: coord_transform updates the x and y values passed to the function. The original values are returned as
        x_save and y_save
    Inputs:
        x: an array of floats
        y: an array of floats
        theta: an angle given in radians
    Return:
         x_save: original x leak coordinates
         y_xave: original y leak coordinates
    """
    # Define sniffer and leak positions with respect to the direction of the wind:

    sin_theta, cos_theta = np.sin(theta), np.cos(theta)
    x_save, y_save = copy.copy(x), copy.copy(y)
    for ind in range(0, len(x)):
        x[ind] = x_save[ind]*sin_theta + y_save[ind]*cos_theta
        y[ind] = -x_save[ind]*cos_theta + y_save[ind]*sin_theta
    return x_save, y_save


def replacement_cap(time, tech):
    """
    Enters the replacement value of the technology into tech.capital at the end of its lifetime.
    Inputs:
        time        Object defining time parameters for the simulation
        tech        LDAR technology object
    """
    lifetimes = 0
    while lifetimes < time.end_time:
        tech.capital[max(1, round(lifetimes / time.delta_t))] = tech.capital_0
        lifetimes += tech.lifetime
