import copy
import random
import numpy as np
from abc import ABCMeta
from GeneralClassesFunctions.simulation_functions import set_kwargs_attrs
from GeneralClassesFunctions.simulation_functions import sample_wr


class DetectionMethod(metaclass=ABCMeta):
    """
    DetectionMethod is an abstract super class that defines the form required for all detection methods
    """
    def __init__(self, time, gas_field, notes=None, **kwargs):
        """
        Inputs:
            time         a time object (Defined in simulation_classes)
            gas_field    a gas_field object (Defined in simulation_classes)
            notes        a description of the object created
            kwargs       optional input dict that will override default parameters
        """
        self.notes = notes
        self.null_repaired = []
        self.repair_cost = [0] * time.n_timesteps
        self.capital = [0] * time.n_timesteps
        self.find_cost = [0] * time.n_timesteps
        self.leaks_found = []
        # leaks will be updated throughout simulations. initial_leaks should remain constant, so a copy is needed.
        self.leaks = copy.deepcopy(gas_field.initial_leaks)
        self.count_found = []
        self.leakage = []
        # Set all attributes defined in kwargs, regardless of whether they already exist
        set_kwargs_attrs(self, kwargs, only_existing=True)
        
    def null_detection(self, time, gas_field):
        """
        Every detection method shares the same null_detection method defined here
        Inputs:
            time         a time object (Defined in simulation_classes)
            gas_field    a gas_field object (Defined in simulation_classes)    
        """
        n_repaired = np.random.poisson(len(self.leaks.flux) * gas_field.null_repair_rate * time.delta_t)
        index_repaired = random.sample(list(range(0, len(self.leaks.flux))), n_repaired)
        for ind in index_repaired:
            self.null_repaired.append(self.leaks.flux[ind])
        self.leaks.delete_leaks(index_repaired)
        self.repair_cost[time.time_index] += sum(sample_wr(gas_field.repair_cost_dist.repair_costs, n_repaired))
