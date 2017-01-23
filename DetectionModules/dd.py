"""
DD simulates a Distributed Detector LDAR program.
"""

from DetectionModules.abstract_detection_method import DetectionMethod
from GeneralClassesFunctions.simulation_functions import sample_wr
from GeneralClassesFunctions.simulation_functions import set_kwargs_attrs
from GeneralClassesFunctions import simulation_functions
import numpy as np
from DetectionModules import helper_functions


class DD(DetectionMethod):
    def __init__(self, time, gas_field, **kwargs):
        """ Inputs:
       gas_field    a gas_field object (Defined in feast_classes)
       time         a time object (Defined in feast_classes)
       kwargs       optional input dicitionary that will override default parameters
        """
        DetectionMethod.__init__(self, time, gas_field)
        # -------------- Hardware properties --------------
        # Sensitivity is the trigger sensitivity for the sniffer
        self.sensitivity = .01  # g/m^3 (multiply by 1500 to get ppm)
        # x, y, and z define the positions of the sensors
        # x is distance East from the center of the well pad
        # y is distance North from the center of the well pad
        self.x = [-3.5, 0, 3.5, 0]
        self.y = [-10, -10, -10, 10]
        self.z = [2, 3, 4, 4]
        self.lifetime = 5*365  # days

        # -------------- Process variables --------------
        self.repair_interval = 50  # days
        self.survey_speed = 500  # component/person-hour
        self.drive_speed = 15  # m/s
        self.setup_time = 0.5  # hours
        self.work_time = 5  # hours/day on average...5 hours per day on average corresponds to 35 hours per week.

        # -------------- Financial --------------
        self.sniffer_cost = 500  # $ per sniffer
        self.camera_capital = 90000  # $
        self.truck_capital = 30000  # $
        self.labor = 100  # $/hour
        self.maintenance_factor = 0.1

        # -------------- Override default parameters with kwargs --------------
        set_kwargs_attrs(self, kwargs)
        
        # Calculated parameters
        self.n_sniff = len(self.x)  # sniffers per well
        self.sniffer_capital = self.sniffer_cost*self.n_sniff  # $ per well
        self.location_time = (gas_field.components_per_site/self.survey_speed/self.n_sniff+gas_field.site_spacing /
                              self.drive_speed/3600+self.setup_time)/self.n_sniff  # hours to locate one leak
        self.time_factor = self.location_time*gas_field.site_count/(self.repair_interval*self.work_time)
        self.capital_0 = self.sniffer_capital*gas_field.site_count+self.time_factor * (self.camera_capital +
                                                                                       self.truck_capital)  # $
        self.maintenance0 = self.capital_0*self.maintenance_factor
        self.maintenance = [self.maintenance0*time.delta_t/365]*time.n_timesteps  # $/year
        self.location_cost = self.labor*self.location_time  # $/leak
            
        #  Calculate capital costs vector
        self.capital = np.zeros(time.n_timesteps)  # dollars
        helper_functions.replacement_cap(time, self)

        #  The finding np.costs must be calculated as leaks are found. The vector is
        #  initialized here.
        self.find_cost = [0]*time.n_timesteps

    def detection(self, time, gas_field, atm):
        """
        Determines which leaks are detected at each timestep
        Inputs:
        time:        object of type Time defining the time variables in the simulation
        gas_field:  object of type GasField defining gas field parameters
        atm:        object of type Atmosphere defining wind speed, direction and atmospheric stability
        """
        self.null_detection(time, gas_field)
        # Periodically repair all detected leaks
        if time.current_time % self.repair_interval < time.delta_t and time.current_time != 0:
            # countFound is a temporary variable to ease notation
            count_found = sum(self.leaks.leaks_detected)
            # account for the np.cost of locating the leaks
            self.find_cost[time.time_index] = self.location_cost*count_found
            # calculate the repair np.cost
            self.repair_cost[time.time_index] += sum(sample_wr(gas_field.repair_cost_dist.repair_costs,
                                                               int(count_found)))
            # update leaksFound and leakStruct
            if count_found > 0:
                # Identify indexes of leaks to delete
                ind = np.where(self.leaks.leaks_detected == 1)[0]
                self.leaks_found.extend(self.leaks.flux[ind])
                self.leaks.delete_leaks([ind])
        # Convert coordinates so that x is measured in the direction of the wind
        wind_dir_rad = np.radians(atm.wind_direction[time.time_index])
        detector_x, detector_y = helper_functions.coord_transform(wind_dir_rad, self.x, self.y)
        x_save, y_save = helper_functions.coord_transform(wind_dir_rad, self.leaks.x, self.leaks.y)
        # At each time step, iterate through each leak to see if it can be detected
        indexes = np.where(self.leaks.leaks_detected == 0)[0]
        # check if each detector detects the leak
        phi = simulation_functions.gauss_leak_model(self.x, self.y, self.z, self.leaks, atm, time.time_index,
                                                    index=indexes)
        for ind in range(0, self.n_sniff):
            found_indexes = np.where(phi[ind, :] > self.sensitivity)[0]
            self.leaks.leaks_detected[indexes[found_indexes]] = 1
        self.leaks.x, self.leaks.y, self.x, self.y = x_save, y_save, detector_x, detector_y
