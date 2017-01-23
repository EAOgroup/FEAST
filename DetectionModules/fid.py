"""
FID simulates a Flame Ionization Detector LDAR program. The module defines a new subclass of DetectionMethod.
The class has two methods: an initialization method and a detection method.
"""
from DetectionModules.abstract_detection_method import DetectionMethod
from GeneralClassesFunctions.simulation_functions import sample_wr
from GeneralClassesFunctions.simulation_functions import set_kwargs_attrs
import numpy as np
from DetectionModules import helper_functions


class FID(DetectionMethod):
    """
    Objects of this class contain the specifications for an FID detection program, as well as a method for detecting 
    leaks
    """

    def __init__(self, time, gas_field, **kwargs):
        """
        Inputs:
              time         a time object (Defined in feast_classes)
              gas_field    a gas_field object (Defined in feast_classes)
              kwargs       optional input dictionary that will override default parameters

        """
        DetectionMethod.__init__(self, time, gas_field)
        # -------------- Hardware variables --------------
        self.lifetime = 10 * 365  # days

        # -------------- Process Variables --------------
        self.survey_interval = 100  # days
        self.survey_speed = 150  # components/person-hour
        self.drive_speed = 15  # m/s
        self.setup_time = 0.1  # hours
        self.work_time = 5  # hours/day on average...5 hours per day on average corresponds to 35 hours per week.
        self.labor = 100  # dollars/hour

        set_kwargs_attrs(self, kwargs)
        # -------------- Set calculated parameters --------------
        self.survey_time = (gas_field.components_per_site / self.survey_speed +
                            gas_field.site_spacing / self.drive_speed / 3600 + self.setup_time)  # hours
        # time_factor accounts for the finite simulation size. The effective capital cost is
        # reduced in the simulation based on the ratio of the wells in the
        # simulation to the number of wells that a single FID could survey.
        self.time_factor = self.survey_time * gas_field.site_count / (self.survey_interval * self.work_time)

        # -------------- Financial Properties --------------
        self.capital_0 = 35000 * self.time_factor  # dollars (covers 5k for FID and 30k for truck)
        self.maintenance_0 = self.capital_0 * 0.1  # dollars/year

        self.capital = np.zeros(time.n_timesteps)
        helper_functions.replacement_cap(time, self)

        # maintenance costs are estimated as 10% of capital per year
        self.maintenance = [self.maintenance_0 * time.delta_t / 365, ] * time.n_timesteps  # $

        # survey_cost is the cost to survey all wells in the natural gas field
        self.survey_cost = self.labor * gas_field.site_count * self.survey_time

        # find_cost is the cost of searching for leaks
        for ind in range(0, time.n_timesteps):
            curr_time = ind * time.delta_t
            if curr_time % self.survey_interval < time.delta_t:
                self.find_cost[ind] = self.survey_cost

    def detection(self, time, gas_field, atm):
        """
        Inputs:
            time        an object of type Time (defined in feast_classes)
            gas_field   an object of type GasField (defined in feast_classes)
            atm         an object of type Atmosphere (defined in feast_classes)
                        Note: atm is not currently used by this method, but it is accepted as an argument for
                        consistency with other detection methods
        """
        self.null_detection(time, gas_field)
        # Take no action unless a survey interval is ending
        if time.current_time % self.survey_interval < time.delta_t:
            # Repair all leaks
            self.repair_cost[time.time_index] += \
                sum(sample_wr(gas_field.repair_cost_dist.repair_costs, len(self.leaks.flux)))
            self.leaks_found.extend(self.leaks.flux)
            self.leaks.delete_leaks(indexes_to_delete='all')
        return
