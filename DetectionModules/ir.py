"""
IR simulates an infrared camera LDAR program. The module defines a new subclass of DetectionMethod.
The class has two methods: an initialization method and a detection method.
"""

from DetectionModules.abstract_detection_method import DetectionMethod
from DetectionModules import helper_functions
from GeneralClassesFunctions.simulation_functions import sample_wr
from GeneralClassesFunctions.simulation_functions import set_kwargs_attrs
import numpy as np
from scipy import integrate
import pickle

# Constants
H = 6.626e-34  # Planck's constant [J-s]
SIGMA = 5.67e-8  # Stefan-Boltzmann constant [W/m^2-K^4]
c = 3e8  # Speed of light [m/s]
K = 1.38e-23  # Boltzmann's constant [J/K]
METHANE_PNNL_FILE = 'InputData/DataObjectInstances/pnnl_methane.p'
# Constants
N_A = 6.023e23  # Avogadro's number
ROOT2PI = np.sqrt(2*np.pi)


class IR(DetectionMethod):
    """
    Objects of this class contain the specifications for an IR detection program, as well as a method for detecting
    leaks
    """
    def __init__(self, time, gas_field, **kwargs):
        """
        Inputs:
              time         a time object (Defined in feast_classes)
              gas_field    a gas_field object (Defined in feast_classes)
              kwargs       optional input dictionary that will override default parameters
        """
        # -------------- Define general parameters and methods --------------
        DetectionMethod.__init__(self, time, gas_field)
        # -------------- Process attributes --------------
        self.work_time = 5  # hours per day (maximum worktime per day per camera)
        self.labor = 100
        # -------------- Technology attributes --------------
        # clustering of pixels...for example, a comp_res_factor of 4 means that signals are calculated at 1/4 of the
        # pixels on the x and z axes, and 1/16 of all pixels
        self.comp_res_factor = 5
        self.pixels_x = 320
        self.pixels_z = 240
        self.spectrum_data = 'pnnl'
        self.fov1 = np.deg2rad(24)  # radians
        self.fov2 = np.deg2rad(18)  # radians
        self.detection_limit = 2000  # pixels
        self.orientation = [0, 0, 1]
        # minimum and maximum wavelength detected by the camera
        self.lambda1 = 3.2e-6  # [m]
        self.lambda2 = 3.4e-6  # [m]
        self.netd = 0.015  # Kelvin
        self.f_number = 1.5
        self.a_d = 9e-10  # pixel area (m^2)
        # -------------- Attributes to be calculated after custom settings from subclass --------------
        # Find cost is the cost of searching for leaks
        self.capital_0 = None
        self.distance = None
        self.default_pnnl = None
        self.survey_time = None
        self.survey_interval = None
        self.find_cost = np.zeros(time.n_timesteps)
        self.time_factor = None
        self.maintenance_0 = None
        self.xpoints = None
        self.zpoints = None
        self.point_ratio = None
        self.wx = None
        self.wz = None
        self.delta_wx = None
        self.delta_wz = None
        self.x = None
        self.y = None
        self.z = None
        self.maintenance = None
        self.survey_cost = None

    def set_attrs(self, time, gas_field):
        """
        Sets claculated attributes of the ir class.
        Inputs:
              time         a time object (Defined in feast_classes)
              gas_field    a gas_field object (Defined in feast_classes)
        """
        # timeFactor accounts for the finite simulation size. The effective capital cost is
        # reduced in the simulation based on the ratio of the wells in the
        # simulation to the number of wells that a single camera could survey.
        self.time_factor = self.survey_time * gas_field.site_count / (self.survey_interval * self.work_time)
        self.capital_0 *= self.time_factor
        self.maintenance_0 = self.capital_0*0.1  # dollars/year
        # Pixels are grouped into clusters of 4 for computational efficiency
        self.xpoints = int(self.pixels_x/self.comp_res_factor)
        self.zpoints = int(self.pixels_z/self.comp_res_factor)
        self.point_ratio = self.pixels_x / self.xpoints * self.pixels_z / self.zpoints
        # width of frame in x and z at the location of the plume
        self.wx = 2 * self.distance * np.tan(self.fov1/2)
        self.wz = 2 * self.distance * np.tan(self.fov2/2)
        # Distance between computation points at the location of the plume
        self.delta_wx = self.wx / self.xpoints
        self.delta_wz = self.wz / self.zpoints
        self.x = np.linspace(1, self.xpoints, self.xpoints)*self.delta_wx
        self.y = np.zeros(self.xpoints)
        self.z = np.linspace(0, self.zpoints//2, self.zpoints//2+1)*self.delta_wz

        # Set self.capital to capital_0 each time a camera reaches it's lifetime.
        self.capital = np.zeros(time.n_timesteps)  # dollars
        helper_functions.replacement_cap(time, self)
        # Maintenance costs are estimated as 10% of capital costs distributed evenly throughout each year
        self.maintenance = np.ones(time.n_timesteps)*self.maintenance_0 * time.delta_t / 365  # $
        
        # surveyCost is the cost to survey all wells in the natural gas field
        self.survey_cost = self.labor * gas_field.site_count * self.survey_time
        
        indexes = np.where(np.linspace(0, time.end_time, time.n_timesteps) % self.survey_interval < time.delta_t)
        self.find_cost[indexes] = self.survey_cost

    def detection(self, time, gas_field, atm):
        """
        detection assumes that all leaks are viewed from the side from the same camera position relative
        to the leak source. Reflection off the ground is neglected and plume rise is assumed to be small
        compared to the height of the camera.
        Inputs:
            time        an object of type Time (defined in simulation_classes)
            gas_field   an object of type GasField (defined in simulation_classes)
            atm         an object of type Atmosphere (defined in simulation_classes)
        """
        self.null_detection(time, gas_field)
        t_index = time.time_index
        if time.current_time % self.survey_interval < time.delta_t:
            # Define new leak object based on self.leaks with each leak at a standard height and sorted by flux
            n_l = self.leaks.n_leaks
            sort_perm = np.argsort(self.leaks.flux)
            flux = np.sort(self.leaks.flux)
            # calculate absorption values for the plume. If e_a, e_g, t_g t_plume and t_a are the same on consecutive
            # timesteps, methane_pnnl returns the same results without repeating any calculations
            e_g, e_a, t_g = atm.e_g[time.time_index], atm.e_a[time.time_index], atm.ground_temp[time.time_index]
            t_a, t_plume = atm.a_temp[time.time_index], gas_field.t_plume
            # Use the specified spectrum data
            if self.spectrum_data.lower() == 'pnnl':
                k_av, nep, tec = self.methane_pnnl(e_g, e_a, t_g, t_plume, t_a)
            else:
                raise ValueError("spectrum_data must be 'pnnl'")
            self.default_pnnl = [e_g, e_a, t_g, t_a, t_plume, k_av, nep, tec]
            high = n_l - 1
            low = -1
            # Binary search for the smallest detected leak
            counter = 0
            if self.pixel_detection(atm, flux[-1], k_av, nep, tec, t_index, e_a) < self.detection_limit:
                detect = []
            else:
                while high - low > 1:
                    counter += 1
                    leak_ind = low + (high-low)//2
                    total_pixels = self.pixel_detection(atm, flux[leak_ind], k_av, nep, tec, t_index, e_a)
                    if total_pixels < self.detection_limit:
                        low = leak_ind
                    else:
                        high = leak_ind
                    if counter > 100:
                        raise ValueError('Aborted due to over 100 iterations in ir.detection')
                detect = sort_perm[high:]
            self.repair_cost[time.time_index] += \
                sum(sample_wr(gas_field.repair_cost_dist.repair_costs, len(detect)))
            self.leaks_found.extend(self.leaks.flux[detect])
            # Delete found leaks
            self.leaks.delete_leaks(detect)

    def pixel_detection(self, atm, flux, k_av, nep, tec, t_index, e_a):
        """
        Calculate the number of pixels that detect a plume
        Inputs:
            atm     atmosphere object
            flux    flux from one leak
            k_av    average absorption coefficient
            nep     noise equivalent power
            tec     temperature emissivity coefficient
            t_index current time index
            e_a     atmospheric emissivity
        Return:
            total_pixels    The number of pixels that detect the leak
        """
        pixels = 0
        pixels0 = 0
        # Store the camera orientation as a bool (This structure allows for other orientations to be supported in future
        # releases)
        if self.orientation[1] == 1:
            isy = True
        elif self.orientation[2] == 1:
            isy = False
        else:
            raise ValueError("ir.orientation must be [0,0,1] or [0,1,0]")
        for j in range(0, self.xpoints):
            # Only the width of the plume perpendicular to the imaging direction matters
            if isy:
                # Calculate sigma_z
                sigma = atm.l[t_index]*self.x[j]/(1+self.x[j]/atm.a[t_index])**atm.q[t_index]
            else:
                # calculate sigma_y
                sigma = atm.k[t_index]*self.x[j]/(1+self.x[j]/atm.a[t_index])**atm.p[t_index]
            for m in range(0, self.zpoints//2+1):
                if m == 0:
                    f = 1
                else:
                    f = np.exp(-self.z[m]**2/(2*sigma**2))
                final = flux * f / (ROOT2PI * atm.wind_speed[t_index] * sigma)
                attn = final * k_av * 0.434 * N_A * 1e-4
                factor = 1-10**(-attn)
                contrast = abs(tec) * factor * (1-e_a)
                if contrast >= nep:
                    pixels += 1
                    if m == 0:
                        pixels0 += 1
        total_pixels = (2*pixels-pixels0) * self.point_ratio
        return total_pixels

    def methane_pnnl(self, e_g, e_a, t_g, t_plume, t_a):
        """
        Calculates the light absorption properties of a methane plume using the Pacific Northwest National Lab library
        Inputs:
            ir      an ir detection method
            e_g     emissivity of the ground
            e_a     emissivity of the atmosphere
            t_g     temperature of the ground
            t_plume temperature of the plume
            t_a     temperature of the atmosphere
        Return:
            kav     average k value over the range of wavelengths from 3.2 to 3.4 microns for methane
            nep     noise equivalent power (used to determine the minimum detectable change in power) [W]
            tec     temperature emissivity contrast [W]
        """
        if self.default_pnnl is not None and self.default_pnnl[:5] == [e_g, e_a, t_g, t_plume, t_a]:
            return self.default_pnnl[5], self.default_pnnl[6], self.default_pnnl[7]
        with open(METHANE_PNNL_FILE, 'rb') as f:
            pnnl = pickle.load(f)
        kav = np.average(pnnl.km)
        nep, tec = self.bb_rad_params(t_g, t_plume, t_a, e_g, e_a)
        return kav, nep, tec

    def bb_rad_params(self, t_g, t_plume, t_a, e_g, e_a):
        """
        Returns black body radiation parameters related to the plume
        Inputs:
            t_g     temperature of the ground
            t_plume temperature of the plume
            t_a     temperature of the atmosphere
            e_g     emissivity of the ground
            e_a     emissivity of the atmosphere
        Return:
            nep     noise equivalent power (used to determine the minimum detectable change in power) [W]
            tec     temperature emissivity contrast [W]
        """

        # -------------- Calculate the noise equivalent power (NEP) --------------
        # based on environmental conditions and camera properties
        # w1g and w2g are the frequency limits over which the camera is sensitive (multiplied by H/K t_g to make them
        # non-dimensional)
        w1g = H * c / (self.lambda2 * K * t_g)
        w2g = H * c / (self.lambda1 * K * t_g)
        # The product n1 * y1 is the integral of the derivative of blackbody radiation with respect to temperature
        # over the range of frequencies that the camera is sensitive to multiplied by a constant.
        n1 = 2 * np.pi * K**4 * t_g**3 / (H**3 * c**2)
        temp_y1 = -np.exp(-w1g) * (720+720 * w1g + 360 * w1g**2+120 * w1g**3 + 30 * w1g**4 + 6 * w1g**5 + w1g**6)
        temp_y2 = -np.exp(-w2g) * (720+720 * w2g + 360 * w2g**2+120 * w2g**3 + 30 * w2g**4 + 6 * w2g**5 + w2g**6)
        y1 = temp_y2 - temp_y1
        y = y1 * n1
        # Noise equivalent power
        nep = y * self.netd * self.a_d / (4 * self.f_number**2)

        # Calculate the power incident on a pixel from a blackbody at the ground temp, plume temp and atmosphere temp.
        ppixelg = self.pixel_power(t_g)
        ppixelp = self.pixel_power(t_plume)
        ppixela = self.pixel_power(t_a)
        # tec is the signal in a pixel caused by an infinitely thick plume in an atmosphere with an emissivity of zero. 
        # tec is a useful value because the signal due to a finite plume in a real atmosphere is proportional to tec. 
        tec = ppixelp - e_g * ppixelg - e_a * (1-e_g) * ppixela

        return nep, tec

    def pixel_power(self, temp):
        """
        Calculate the the power incident on a pixel from an infinite blackbody emitter at a given temperature.
        Inputs:
            temp    Temperature of the emitter (K)
        Return:
            pixel_power   power incident on the pixel (W)
        """
        # Calculate the non-dimensional frequency limits of the sensor
        w1 = H * c / (self.lambda2 * K * temp)
        w2 = H * c / (self.lambda1 * K * temp)
        # Integrate the blackbody radiation over the frequency range
        temp_int = integrate.quad(lambda x: x**3 / (np.exp(x)-1), w1, w2)
        # calculate the power incident on one camera pixel
        frac = temp_int[0] / (np.pi**4 / 15)
        sblaw = SIGMA * temp**4 * self.a_d
        power = (4/np.pi) * sblaw * np.tan(self.fov1 / 2) * np.tan(self.fov2 / 2)
        pixel_power = power * frac
        return pixel_power


class MIR(IR):
    """
    Class for a manually operated infrared detection method. The plume is assumed to be observed from the side in the
    MIR class and from above in the AIR class. Otherwise, the two classes are nearly identical with different default
    values.
    """
    def __init__(self, time, gas_field, **kwargs):
        IR.__init__(self, time, gas_field)
        # -------------- Process details --------------
        self.survey_interval = 100  # days
        self.survey_speed = 500  # component/person-hour
        self.drive_speed = 15  # m/s
        self.setup_time = 0.5  # hours
        self.work_time = 5  # hours/day on average...5 hours per day on average corresponds to 35 hours per week.
        # -------------- Detection technology properties --------------
        self.distance = 2  # meters from leak
        # capital_0 is spread over the maximum possible number of wells
        self.capital_0 = 120000  # dollars (covers 90k for camera and 30k for truck)
        self.labor = 100  # dollars/hour
        self.lifetime = 10*365  # days
        self.orientation = np.array([0, 1, 0])  # Camera orientation relative to the wind (sign is irrelevant)
        self.spectrum_data = 'pnnl'  # database to use for spectral properties
        # -------------- Override default parameters with kwargs --------------
        set_kwargs_attrs(self, kwargs)
        # -------------- Calculated parameters --------------
        self.survey_time = (gas_field.components_per_site/self.survey_speed+gas_field.site_spacing/self.drive_speed /
                            3600 + self.setup_time)  # hours
        self.set_attrs(time, gas_field)


class AIR(IR):
    """
    Class for a manually operated infrared detection method. The plume is assumed to be observed from above in the AIR
    class and from the side in the MIR class . Otherwise,the two classes are nearly identical with different default
    values.
    """
    def __init__(self, time, gas_field, **kwargs):
        IR.__init__(self, time, gas_field, **kwargs)
        # -------------- Process details --------------
        self.survey_interval = 14  # days
        self.survey_speed = 5  # m/s
        self.drive_speed = 15  # m/s
        self.setup_factor = 1.3  # allocates 30% of flight time to landing, refueling, and launching quadcopter
        self.work_time = 5  # hours/day on average...5 hours per day on average corresponds to 35 hours per week.
        # -------------- Detection technology properties --------------
        self.distance = 20  # meters
        self.detection_limit = 2000  # pixels
        # capital_0 is spread over the maximum possible number of wells
        self.capital_0 = 193000  # dollars (covers 90k for camera, 30k for truck, 50K for UAV, 10k for optics and 2500
        # for computers)
        self.labor = 100  # dollars/hour
        self.lifetime = 10*365  # days
        # camera orientation relative to the wind (x is with the wind, sign is irrelevant)
        self.orientation = np.array([0, 0, 1])
        self.spectrum_data = 'pnnl'  # database to use for spectral properties
        # -------------- Override default parameters with kwargs --------------
        set_kwargs_attrs(self, kwargs)
        # -------------- Calculated Parameters --------------
        # The camera is assumed to be oriented with its z axis perpendicular to the wind direction
        self.image_width = self.distance * np.tan(self.fov2/2) * 2
        drive_time = gas_field.site_spacing / self.drive_speed
        site_time = gas_field.well_length**2 / (self.image_width * self.survey_speed / 2)
        self.survey_time = (site_time + drive_time) * self.setup_factor / 3600
        self.set_attrs(time, gas_field)
