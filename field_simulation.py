from GeneralClassesFunctions import simulation_classes
# DetectionModules is used within an eval() function call.
import DetectionModules as Dm
from GeneralClassesFunctions.simulation_functions import new_leak_count, save_results
import copy


def field_simulation(gas_field=None, atm=None, input_leaks=None, dir_out='Results', time=None,
                     econ_set=None, tech_dict=None, detection_techs=None, display_status=True):
    """
    field_simulation generates a single realization of scenario. The scenario is defined by the input values.
    gas_field           a GasField object
    atm                 an Atmosphere object
    input_leaks         a list of leaks to be generated at each timestep
    dir_out             directory name in which to save results
    time                a Time object
    econ_set            a FinanceSettings object
    tech_dict           a dict of detection technology objects
    detection_techs     a list of detection technology identifying strings
    """
    # -------------- Define settings --------------
    # time defines parameters related to time in the model. Time units are days.
    if time is None:
        time = simulation_classes.Time()

    if gas_field is None:
        gas_field = simulation_classes.GasField()
    leak_list = copy.deepcopy(gas_field.initial_leaks)

    # Note: econ_settings are not used during the simulation, but are saved for use in post simulation data processing
    if econ_set is None:
        econ_set = simulation_classes.FinanceSettings()

    if atm is None:
        atm = simulation_classes.Atmosphere(time.n_timesteps)

    if detection_techs is None:
        detection_techs = ['fid.FID', 'null.Null', 'ir.AIR', 'ir.MIR', 'dd.DD']

    new_leaks, no_repair_leakage = [], []
    if tech_dict is None:
        tech_dict = dict()
        for tech in detection_techs:
            [tech_mod, tech_class] = tech.split(sep='.')
            tech_dict[tech_class] = eval('Dm.' + tech_mod.lower() + '.' + tech_class + '(time, gas_field)')
    elif 'Null' not in tech_dict:
        tech_dict['Null'] = Dm.null.Null(time, gas_field)
    # -------------- Run the simulation --------------
    n_leaks = []
    for time.time_index in range(0, time.n_timesteps):
        time.current_time += time.delta_t
        if display_status and time.current_time % int(time.end_time/10) < time.delta_t:
            print("Currently evaluating time step " + str(time.time_index) + " of " + str(time.n_timesteps))
        if input_leaks is None:
            new_leaks.append(gas_field.leak_size_maker(new_leak_count(time, gas_field), gas_field))
        else:
            new_leaks.append(input_leaks[time.time_index])
        no_repair_leakage.append(sum(leak_list.flux))
        leak_list.extend(new_leaks[-1])
        # Loop through each LDAR program:
        for tech_obj in tech_dict.values():
            n_leaks.append(tech_obj.leaks.n_leaks)
            tech_obj.leakage.append(sum(tech_obj.leaks.flux))
            tech_obj.detection(time, gas_field, atm)
            tech_obj.leaks.extend(new_leaks[-1])
    # -------------- Save results --------------
    results = simulation_classes.Results(time, gas_field, tech_dict, leak_list, no_repair_leakage, atm, econ_set,
                                         new_leaks, n_leaks)
    save_results(dir_out, results)
