import matplotlib.pyplot as plt
from pickle import load
from os import listdir
from os.path import isfile, join
import numpy as np
from . import results_analysis_functions

# Stanford bright color palette (https://identity.stanford.edu/overview/color)
color_set = [[140, 21, 21], [0, 152, 219], [0, 155, 118],  [178, 111, 22], [83, 40, 79],  [0, 0, 0], [77, 79, 83]]
for i in range(0, len(color_set)):
    for j in range(0, len(color_set[i])):
        color_set[i][j] /= 255


def display_save(display, save, file_out):
    """
    display_save can display and/or save a plot
    Inputs
        display         choose whether to display the plot or not
        save            choose whether to save the plot or not
        file_out        define a path at which to save the file
    """
    if save:
        if file_out is None:
            file_out = 'Figures/time_series' + str(len(listdir('Figures'))) + '.pdf'
        plt.savefig(file_out)

    if display:
        plt.show()
    else:
        plt.close()
    plt.ioff()


def display_settings(interactive=True, inline=True):
    """
    If running in IPYTHON, display_settings uses magic commands to plot interactively
    interactive         determines whether the plot is interactive
    """
    if interactive:
        # display settings depend on the python environment being used
        try:
            __IPYTHON__
            ipy = get_ipython()
            if inline is True:
                try:
                    ipy.magic("matplotlib inline")
                except:
                    ipy.magic("matplotlib")
            else:
                ipy.magic("matplotlib")
        except:
            plt.ion()


def time_series(results_file, display=True, save=False, file_out=None, interactive=True):
    """
    Display a time series of leakage from each detection method in a results file
    Inputs:
        results_file    path to a results file
        display         choose whether to display the plot or not
        save            choose whether to save the plot or not
        file_out        define a path at which to save the file
        interactive     choose whether to make the plot interactive
    """
    display_settings(interactive)
    results = load(open(results_file, 'rb'))
    tech_dict = results.tech_dict
    fig = plt.figure()
    ax = fig.add_subplot(111)
    counter = -1
    for tech in tech_dict.keys():
        counter += 1
        tech_dict[tech].line = ax.plot(np.array(range(0, results.time.n_timesteps)) * results.time.delta_t,
                                       np.array(tech_dict[tech].leakage)/results.gas_field.site_count, label=tech,
                                       color=color_set[counter])
    plt.xlabel('Time [days]')
    plt.ylabel('Leakage per well [g/s]')
    plt.legend()
    display_save(display, save, file_out)


def summary_plotter(directory, display=True, save=False, file_out=None, interactive=True):
    """
    The NPV for each realization stored in 'directory is calculated and displayed in a stacked bar chart. Each component
    of the NPV is displayed separately in the chart.
    Inputs:
        directory    path to a directory containing results files
        display      boolean value that determines whether or not the plot is displayed
        save         boolean value that determines whether the plot is saved
        file_out     file_out path to a file in which to save the figure if save is True
        interactive  determines whether the plot is produced interactively
    """
    # load data
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
    sample = load(open(directory + '/' + files[0], 'rb'))
    _, npv_in, _, leakage = results_analysis_functions.results_analysis(directory)
    n_realizations = len(npv_in['Total'][0, :])
    tech_keys = []
    # Get a list of indexes that do not have data related to the null detection method
    for key in sample.tech_dict.keys():
        if key != 'Null':
            tech_keys.append(key)
    n_tech = len(sample.tech_dict)
    npv_std = np.std(npv_in['Total'], 1)
    totals = np.average(npv_in['Total'], 1)
    npv = dict()
    for key in npv_in.keys():
        npv[key] = np.sum(npv_in[key], 1)/n_realizations

    # create figure
    display_settings(interactive)
    plt.figure()
    width = 0.35
    ind = np.linspace(0, n_tech-2, n_tech-1)
    counter = 0
    neg_bottoms = [0]*(n_tech-1)
    pos_bottoms = [0]*(n_tech-1)
    key_list = list(npv.keys())
    key_list.sort()
    for key in key_list:
        if key == 'Total' or key == 'Gas':
            continue
        else:
            for k in ind:
                m = int(k)
                if npv[key][m] >= 0:
                    if m == 0:
                        plt.bar(m, -npv[key][m], width, bottom=neg_bottoms[m], color=color_set[counter], label=key)
                    else:
                        plt.bar(m, -npv[key][m], width, bottom=neg_bottoms[m], color=color_set[counter])
                    neg_bottoms[m] -= npv[key][m]
                else:
                    if m == 0:
                        plt.bar(m, -npv[key][m], width, bottom=pos_bottoms[m], color=color_set[counter], label=key)
                    else:
                        plt.bar(m, -npv[key][m], width, bottom=pos_bottoms[m], color=color_set[counter])
                    pos_bottoms[m] -= npv[key][m]
        counter += 1
    plt.bar(ind, npv['Gas'], width, bottom=pos_bottoms, color=color_set[counter], label='Gas')
    npv_std /= np.sqrt(n_realizations)
    plt.errorbar(ind + width / 2, totals, yerr=npv_std, linewidth=5, fmt='o', ms=1, label='Total NPV')
    plt.ylabel("NPV (k$/well)")
    plt.xticks(np.linspace(width/2, len(tech_keys) - 1 + width/2, len(tech_keys)), tech_keys)
    ax = plt.gca()
    ax.legend(bbox_to_anchor=(1.4, 1))
    display_save(display, save, file_out)


def hist_plotter(directory, display=True, save=False, file_out=None, interactive=True, inline=False, bins=None):
    """
    Plots histogram of leaks found by each technology based on results in the folder 'Directory'
    Inputs:
        directory     Path to a results folder to be analyzed
        display      boolean value that determines whether or not the plot is displayed
        save         boolean value that determines whether the plot is saved
        file_out     file_out path to a file in which to save the figure if save is True
        interactive  determines whether the plot is produced interactively
        inline       determines whether the plot is displayed inline if possible
        bins         defines the bins used in the histograms
    Return: None
    """

    display_settings(interactive, inline=inline)
    _, _, found, _ = results_analysis_functions.results_analysis(directory)
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
    leaks_made = []
    sample = load(open(directory + '/' + files[0], 'rb'))
    techs = sample.tech_dict.keys()
    tech_leaks_found = {}
    well_count = 0
    # Load data from files
    for tech in techs:
        tech_leaks_found[tech] = []
    for file in files:
        with open(directory + '/' + file, 'rb') as f:
            sample = load(f)
            leaks_made.extend(sample.gas_field.initial_leaks.flux)
            well_count += sample.gas_field.site_count
            for leak_list in sample.new_leaks:
                leaks_made.extend(leak_list.flux)
            for tech in techs:
                tech_leaks_found[tech].extend(sample.tech_dict[tech].leaks_found)
    leaks_made = np.array(leaks_made)
    n_bins = 30
    if bins is None:
        bins = np.linspace(0, max(leaks_made), n_bins+1)
    cumsum_points = 100
    counts, bins, patches = plt.hist(leaks_made, bins, normed=0, alpha=0.75)
    plt.close()
    # Calculative the cumulative fraction of leakage above a range of leak fluxes
    cumsum = np.zeros(cumsum_points)
    cumsumx = np.linspace(0, max(bins), cumsum_points)
    for ind in reversed(range(0, cumsum_points-1)):
        cumsum[ind] = sum(leaks_made[leaks_made > cumsumx[ind]])
    cumsum /= cumsum[0]
    centers = 0.5*(bins[1:] + bins[:-1])
    ind = 0
    delta = bins[1] - bins[0]
    # Create a separate plot for each technology
    for tech in techs:
        if tech.lower() == 'null':
            continue
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        counts_temp, bins_temp, patches = ax1.hist(tech_leaks_found[tech], bins, normed=0, alpha=0.75)
        ax1.cla()
        counts_temp /= well_count
        ax1.bar(bins_temp[0:n_bins], counts_temp, width=delta, color=color_set[ind+1], label="Leaks found")
        ax1.plot(centers, counts/well_count, 'o', color=color_set[0], label="Leaks generated")
        ax2.plot(cumsumx, cumsum, label="Cumulative leakage fraction")
        plt.title(tech)
        ax1.set_xlabel("Leak flux (g/s)")
        ax1.set_ylabel("Leaks found per well")
        ax2.set_ylabel("Cumulative leakage fraction")
        ax1.set_ylim(0, 2)
        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        handles = handles1
        handles.extend(handles2)
        labels = labels1
        labels.extend(labels2)
        fig.legend(handles, labels, bbox_to_anchor=(0.9, 0.7))
        ind += 1
        display_save(display, save, file_out)
