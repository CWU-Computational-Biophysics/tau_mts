# plotter.py
# generates all standard figures from a directory

# imports
import os
from os import PathLike
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps

from configurations import DEF_VAR_LIST, DEF_REM_COL_LIST, CSTYLE_FILE_STR
from grouped_bar_plot import grouped_bar_plot
from processing import load_mat_dir


def main(data_dir: PathLike):
    # configure the figure save directory
    FIGURE_DIR = Path("figures", data_dir.parts[-1])
    os.makedirs(FIGURE_DIR, exist_ok=True)

    # set plotting style
    plt.style.use(["default", CSTYLE_FILE_STR])

    # set the colormap
    DEFAULT_COLORMAP = list(colormaps["tab10"].colors)

    # get the data dictionary from the data directory
    data_dict = load_mat_dir(data_dir, DEF_VAR_LIST, DEF_REM_COL_LIST)

    # add the final MT length to the sim_df
    data_dict["sims"].keys()
    length_vals = {}
    for sim_name, sim_dict in data_dict["sims"].items():
        # extract the final length as the last entry in first col of mt_length
        length_vals[sim_name] = sim_dict["mt_length"][-1, 0]

    # add the dictionary to the df
    data_dict["df"]["final_length_units"] = data_dict["df"].index.map(length_vals)

    # add final length column
    data_dict["df"]["final_length"] = data_dict["df"]["final_length_units"]*data_dict["df"]["dx"]

    # create an averages df grouped by tm_ratio
    data_avg_df = data_dict["df"].groupby("tm_ratio").mean()


    # define the grouped data and group labels
    group_data = {"L": data_avg_df["final_length"]}
    group_labels = data_avg_df.index

    # define the figure
    fig, ax = plt.subplots()

    # plot the grouped bar plot
    fig, ax = grouped_bar_plot(
        fig, ax,
        group_data,
        group_labels,
        multiplier=1,
        colors=DEFAULT_COLORMAP)

    # add a margin via xmargin to the x axis
    ax.margins(x=0.1)

    # set the labels
    ax.set_xlabel(r"Ratio of tau binding rate to map6 binding rate")
    ax.set_ylabel(r"Final microtubule length [$\qty{}{\micro\meter}$]")
    ax.set_title(r"Effect of binding ratio on microtubule growth")

    # save the figure
    fig.savefig(FIGURE_DIR / r"growth_effect")


    # make a grouped bar plot of binding asymmetry
    # on the horizontal: ratio of tau binding rate to map6 binding rate
    # on the vertical:
    # 1. tau plus end density to stable domain density ratio (tauasymm)
    # 2. map6 plus end density to stable domain density ratio (mapasymm)

    # define the group data and group labels
    group_data = {"tau asymmetry": data_avg_df["tau_plus_end_asym"], "map6 asymmetry": data_avg_df["map6_plus_end_asym"]}
    group_labels = data_avg_df.index

    # define the figure
    fig, ax = plt.subplots()

    # plot the data
    fig, ax = grouped_bar_plot(
        fig, ax,
        group_data,
        group_labels,
        multiplier=1,
        colors=DEFAULT_COLORMAP)

    # add a margin via xmargin to the x axis
    ax.margins(x=0.1)

    # set the labels
    ax.set_xlabel(r"Ratio of tau binding rate to map6 binding rate")
    ax.set_ylabel(r"Plus end density to stable domain density ratio")
    ax.set_title(r"Final binding asymmetry")

    # set the legend
    ax.legend(loc="upper left")

    # save the figure
    fig.savefig(FIGURE_DIR / r"asymm_effect")


    # make a grouped bar plot
    # on the horizontal: ratio of tau binding rate to map6 binding rate
    # on the vertical:
    # 1. plus-end tau fraction (taufractip)
    # 2. plus-end map6 fraction (mapfractip)
    # 3. stable domain tau fraction (taufraclength)
    # 4. stable domain map6 fraction (mapfraclength)

    # define the group data and group labels (T0/M0)
    group_data = {
        "tau tip fraction": data_avg_df["tau_frac_tip"],
        "map6 tip fraction": data_avg_df["map6_frac_tip"],
        "tau length fraction": data_avg_df["tau_frac_length"],
        "map6 length fraction": data_avg_df["map6_frac_length"],
    }
    group_labels = data_avg_df.index

    # define the colors
    colors = ["orange", "blue", "grey", "green"]

    # define the figure
    fig, ax = plt.subplots()

    # plot the data
    fig, ax = grouped_bar_plot(
        fig, ax,
        group_data,
        group_labels,
        multiplier=-1/2,
        colors=DEFAULT_COLORMAP)

    # add a margin via xmargin to the x axis
    ax.margins(x=0.1)

    # set the labels
    ax.set_xlabel(r"Ratio of tau binding rate to map6 binding rate (map6 rate constant)")
    ax.set_ylabel(r"Fraction of binding sites bound with this MAP")
    ax.set_title(r"Tau and Map6 distribution on MT: dependence on tau on rate")

    # add a legend
    ax.legend()

    # save the figure
    fig.savefig(FIGURE_DIR / r"binding_distribution_grouped")


    # plot the length of the mt over time
    # define the figure
    fig, ax = plt.subplots()

    # define a dictionary to color similar tau_ratios
    tau_ratio_dict = {}

    # iterate over each simulation in sim_dict and get the time and length array
    sim_names = list(data_dict["sims"].keys())
    sim_names.sort(key=lambda x: data_dict["df"].loc[x]["tm_ratio"])
    for index, sim_name in enumerate(sim_names):
        # get the arrays
        mt_length_array = data_dict["sims"][sim_name]["mt_length"]
        length_list = mt_length_array[:, 0]

        # get the time information
        dt = data_dict["df"].loc[sim_name]["dt"]
        ttot = data_dict["df"].loc[sim_name]["ttot"]
        time_list = np.arange(0, ttot, dt)

        # get the tau binding rate T0 ratio
        tau_ratio = data_dict["df"].loc[sim_name]["tm_ratio"]

        # define the label
        label = r"$\tau_{\text{on}}/\text{map6}_\text{on}=" + f"{tau_ratio}" + r"\ $"

        # get a unique color for each tau ratio
        # if the color has been used, make the label None
        if tau_ratio not in tau_ratio_dict.keys():
            tau_ratio_dict[tau_ratio] = colors[index]
        else:
            label = None

        # plot the length over time for this mt
        ax.plot(
            time_list,
            length_list,
            label=label,
            color=tau_ratio_dict[tau_ratio]
        )

    # set the axis limits
    # the x limit should be the final time plus the second time
    # this has the effect of showing an "even" value
    ax.set_xlim(0, time_list[-1] + time_list[1])

    # set the labels
    ax.set_xlabel(r"Time $\left[\qty{}{\second}\right]$")
    ax.set_ylabel(r"Microtubule length $\left[\qty{}{\micro\meter}\right]$")
    ax.set_title(r"Microtubule growth over time")

    # add a legend
    ax.legend()

    # save the figure
    fig.savefig(FIGURE_DIR / r"mt_growth_over_time")
