# distributions.py
# store functions for plotting distributions

# imports
from os import PathLike
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import colormaps
from matplotlib.ticker import MaxNLocator
from matplotlib.animation import FuncAnimation

from tqdm import tqdm


# define a function to extract distribution information
def extract_distribution(sim_name: str, time_step: int, data_dict: dict) -> tuple:
    """ Extracts the distribution information from a simulation at a given time step.

    Arguments:
        sim_name {str} -- sim_name of the simulation in the data_dict
        time_step {int} -- time step of the simulation
        data_dict {dict} -- dictionary of simulation data

    Returns:
        tuple -- tuple of time_list, length_value, length_units, mt_state, dx
    """

    # extract time information
    time_list = data_dict["sims"][sim_name]["time"]

    # extract spatial information
    dx = data_dict["df"].loc[sim_name]["dx"]
    length_units = data_dict["sims"][sim_name]["mt_length"][time_step][0]
    length_value = length_units * dx

    # get this row of the mt_grid
    mt_state = data_dict["sims"][sim_name]["mt_grid"][time_step]
    mt_state = mt_state[0:length_units]

    # calculate grid positions
    # these are the positions of the left edges and the right edge
    grid_left_edge = np.arange(0, length_units) * dx

    # join grid positions and mt state as columns in np array
    mt_state = np.column_stack((grid_left_edge, mt_state))

    return time_list, length_value, length_units, mt_state, dx


# define a function to extract sequences
def extract_sequences(sim_name: str, time_step: int, data_dict: dict) -> pd.DataFrame:
    """ Extracts the protein sequences from a simulation at a given time step.

    Arguments:
        sim_name {str} -- sim_name of the simulation in the data_dict
        time_step {int} -- time step of the simulation
        data_dict {dict} -- dictionary of simulation data

    Returns:
        pd.DataFrame -- DataFrame of the protein sequences

    """

    # extract the distribution information
    _, _, _, mt_state_array, _ = extract_distribution(sim_name, time_step, data_dict)

    # define the sequence dictionary
    sequence_dict = {
        "start_pos": [],
        "width": [],
        "end_pos": [],
        "center_pos": [],
        "num_proteins": [],
        "protein_type": []
    }

    # track the start position of the current sequence and protein type
    start_pos = mt_state_array[0][0]
    start_index = 0
    protein_type = mt_state_array[0][1]

    # iterate over the mt_state_array
    for i, (pos, protein) in enumerate(mt_state_array):
        # check if the protein type has changed or its the last entry
        if protein != protein_type or i == len(mt_state_array) - 1:
            # calculate the sequence width and center
            domain_width = pos - start_pos
            domain_center = start_pos + domain_width / 2

            # append the sequence to the sequence_dict
            sequence_dict["start_pos"].append(start_pos)
            sequence_dict["width"].append(domain_width)
            sequence_dict["end_pos"].append(pos)
            sequence_dict["center_pos"].append(domain_center)
            sequence_dict["num_proteins"].append(int(i - start_index))
            sequence_dict["protein_type"].append(int(protein_type))

            # update the start position and protein type
            start_pos = pos
            start_index = i
            protein_type = protein

    # convert the sequence_dict to a dataframe
    sequence_df = pd.DataFrame(sequence_dict)

    # return
    return sequence_df


# define a function to plot clusters
def plot_clusters(
    fig, ax,
    sim_name: str,
    data_dict: dict,
    time_step: int = -1,
    colors: list = None,
    invert_map6: bool = True,
    protein_points_size: int = 0,
    horizontal_line: bool = True,
    binding_ticks: bool = False,
    legend_loc: str = "upper left",
    yabs_max_override: int = None,
    ):
    """Plots the protein clusters along a microtubule

    Arguments:
        fig -- figure object
        ax -- axis object
        sim_name {str} -- sim_name of the simulation in the data_dict
        data_dict {dict} -- dictionary of simulation data

    Keyword Arguments:
        time_step {int} -- time step of the simulation (default: {-1})
        colors {list} -- list of colors for the protein types (default: {None})
        invert_map6 {bool} -- invert the map6 axis for readability (default: {True})
        protein_points_size {int} -- size of the protein points on the horizontal axis (default: {0})
        horizontal_line {bool} -- draw a horizontal line at y=0 (default: {True})
        binding_ticks {bool} -- show ticks for each binding site on the horizontal axis (default: {False})
        legend_loc {str} -- location of the legend (default: {"upper left"})
        yabs_max_override {int} -- override the y-axis maximum value (default: {None})
    """

    # extract the distribution information
    _, length_value, _, mt_state_array, grid_spacing = extract_distribution(sim_name, time_step, data_dict)

    # extract the sequence information
    # this was previously called the domain information
    sequence_df = extract_sequences(sim_name, time_step, data_dict)

    # get a list of colors
    if colors is None:
        colors = colormaps["tab10"].colors

    # define tau and map6 colors
    tau_color = colors[0]
    map6_color = colors[1]

    # use the sequence_df to plot the sequence information
    for (_, start_pos, width, end_pos, center_pos, num_proteins, protein_type) in sequence_df.itertuples():
        # skip protein_type 0
        if protein_type == 0: continue

        # use the protein type to determine the color
        if protein_type == 1: color = tau_color
        elif protein_type == 2: color = map6_color

        # get height_scale and check for map6 inversion
        height_scale = 1
        if protein_type == 2 and invert_map6:
            height_scale = -1

        # use a point with a horizontal error bar
        ax.errorbar(
            center_pos,
            num_proteins*height_scale,
            xerr=width / 2,
            fmt="o",
            color="black",
            markersize=0
        )

        # shade the region below the bar with no outline
        ax.fill_between(
            [start_pos, end_pos],
            0,
            num_proteins*height_scale,
            color=color,
            alpha=1,
            linewidth=0,
            zorder=0,
        )

    # plot each protein on the horizontal if protein points
    if protein_points_size != 0:
        for (pos, protein) in mt_state_array:
            # protein = 1 means tau
            # protein = 2 means map6
            # check that the protein value is not 0
            if protein == 0: continue

            # use the protein value to determine the color
            if protein == 1: color = tau_color
            elif protein == 2: color = map6_color

            # plot the protein as a dot with a black outline
            protein_pos = pos + grid_spacing/2
            ax.scatter(protein_pos, 0, color=color, s=protein_points_size, edgecolors="black", zorder=2)

    # plot a horizontal line at y=0
    # the line should be the length of the mt at the current time step
    # plot the line in the background
    if horizontal_line:
        # ax.axhline(0, color="black", zorder=1)
        ax.plot([0, length_value], [0, 0], color="black", zorder=1)

    # get the largest value on the y axis after plotting
    if yabs_max_override is not None:
        yabs_max = yabs_max_override
    else:
        yabs_max = abs(max(ax.get_ylim(), key=abs))

    # force largest value to be the next highest multiple of 5
    yabs_max = 5 * (yabs_max // 5 + 1)

    # set the y axis limits
    ax.set_ylim(ymin=-yabs_max, ymax=yabs_max)

    # relabel the y-axis in absolute value using a formatter
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{abs(x):.0f}"))

    # plot small vertical bars for each binding site
    # the vertical bars should be centered on y=0
    # they vertical bars should extend vertically just a small amount
    if binding_ticks:
        # define the tick height
        height = 2 * ax.get_ylim()[0] / 100

        # draw a tick at the left edge of each binding site
        for pos in mt_state_array[:, 0]:
            ax.vlines(pos, -height, height, color="black")

        # draw one more at the end (at pos + grid_spacing)
        ax.vlines(pos + grid_spacing, -height, height, color="black")

    # set the xlim to (0, mt length at current time)
    ax.set_xlim(0, length_value)

    # set the labels
    ax.set_xlabel(r"Position along microtubule $\left[\qty{}{\micro\meter}\right]$")
    ax.set_ylabel(r"Number of proteins in uninterupted sequence")
    ax.set_title(r"Protein cluster distribution along microtubule")

    # add a legend with a label for map6 and tau using artists
    # the label should be a colored rectangle
    tau_label = plt.Rectangle((0, 0), 1, 1, fc=tau_color, edgecolor="black")
    map6_label = plt.Rectangle((0, 0), 1, 1, fc=map6_color, edgecolor="black")

    # define artists for the protein dots
    # include the outline
    tau_dot = plt.Line2D(
        [0], [0],
        marker="o",
        color="w",
        markerfacecolor=tau_color,
        markersize=10,
        markeredgewidth=1,
        markeredgecolor="black")
    map6_dot = plt.Line2D(
        [0], [0],
        marker="o",
        color="w",
        markerfacecolor=map6_color,
        markersize=10,
        markeredgewidth=1,
        markeredgecolor="black")

    # add legend to the upper-left depening on protein print status
    if protein_points_size != 0:
        ax.legend(
            [tau_label, tau_dot, map6_label, map6_dot],
            ["Tau Cluster", "Tau", "MAP6 Cluster", "MAP6"],
            loc=legend_loc)
    else:
        ax.legend(
            [tau_label, map6_label],
            ["Tau Sequence", "MAP6 Sequence"],
            loc=legend_loc)


# define the animation generation function
def create_protein_animation(
    save_path: PathLike,
    sim_name: str,
    data_dict: dict,
    frame_rate: int = 60,
    anim_time: int = 10,
    overwrite: bool = False,
    protein_points_size: int = 0,
    binding_ticks: bool = False,
    override_yaxis: bool = True):
    """Creates an animation of the protein clusters along the microtubule

    Arguments:
        save_path {PathLike} -- path to save the animation
        sim_name {str} -- name of the simulation in the sim_dict
        sim_dict {dict} -- dictionary of simulation data

    Keyword Arguments:
        frame_rate {int} -- frame rate of the animation in frames per second (default: {60})
        anim_time {int} -- desired runtime of the animation in seconds (default: {10})
        overwrite {bool} -- whether to overwrite the save_path if it already exists (default: {False})
        protein_points_size {int} -- size of the protein points on the horizontal axis (default: {0})
        binding_ticks {bool} -- show ticks for each binding site on the horizontal axis (default: {False})
        override_yaxis {bool} -- override the y-axis maximum value (default: {True})
    """

    # validate the save path (check that it is a path and exists and is empty)
    # the save path could be a string or other pathlike but it will be recast
    if isinstance(save_path, str):
        save_path = Path(save_path)

    if not isinstance(save_path, Path):
        raise ValueError("save_path must be a Path object.")

    if save_path.exists() and not overwrite:
        raise FileExistsError(f"File '{save_path}' already exists. Set 'overwrite=True' to overwrite the file.")

    # check that the sim dict has sims as a key
    if "sims" not in data_dict.keys():
        raise ValueError("The sim_dict must have a 'sims' key.")

    # validate the sim name
    if sim_name not in data_dict["sims"].keys():
        raise ValueError(f"Simulation '{sim_name}' not found in sim_dict. Available simulations: {list(data_dict['sims'].keys())}.")

    # get the frame count
    frame_count = anim_time * frame_rate

    # get time step information
    max_time_step = data_dict["df"].loc[sim_name]["steps"]
    time_steps = np.linspace(0, max_time_step-1, frame_count, dtype=int)
    time_list = data_dict["sims"][sim_name]["time"]

    # generate an animation from the frames
    # get the max mt length for the simulation
    max_mt_length = data_dict["df"].loc[sim_name]["final_length"]

    # get the largest cluster size
    # iterate over each time step and extract sequence information
    # then save the largest num_proteins
    max_protein_count = 0
    if override_yaxis:
        for time_step in time_steps:
            sequence_df = extract_sequences(sim_name, time_step, data_dict)
            max_protein_count = max(max_protein_count, sequence_df["num_proteins"].max())

    # if the max_protein_count is 0 then set it to none to let the drawing function handle it
    if max_protein_count == 0:
        max_protein_count = None

    # define the figure
    fig, ax = plt.subplots(figsize=(16, 9), layout="constrained")
    # fig.tight_layout()

    # define the animation
    def animate(i):
        # plot the clusters
        ax.clear()
        plot_clusters(
            fig, ax,
            sim_name=sim_name,
            data_dict=data_dict,
            time_step=time_steps[i],
            protein_points_size=protein_points_size,
            binding_ticks=binding_ticks,
            yabs_max_override=max_protein_count)

        # set the x-axis limits to the max mt length
        ax.set_xlim(0, max_mt_length)

        # add the frame number and time to the top center of the plot
        frame_time = time_list[time_steps[i]]
        step_str = f"Step: {time_steps[i]}/{max_time_step-1}"
        time_str = r"Time: \qty{" + f"{frame_time:.2f}" + r"}{\second}"
        ax.text(0.5, 0.05,
                time_str + r"\(\qquad\qquad\)" + step_str,
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes)

    # calculate the time interval from desired runtime and
    # frame count [second / frame]
    interval = anim_time / frame_count
    # [millisecond / frame]
    scaled_interval = interval * 1000

    # create the animation
    with tqdm(total=frame_count, desc="Generating Frames") as pbar:
        def update(*args):
            pbar.update(1)

        # generate the animation function
        animation = FuncAnimation(
            fig,
            animate,
            frames=frame_count,
            interval=scaled_interval,
            repeat=False)

        # iterate over the animation function and save
        pbar.set_description(f"Generating Animation: '{save_path.name}'")
        animation.save(save_path, writer="ffmpeg", progress_callback=update, dpi=100)

    # close the figure
    plt.close(fig)
