# imports
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from simulation import Simulation
from tqdm import tqdm

# set plotting style
plt.style.use(['default', 'biophysics.mplstyle'])

# define the kymograph function
def kymo_plot(sim: Simulation, save_file: Path):
    # make a plot that uses cell coloring to analyze a 2D array
    fig, ax = plt.subplots()

    # get the raw grid and time vec
    time = sim.get_time_vec()
    grid = sim.get_raw_grid()

    # create custom cmap from color_dict
    color_dict = sim.get_color_dict()
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'custom-cmap',
        [color_dict[0], color_dict[1], color_dict[2], color_dict[3]]
    )

    # generate the grid plot
    # force a square ratio by calculating the aspect ratio
    aspect = grid.shape[1] / grid.shape[0]
    ax.imshow(grid, cmap=cmap, interpolation='none', aspect=aspect)

    # change axis limits
    ax.set_xlim(0, grid.shape[1])
    ax.set_ylim(grid.shape[0], 0)

    # rescale the y-axis from steps to time
    tick_interval = 100
    # get the nearest floored integer of length thats a multiple of the tick count
    tick_count = int(np.floor(time[-1] / tick_interval))
    # generate the labels
    labels = np.arange(0, tick_count+1) * tick_interval
    # calculate the labels position based off the max time
    label_positions = labels / time[-1] * grid.shape[0]
    # set the ticks and labels
    ax.set_yticks(ticks=label_positions, labels=labels)

    # rescale the x-axis from domains to position
    tick_interval = 5
    # get the nearest floored integer of length thats a multiple of the tick count
    tick_count = int(np.floor(sim.get_max_length() / tick_interval))
    # generate the labels
    labels = np.arange(0, tick_count+1) * tick_interval
    # calculate the labels position based off the max length
    label_positions = labels / sim.get_max_length() * grid.shape[1]
    # set the ticks and labels
    ax.set_xticks(ticks=label_positions, labels=labels)


    # add plot labels
    ax.set_xlabel(r"Position along MT $\left[\qty{}{\micro\meter}\right]$")
    ax.set_ylabel(r"Time $\left[\qty{}{\second}\right]$")
    ax.set_title(f"Simulation image of ``{sim.get_name()}\'\'")

    # add a legend
    handles = [
        mpl.patches.Patch(facecolor=color_dict[sim.get_grid_type('tau')], label='Tau', edgecolor='black'),
        mpl.patches.Patch(facecolor=color_dict[sim.get_grid_type('map6')], label='MAP6', edgecolor='black'),
        mpl.patches.Patch(facecolor=color_dict[sim.get_grid_type('empty')], label='Empty', edgecolor='black'),
        mpl.patches.Patch(facecolor=color_dict[sim.get_grid_type('notexist')], label='Not exist', edgecolor='black'),
    ]
    ax.legend(handles=handles, loc='upper right')

    # save the figure
    save_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_file)


# given a directory, 'paper_data2', find all .mat files in the directory
# also including child directories
# return the parent directory and the file name
def find_files(directory: Path):
    # get the list of all files in the directory
    files = directory.rglob('*.mat')
    # iterate through the files
    for file in files:
        # get the parent directory
        parent = file.parent
        # get the file name
        name = file.name
        # yield the parent and name
        yield parent, name


# iterate over the directory 'paper_data2' and generate kymographs of each .mat file
# save each kymograph to paper_figures/parent/"sim_kymo_{name}.pdf"
def main():
    # get the parent directory
    parent = Path('paper_data2')
    # iterate over the files
    for parent, name in tqdm(list(find_files(parent)), desc='Generating kymographs'):
        # create the simulation object
        sim = Simulation(mat_file=parent / name, skip_check=True)
        sim.set_color_dict({
            sim.get_grid_type('tau'): 'tab:blue',
            sim.get_grid_type('map6'): 'tab:orange',
            sim.get_grid_type('empty'): 'white',
            sim.get_grid_type('notexist'): 'grey',
        })
        # generate the kymograph
        kymo_plot(sim, parent / f"sim_kymo_{name}.pdf")


if __name__ == '__main__':
    main()
