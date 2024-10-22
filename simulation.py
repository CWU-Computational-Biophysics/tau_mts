# simulation.py

# A class representing simulation data from a .mat file.
# TODO: Add docstrings

# imports
from os import PathLike
from pathlib import Path
from typing import Iterator

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.io as sio

import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib import colormaps
from matplotlib.ticker import MaxNLocator
from matplotlib.animation import FuncAnimation

from sequence import Sequence


class Simulation:

    def __init__(self, mat_file: PathLike, color_dict: dict = None):
        # convert to Path and check that file exists
        mat_file = Path(mat_file)
        if not mat_file.exists():
            raise FileNotFoundError(f"File not found: {mat_file}")

        # load the data file with scipy
        data = sio.loadmat(
            mat_file,
            squeeze_me=True,
            simplify_cells=True
        )

        # save the three data sets with lowercase keys for enums
        self.param_dict = data['params']
        self.grid_dict = {k.lower():v for k,v in data['grids'].items()}
        self.growth_dict = {k.lower():v for k,v in data['growths'].items()}

        # save the vector arrays
        self.length_vec = np.array(data['mt_length'])
        self.state_vec = np.array(data['mt_state'])

        # save the MT grid
        self.raw_grid = np.array(data['mt_grid'])

        # calculate and save the number of simulation steps
        self.nsteps = len(self.length_vec)

        # make and verify trimmed grids
        self.trim_grid = self._gen_trimmed_grids()

        # check if the color dict has the same keys as grid_dict
        # if there are missing keys, add them with value None
        if color_dict is not None:
            for key in self.grid_dict:
                if key not in color_dict:
                    color_dict[key] = None

        # store the color dict
        self.color_dict = color_dict


    def _gen_trimmed_grids(self) -> list[npt.ArrayLike]:
        # make placeholder list
        trim_grid = []

        # iterate over each grid
        for si in np.arange(self.nsteps):
            # get the grid at the current step
            grid = self.get_grid_at(si)

            # find occurences of NOTEXIST
            notexists = np.where(grid == self.get_grid_type('notexist'))[0]

            # if no notexist, then MT fully exists
            if notexists.size == 0:
                trim_grid.append(grid)
            else:
                # predict the total number of notexist indices
                notexist_predict = notexists[-1] - notexists[0] + 1

                # compare to the number of notexist indices
                if notexist_predict != len(notexists):
                    # print a warning
                    print(f"Warning: Inconsistent NOTEXIST chain at step {si}")
                else:
                    # append a trimmed grid
                    first_notexist = notexists[0]
                    trim_grid.append(grid[:first_notexist])

        return trim_grid


    def _valid_step(self, step: int) -> bool:
        return step >= 0 and step < self.get_nsteps()


    def _gen_sequence_at(self, step: int):
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # get the trimmed grid at the current step
        grid = self.get_trimmed_grid_at(step)

        # verify that all entries of grid
        # are present as values in self.grid_dict
        if not np.all(np.isin(grid, list(self.grid_dict.values()))):
            raise ValueError(f"grid at step {step} contains unkown values")

        # define the list to hold sequence objects
        sequence_list = []

        # initialize the sequence finder with the first index and protein
        start_index = 0
        start_protein = grid[0]
        print(f"Start protein: {start_protein}, index: {start_index}")

        # iterate over the proteins in grid
        for i, protein in enumerate(grid):
            # if the protein changes, then add the sequence
            if protein != start_protein:
                # append the sequence object to the list
                sequence = Sequence(
                    start_index,
                    i,
                    start_protein,
                    self.grid_dict,
                )

                sequence_list.append(sequence)

                # update the start index and protein
                start_index = i
                start_protein = protein

        # append the last sequence
        sequence = Sequence(
            start_index,
            len(grid),
            start_protein,
            self.grid_dict,
        )

        sequence_list.append(sequence)

        # return the sequence list
        return sequence_list


    def get_nsteps(self) -> int:
        return self.nsteps


    def get_param(self, param: str) -> int | float | bool:
        # ensure parameter is valid
        if param not in self.param_dict:
            raise ValueError(f"Invalid parameter: {param}")

        return self.param_dict[param]


    def get_grid_type(self, grid: str) -> int:
        # lowercase grid
        grid = grid.lower()

        # ensure grid is valid
        if grid not in self.grid_dict:
            raise ValueError(f"Invalid grid: {grid}")

        return self.grid_dict[grid]


    def get_growth_type(self, growth: str) -> int:
        # lowercase growth
        growth = growth.lower()

        # ensure growth is valid
        if growth not in self.growth_dict:
            raise ValueError(f"Invalid growth: {growth}")

        return self.growth_dict[growth]


    def get_grid_at(self, step: int) -> npt.ArrayLike:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        return self.raw_grid[step]


    def get_trimmed_grid_at(self, step: int) -> npt.ArrayLike:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        return self.trim_grid[step]


    def get_length_at(self, step: int) -> float:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        return self.length_vec[step]


    def get_length_units_at(self, step: int) -> int:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        return len(self.get_trimmed_grid_at(step))


    def get_max_length_units(self) -> int:
        return len(self.get_grid_at(0))


    def get_max_length(self) -> float:
        return np.max(self.length_vec)


    def get_state_at(self, step: int) -> float:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        return self.state_vec[step]


    def get_sequence_at(self, step: int) -> npt.ArrayLike:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        return self._gen_sequence_at(step)


    def _seq_iter(self, step: int) -> Iterator[Sequence]:
        return iter(self.get_sequence_at(step))


    def get_protein_plot_points_at(self, step: int) -> pd.DataFrame:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # iterate over each sequence object
        data = {
            'index': [],
            'type': [],
        }
        for seq in self._seq_iter(step):
            # get the sequence information
            start_index = seq.get_start_index()
            end_index = seq.get_end_index()
            type = seq.get_type()

            # append this information a number of times
            data['index'].extend(range(start_index, end_index))
            data['type'].extend([type] * seq.get_length())

        # convert to dataframe
        df = pd.DataFrame(data)

        # add a column 'pos' equal to 'index' * 'dx'
        dx = self.get_param('dx')
        df['domain_start_pos'] = df['index'] * dx
        df['protein_draw_pos'] = df['domain_start_pos'] + (dx/2)

        # check that 'index' col is equal to the index of the df
        if not np.all(df.index == df['index']):
            raise ValueError("Index column does not match dataframe index, possible missing proteins")

        # remove the index col
        df.drop(columns=['index'], inplace=True)

        return df


    def get_sequence_plot_points_at(self, step: int) -> pd.DataFrame:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # generate a data list from a sequence
        data_list = [
            {'start_index': seq.get_start_index(),
             'end_index': seq.get_end_index(),
             'length_units': seq.get_length(),
             'type': seq.get_type()} for seq in self._seq_iter(step)
        ]

        # convert to dataframe
        df = pd.DataFrame(data_list)

        # add columns for pos from index * dx
        df['start_pos'] = df['start_index'] * self.get_param('dx')
        df['end_pos'] = df['end_index'] * self.get_param('dx')
        df['length'] = df['length_units'] * self.get_param('dx')

        return df


    def calc_plot_y_lims_at(self, step: int, override: float = None) -> float:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # get the largest value on the y-axis
        # this is the largest length of a sequence at this step
        longest_seq = max([seq.get_length() for seq in self._seq_iter(step)])

        # calculate absolute ymax
        abs_ymax = longest_seq
        if override:
            abs_ymax = override

        # force largest value to be the next highest multiple of 5
        ymax = np.ceil(abs_ymax / 5) * 5

        return ymax


    def plot_sequence_at(self, ax: Axes, step: int, one_side: bool = False, max_len: bool = True) -> None:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # set the plot x-limits by length
        # either by true max or current max
        if max_len:
            ax.set_xlim(0, self.get_max_length())
        else:
            ax.set_xlim(0, self.get_length_at(step))

        # set the y-limits by calculation
        ymax = self.calc_plot_y_lims_at(step)
        if not one_side:
            # symmetric limits
            ax.set_ylim(ymin=-ymax, ymax=ymax)

            # relabel the y-axis in absolute value using a formatter
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{abs(x):.0f}"))
        else:
            # bars are all vertical
            ax.set_ylim(ymin=0, ymax=ymax)

        # iterate over the sequences at the current step
        for seq in self.get_sequence_plot_points_at(step).itertuples(index=False):
            # skip empty sequences
            if seq.type == self.get_grid_type('empty'):
                continue

            # error on notexist sequences
            if seq.type == self.get_grid_type('notexist'):
                raise ValueError("Cannot plot NOTEXIST sequences")

            # check that the protein type is in the color dict
            if seq.type not in self.color_dict:
                raise ValueError(f"Color not defined for protein: {seq.type}")
            color = self.color_dict[seq.type]

            # define a height inverter for map6 when not using one_side
            height_inverter = 1
            if (not one_side) and (seq.type == self.get_grid_type('map6')):
                height_inverter = -1

            # collect the points to plot
            start_x = seq.start_pos
            end_x = seq.end_pos
            end_y = seq.length_units * height_inverter

            # plot a horizontal error bar for the top of the sequence bar
            ax.errorbar(
                x=(start_x + end_x) / 2,
                y=end_y,
                xerr=(end_x - start_x) / 2,
                fmt='none',
                ecolor=color,
                capsize=5,
                capthick=2,
            )

            # shade the region below the bar with no outline
            ax.fill_between(
                x=[start_x, end_x],
                y1=0,
                y2=end_y,
                color=color,
                alpha=1,
                linewidth=0,
                zorder=0,
            )

        # draw a horizontal black line at y=0
        ax.axhline(0, color='black', linewidth=1)

        # set the labels
        ax.set_xlabel(r"Position along microtubule $\left[\qty{}{\micro\meter}\right]$")
        ax.set_ylabel(r"Number of proteins in uninterupted sequence")
        ax.set_title(r"Protein cluster distribution along microtubule")

        # add a legend with a label for map6 and tau using artists
        # the label should be a colored rectangle
        tau_color = self.color_dict[self.get_grid_type('tau')]
        map6_color = self.color_dict[self.get_grid_type('map6')]
        tau_label = plt.Rectangle((0, 0), 1, 1, fc=tau_color, edgecolor="black")
        map6_label = plt.Rectangle((0, 0), 1, 1, fc=map6_color, edgecolor="black")

        # add legend to the upper-left
        ax.legend(
            [tau_label, map6_label],
            [r"Tau Sequence", r"Map6 Sequence"],
            loc='upper left',
        )
