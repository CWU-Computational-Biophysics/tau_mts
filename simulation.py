# simulation.py

# A class representing simulation data from a .mat file.
# TODO: Add docstrings

# imports
from os import PathLike
from pathlib import Path
from typing import Callable, Iterator

import h5py
import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.io as sio

import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from matplotlib.animation import FuncAnimation

from rich import print
from tqdm import tqdm

import h5py
import numpy as np


def loadmat_v7_3(filepath, squeeze_me=True, simplify_cells=True):
    """
    Load MATLAB v7.3 .mat files using h5py with automatic transposition and
    conversion of 1x1 arrays to scalars.

    Parameters:
        filepath (str): Path to the .mat file.
        squeeze_me (bool): If True, squeeze unit dimensions from arrays.
        simplify_cells (bool): If True, convert MATLAB cell arrays to lists.

    Returns:
        dict: Dictionary containing MATLAB variables.
    """
    def mat_to_dict(mat_obj):
        """Recursively convert MATLAB objects to Python dictionaries/lists."""
        if isinstance(mat_obj, h5py.Dataset):
            # Convert datasets to numpy arrays
            data = mat_obj[()]
            if data.ndim > 1:  # Transpose multi-dimensional arrays
                data = data.T
            if squeeze_me:
                data = np.squeeze(data)
            if data.shape == ():  # Convert single-element arrays to scalars
                data = data.item()
            return data
        elif isinstance(mat_obj, h5py.Group):
            # Convert groups to dictionaries
            return {key: mat_to_dict(mat_obj[key]) for key in mat_obj.keys()}
        else:
            raise TypeError(f"Unsupported MATLAB object type: {type(mat_obj)}")

    def simplify(value):
        """Simplify MATLAB cells to Python lists."""
        if isinstance(value, np.ndarray) and value.dtype.kind == 'O':
            return [simplify(item) for item in value]
        elif isinstance(value, dict):
            return {k: simplify(v) for k, v in value.items()}
        else:
            return value

    with h5py.File(filepath, 'r') as mat_file:
        data = {key: mat_to_dict(mat_file[key]) for key in mat_file.keys()}
        if simplify_cells:
            data = simplify(data)
        return data


class Sequence:

    def __init__(self, start_index: int, end_index: int, type: int, grid_dict: dict):
        # check that protein type is in the values of grid_dict
        if type not in grid_dict.values():
            raise ValueError(f"Invalid protein type: {type}")

        # save the indices
        self.start_index = start_index
        self.end_index = end_index
        self.length = self.end_index - self.start_index

        # ensure length is positive
        if self.length <= 0:
            raise ValueError(f"Invalid start or end index: {start_index}, {end_index}")

        # save the protein type and name
        self.type = type
        self.protein = list(grid_dict.keys())[list(grid_dict.values()).index(self.type)]


    def __str__(self):
        return f"Sequence of {self.length} {self.protein} starting at {self.start_index}"


    def __repr__(self):
        return self.__str__()


    def get_start_index(self) -> int:
        return self.start_index


    def get_end_index(self) -> int:
        return self.end_index


    def get_length(self) -> int:
        return self.length


    def get_type(self) -> int:
        return self.type


class Simulation:

    def __init__(self, mat_file: PathLike, color_dict: dict = None, marker_dict: dict = None, skip_check: bool = False):
        # convert to Path and check that file exists
        mat_file = Path(mat_file)
        if not mat_file.exists():
            raise FileNotFoundError(f"File not found: {mat_file}")

        # set the simulation name based on the file path
        self.sim_name = mat_file.stem
        self.sim_file = mat_file

        # load the data file with scipy
        # data = sio.loadmat(
        #     mat_file,
        #     squeeze_me=True,
        #     simplify_cells=True
        # )

        # try to load the data file with h5py
        # fallback to scipy
        try:
            data = loadmat_v7_3(
                mat_file,
                squeeze_me=True,
                simplify_cells=True
            )
        except OSError:
            # print(f"[yellow]Warning: Could not load '{mat_file}' with h5py, falling back to scipy[/yellow]")
            try:
                data = sio.loadmat(
                    mat_file,
                    squeeze_me=True,
                    simplify_cells=True
                )
            except Exception as e:
                print(f"[red]Error: Could not load '{mat_file}': {e}[/red]")
                raise e


        # save the three data sets with lowercase keys for enums
        self.param_dict = data['params']
        self.grid_dict = {k.lower():v for k,v in data['grids'].items()}
        self.growth_dict = {k.lower():v for k,v in data['growths'].items()}

        # save the vector arrays
        self.length_vec = np.array(data['mt_length'])
        self.state_vec = np.array(data['mt_state'])

        # save the MT grid
        self.raw_grid = np.array(data['mt_grid'])

        # calculate and save the number of simulation steps and time vector
        self.nsteps = len(self.length_vec)
        self.time_vec = np.arange(self.get_nsteps()) * self.get_param('dt')

        # make and verify trimmed grids
        self.trim_grid = self._gen_trimmed_grids()

        # check if the color dict has the same keys as grid_dict
        # if there are missing keys, add them with value None
        self.color_dict = None
        if color_dict is not None:
            self.set_color_dict(color_dict)
        else:
            # create a default color dict
            # print(self.grid_dict)
            color_dict = {self.get_grid_type(k):None for k in self.grid_dict}
            color_dict[self.get_grid_type('tau')] = 'tab:blue'
            color_dict[self.get_grid_type('map6')] = 'tab:orange'
            self.color_dict = color_dict

        # check if the marker dict has the same keys as grid_dict
        # if there are missing keys, add them with value None
        self.marker_dict = None
        if marker_dict is not None:
            self.set_marker_dict(marker_dict)
        else:
            # create a default marker dict
            marker_dict = {self.get_grid_type(k):None for k in self.grid_dict}
            marker_dict[self.get_grid_type('tau')] = 'o'
            marker_dict[self.get_grid_type('map6')] = 's'
            self.marker_dict = marker_dict

        # run sim self check
        if not skip_check:
            print(f"[green]Running self check on '{self.sim_name}'[/green]")
            result = self._sim_self_check()
            if not result:
                print(f"[red]Warning: Simulation '{self.sim_name}' failed self-check[/red]")

        # calculate final array state parameters
        self.tau_frac_tip = 0
        self.tau_frac_length = 0
        self.tau_plus_end_asym = 0
        self.map6_frac_tip = 0
        self.map6_frac_length = 0
        self.map6_plus_end_asym = 0
        self._calc_final_array_states()


    def _calc_final_array_states(self, final_n: int = 10):
        # final_n is the number of grid points to be considered as the tip
        # get the tip array and length array from the last grid
        trimmed_grid = self.get_trimmed_grid_at(self.get_nsteps()-1)

        # check that the grid contains more than final_n points
        if len(trimmed_grid) < final_n:
            print(f"[red]Warning: Final grid at step {self.get_nsteps()-1} has less than {final_n} points. No final array state parameters will be calculated.[/red]")

        # cut the grid into two at end-final_n
        tip_grid = trimmed_grid[-final_n:]
        length_grid = trimmed_grid[:-final_n]

        # sum the protein types at each segment
        tip_tau = np.sum(tip_grid == self.get_grid_type('tau'))
        tip_map6 = np.sum(tip_grid == self.get_grid_type('map6'))
        length_tau = np.sum(length_grid == self.get_grid_type('tau'))
        length_map6 = np.sum(length_grid == self.get_grid_type('map6'))

        # calculate the final binding fraction asymmetry
        self.tau_frac_tip = tip_tau / len(tip_grid)
        self.tau_frac_length = length_tau / len(length_grid)
        self.map6_frac_tip = tip_map6 / len(tip_grid)
        self.map6_frac_length = length_map6 / len(length_grid)

        # handle divide by 0
        self.tau_plus_end_asym = 0
        if self.tau_frac_length != 0:
            self.tau_plus_end_asym = self.tau_frac_tip / self.tau_frac_length

        self.map6_plus_end_asym = 0
        if self.map6_frac_length != 0:
            self.map6_plus_end_asym = self.map6_frac_tip / self.map6_frac_length


    def get_tau_frac_tip(self) -> float:
        return self.tau_frac_tip

    def get_tau_frac_length(self) -> float:
        return self.tau_frac_length

    def get_tau_plus_end_asym(self) -> float:
        return self.tau_plus_end_asym

    def get_map6_frac_tip(self) -> float:
        return self.map6_frac_tip

    def get_map6_frac_length(self) -> float:
        return self.map6_frac_length

    def get_map6_plus_end_asym(self) -> float:
        return self.map6_plus_end_asym


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
        # print(f"Start protein: {start_protein}, index: {start_index}")

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


    def _sim_self_check(self) -> bool:
        # validate that mt grid count only increases when mt state is growing
        # likewise, mt grid count only decrases when mt state is shrinking
        # likewise, mt grid count does not change when mt state is stable

        # iterate over all steps
        has_passed = True
        for si in np.arange(self.get_nsteps() - 1):
            # retrieve the current length, grid count, and state
            length = self.get_length_at(si)
            grids = self.get_length_units_at(si)
            state = self.get_growth_state_at(si)

            # retrieve the next length and state
            next_length = self.get_length_at(si + 1)
            next_grids = self.get_length_units_at(si + 1)
            # next_state = self.get_growth_state_at(si + 1)

            # calculate delta length and growth rate
            delta_length = next_length - length
            delta_grid = next_grids - grids
            growth_rate = delta_length / self.get_param('dt')

            # error if the length decreases when the state is growing
            if state == self.get_growth_type('growing') and delta_length < 0:
                has_passed = False
                print(f"[yellow]Warning: Length decreased ({delta_length}) under growth at step {si}[/yellow]")

            # error if the length increases when the state is shrinking
            if state == self.get_growth_type('shrinking') and delta_length > 0:
                has_passed = False
                print(f"[yellow]Warning: Length increased ({delta_length}) under shrinkage at step {si}[/yellow]")

            # error if the length changes when the state is stable
            if state == self.get_growth_type('stable') and delta_length != 0:
                has_passed = False
                print(f"[yellow]Warning: Length changed ({delta_length}) under stability at step {si}[/yellow]")

            # error if the growth rate is greater than vp with slight margin for floating point error
            if abs(growth_rate) > self.get_param('vp') + 1e-6:
                has_passed = False
                print(f"[yellow]Warning: Length changed ({delta_length}) faster ({growth_rate}) than vp at step {si}[/yellow]")

            # error if the grid count changes by more than 1
            if abs(delta_grid) > 1:
                has_passed = False
                print(f"[yellow]Warning: Grid count changed by {delta_grid} at step {si}[/yellow]")


            # check for inconsistent NOTEXIST sequence
            # find occurences of NOTEXIST
            notexists = np.where(self.get_grid_at(si) == self.get_grid_type('notexist'))[0]

            # ensure that after the first notexist index
            # each subsequent index is one greater than the last
            if not np.all(np.diff(notexists) == 1):
                has_passed = False
                print(f"[yellow]Warning: Inconsistent NOTEXIST chain at step {si}[/yellow]")

        # return the result
        return has_passed


    def set_color_dict(self, color_dict: dict) -> None:
        # check that color_dict has the same keys as grid_dict
        # if there are missing keys, add them with value None
        for key in list(self.grid_dict.values()):
            if key not in color_dict:
                color_dict[key] = None

        # report a warning for extra keys in color dict
        for key in color_dict:
            if key not in list(self.grid_dict.values()):
                print(f"Warning: Extra key in color dict: {key}")

        # store the color dict
        self.color_dict = color_dict


    def set_marker_dict(self, marker_dict: dict) -> None:
        # check that marker_dict has the same keys as grid_dict
        # if there are missing keys, add them with value None
        for key in self.grid_dict:
            if key not in marker_dict:
                marker_dict[key] = None

        # report a warning for extra keys in marker dict
        for key in marker_dict:
            if key not in self.grid_dict:
                print(f"Warning: Extra key in marker dict: {key}")

        # store the marker dict
        self.marker_dict = marker_dict


    def get_nsteps(self) -> int:
        return self.nsteps


    def get_param(self, param: str) -> int | float | bool:
        # ensure parameter is valid
        if param not in self.param_dict:
            raise ValueError(f"Invalid parameter: {param}")

        return self.param_dict[param]


    def get_growth_state_at(self, step: int) -> int:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        return self.state_vec[step]


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


    def get_length_vec(self) -> npt.ArrayLike:
        return self.length_vec


    def get_time_vec(self) -> npt.ArrayLike:
        return self.time_vec


    def get_length_diff_at(self, step: int) -> float:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # get the true length at the step
        true_len = self.get_length_at(step)

        # get the length in domains
        domain_len = self.get_length_units_at(step) * self.get_param('dx')

        # report the difference
        return true_len - domain_len


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
        df['domain_end_pos'] = (df['index'] + 1) * dx
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


    def get_longest_seq(self, step: int = None) -> float:
        # if step is provided, then ensure it is valid
        if step is not None and not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # if step is provided, return the longest sequence at that step
        if step is not None:
            return max([seq.get_length() for seq in self._seq_iter(step)])

        # otherwise get the longest over all steps
        seq_list = [
            max([seq.get_length() for seq in self._seq_iter(si)])
            for si in np.arange(self.get_nsteps())
        ]
        return max(seq_list)


    def calc_plot_y_lims_at(self, step: int, override: float = None) -> float:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # get the largest value on the y-axis
        # this is the largest length of a sequence at this step
        longest_seq = self.get_longest_seq(step=step)

        # calculate absolute ymax
        abs_ymax = longest_seq
        if override:
            abs_ymax = override

        # force largest value to be the next highest multiple of 5
        ymax = np.ceil(abs_ymax / 5) * 5

        return ymax


    def add_plot_elements(self, step: int, ax: Axes, max_len: bool = True) -> None:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # set the plot x-limits by length
        # either by true max or current max
        if max_len:
            ax.set_xlim(0, self.get_max_length())
        else:
            ax.set_xlim(0, self.get_length_at(step))

        # set the labels
        ax.set_xlabel(r"Position along microtubule $\left[\qty{}{\micro\meter}\right]$")
        ax.set_title(r"Protein cluster distribution along microtubule")


    def plot_sequence_at(self, step: int, ax: Axes, one_side: bool = False, max_len: bool = True, plot_points: bool = False, point_domain_ticks: bool = False, ymax_override: float = None) -> None:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # add default plot elements
        self.add_plot_elements(step, ax, max_len)

        # set the y-limits by calculation
        ymax = self.calc_plot_y_lims_at(step, override=ymax_override)
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
                fmt='o',
                markersize=0,
                color=color,
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

        # draw a horizontal black line at y=0 from x=0 to the current length of the MT
        mt_len = self.get_length_at(step)
        ax.plot(
            [0, mt_len],
            [0, 0],
            color='black',
        )

        # add y-axis label
        ax.set_ylabel(r"Number of proteins in uninterupted sequence")

        # add a legend with a label for map6 and tau using artists
        # the label should be a colored rectangle
        tau_color = self.color_dict[self.get_grid_type('tau')]
        map6_color = self.color_dict[self.get_grid_type('map6')]
        tau_label = plt.Rectangle((0, 0), 1, 1, fc=tau_color, edgecolor="black")
        map6_label = plt.Rectangle((0, 0), 1, 1, fc=map6_color, edgecolor="black")

        # add legend to the upper-left
        ax.legend(
            [tau_label, map6_label],
            [r"Tau", r"MAP6"],
            loc='upper left',
        )

        # possibly call plot_proteins_at
        if plot_points:
            self.plot_proteins_at(
                step=step,
                ax=ax,
                remove_yaxis=False,
                domain_ticks=point_domain_ticks,
                max_len=max_len)


    def plot_proteins_at(self, step: int, ax: Axes, remove_yaxis: bool = False, domain_ticks: bool = False, max_len: bool = True) -> None:
        # ensure step is valid
        if not self._valid_step(step):
            raise ValueError(f"Invalid step: {step}")

        # add default plot elements
        self.add_plot_elements(step, ax, max_len)

        # check for remove_yaxis
        if remove_yaxis:
            ax.yaxis.set_visible(False)

        # define the binding tick height as 2% of the y_lim
        ymax = self.calc_plot_y_lims_at(step)
        binding_tick_height = 0.02 * ymax

        # iterate over proteins
        for protein in self.get_protein_plot_points_at(step).itertuples(index=False):
            # if domain_ticks, add vertical lines at domain_end_pos
            if domain_ticks:
                ax.vlines(
                    x=protein.domain_end_pos,
                    ymin=-binding_tick_height,
                    ymax=binding_tick_height,
                    color='black',
                )

            # skip 'empty' type proteins
            if protein.type == self.get_grid_type('empty'):
                continue

            # error on 'notexist' type proteins
            if protein.type == self.get_grid_type('notexist'):
                raise ValueError(f"Error at step {step}: cannot plot NOTEXIST proteins")

            # check that the protein type is in the color dict
            if protein.type not in self.color_dict:
                raise ValueError(f"Error at step {step}: color not defined for protein: {protein.type}")

            # get the color for this protein
            color = self.color_dict[protein.type]

            # get the marker for this protein
            marker = self.marker_dict[protein.type]

            # plot the protein as a circle with a black edge
            ax.scatter(
                x=protein.protein_draw_pos,
                y=0,
                s=100,
                c=color,
                edgecolor='black',
                zorder=2,
                marker=marker,
            )

        # add a legend with a label for map6 and tau using artists
        # the label should be a colored rectangle
        # the label should respect the marker choice
        tau_color = self.color_dict[self.get_grid_type('tau')]
        map6_color = self.color_dict[self.get_grid_type('map6')]
        tau_marker = self.marker_dict[self.get_grid_type('tau')]
        map6_marker = self.marker_dict[self.get_grid_type('map6')]
        tau_label = plt.Line2D(
            [0], [0],
            marker=tau_marker,
            color='w',
            markerfacecolor=tau_color,
            markersize=10,
            markeredgecolor='black',
        )
        map6_label = plt.Line2D(
            [0], [0],
            marker=map6_marker,
            color='w',
            markerfacecolor=map6_color,
            markersize=10,
            markeredgecolor='black'
        )
        # tau_label = plt.Rectangle((0, 0), 1, 1, fc=tau_color, edgecolor="black")
        # map6_label = plt.Rectangle((0, 0), 1, 1, fc=map6_color, edgecolor="black")

        # add legend to the upper-left
        # append, do not erase
        ax.legend(
            [tau_label, map6_label],
            [r"Tau", r"MAP6"],
            loc='upper left',
        )


    def animate_sequences(self, fig: Figure, anim_func: Callable, save_file: PathLike, progress_bar: bool = True) -> None:
        # ensure the save_file path exists
        save_file = Path(save_file)
        save_file.parent.mkdir(parents=True, exist_ok=True)

        # run the animation generator
        with tqdm(
            total=self.get_nsteps(),
            desc=f"Generating Animation: '{save_file.name}'",
            disable=(not progress_bar)
            ) as pbar:

            # pbar update function
            def update(*args):
                pbar.update(1)

            # generate the animation
            frames = self.get_nsteps()
            animation = FuncAnimation(fig, anim_func, frames=frames, repeat=False)

            # iterate over the animation function and save
            animation.save(save_file, writer='ffmpeg', progress_callback=update, dpi=100, fps=30)


    def get_name(self) -> str:
        return self.sim_name


    def get_color_dict(self) -> dict:
        return self.color_dict


    def get_raw_grid(self) -> npt.ArrayLike:
        return self.raw_grid
