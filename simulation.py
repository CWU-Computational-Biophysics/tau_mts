# simulation.py

# A class representing simulation data from a .mat file.
# TODO: Add docstrings

# imports
from os import PathLike
from pathlib import Path

import numpy as np
import numpy.typing as npt
import scipy.io as sio


class Simulation:

    def __init__(self, mat_file: PathLike):
        # load the data file with scipy
        data = sio.loadmat(
            mat_file,
            squeeze_me=True,
            simplify_cells=True
        )

        # save the three data sets
        self.params = data['params']
        self.grids = data['grids']
        self.growths = data['growths']

        # save the vector arrays
        self.length_vec = data['mt_length']
        self.state_vec = data['mt_state']

        # save the MT grid
        self.grid = data['mt_grid']

        # calculate and save the number of simulation steps
        self.nsteps = len(self.length_vec)

        # make and verify trimmed grids
        self.trim_grid = self._gen_trimmed_grids()


    def _gen_trimmed_grids(self) -> npt.ArrayLike:
        # make placeholder list
        trim_grid = []

        # iterate over each grid
        for si in np.arange(self.nsteps):
            # get the grid at the current step
            grid = self.grid[si]

            # find occurences of NOTEXIST
            notexists = np.where(grid == self.grids['NOTEXIST'])[0]

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


    def _gen_sequences(self) -> npt.ArrayLike:
        pass


    def get_nsteps(self) -> int:
        return self.nsteps


    def get_param(self, param: str) -> int | float | bool:
        # ensure parameter is valid
        if param not in self.params:
            raise ValueError(f"Invalid parameter: {param}")

        return self.params[param]


    def get_grid_at(self, step: int) -> npt.ArrayLike:
        # ensure step is valid
        if step < 0 or step >= self.nsteps:
            raise ValueError(f"Invalid step: {step}")

        return self.grid[step]

    def get_trimmed_grid_at(self, step: int) -> npt.ArrayLike:
        # ensure step is valid
        if step < 0 or step >= self.nsteps:
            raise ValueError(f"Invalid step: {step}")

        return self.trim_grid[step]
