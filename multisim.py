# multisim.py

# A class designed to collect statistics on many simulations.
# TODO: Add docstrings

# imports
from os import PathLike
from pathlib import Path
from typing import Iterable

from rich import print

from simulation import Simulation


class MultiSim:

    def __init__(self, mat_files: list[PathLike] = None, file_dir: PathLike = None, color_dict: dict = None, marker_dict: dict = None, skip_check: bool = True, sort_by: str = None):
        # check that only one of mat_files and file_dir are None
        if mat_files is None and file_dir is None:
            raise ValueError("Either mat_files or file_dir must be provided, not both")

        # save the data to be passed to Simulation objects
        self.color_dict = color_dict
        self.marker_dict = marker_dict
        self.skip_check = skip_check

        # define a list to store the Simulation objects
        self.sim_list = []

        # if mat_files is provided, load the files
        if mat_files:
            self.load_files(mat_files)

        # if file_dir is provided, load the files
        if file_dir:
            self.load_dir(file_dir)

        # if a sorting parameter is provided, call sort by
        if sort_by:
            self.sort_by(sort_by)


    def load_files(self, mat_files: list[PathLike], overwrite: bool = False):
        # if sim_list is not empty AND overwrite = False, error
        if self.sim_list and not overwrite:
            print("[red]Error:[/red] sim_list is not empty")
            return

        # check for overwrite
        if overwrite:
            self.sim_list = []

        # loop over each file and attempt to create a Simulation object
        for file in mat_files:
            # convert the file to a Path
            file_path = Path(file)

            # check that the file exists with the .mat extension
            if not file_path.exists():
                print(f"[red]Error:[/red] {file} does not exist")
                continue
            elif file_path.suffix != ".mat":
                print(f"[red]Error:[/red] {file} is not a .mat file")
                continue

            # attempt to create a Simulation object
            sim = Simulation(
                mat_file=file_path,
                color_dict=self.color_dict,
                marker_dict=self.marker_dict,
                skip_check=self.skip_check
            )

            # append to the sim_list
            self.sim_list.append(sim)


    def load_dir(self, file_dir: PathLike, overwrite: bool = False):
        # if sim_list is not empty AND overwite = False, error
        if self.sim_list and not overwrite:
            print("[red]Error:[/red] sim_list is not empty")
            return

        # get all the files in the file_dir with the extension .mat
        file_dir = Path(file_dir)
        mat_files = list(file_dir.glob("*.mat"))

        # if mat_files is empty, error
        if not mat_files:
            print(f"[red]Error:[/red] no .mat files found in {file_dir}")
            return

        # pass the files to load_files
        self.load_files(mat_files=mat_files, overwrite=overwrite)


    def remove_mask(self, mask: list) -> None:
        # check that mask is the same length as sim_list
        if len(mask) != len(self.sim_list):
            print("[red]Error:[/red] mask is not the same length as sim_list")
            return

        # remove the Simulation objects from sim_list
        self.sim_list = [sim for sim, m in zip(self.sim_list, mask) if not m]


    def get_iter(self) -> Iterable[Simulation]:
        return iter(self.sim_list)


    def sort_by(self, param: str) -> None:
        # check that param is a valid parameter in each simulation
        if not all(param in sim.param_dict for sim in self.get_iter()):
            print("[red]Error:[/red] param is not a valid parameter for all simulations")
            return

        # sort the sim_list by the parameter
        self.sim_list.sort(key=lambda x: x.get_param(param))


    def set_order(self, order: list) -> None:
        # check that order is the same length as sim_list
        if len(order) != len(self.sim_list):
            print("[red]Error:[/red] order is not the same length as sim_list")
            return

        # sort the sim_list by the order
        self.sim_list = [sim for _, sim in sorted(zip(order, self.sim_list))]
