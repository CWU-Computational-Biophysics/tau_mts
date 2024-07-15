# processing.py
# data processing for 'tau_map6.m' simulation '.mat' files

# import packages
import os
import warnings
from os import PathLike
from pathlib import Path

import pandas as pd
import numpy as np
import scipy.io as sio

from rich import print

from configurations import DEF_VAR_LIST


# define functions
# define a function to extract data from a .mat file
def extract_mat_data(mat_file: PathLike, var_list: list = None) -> dict:
    """
    Extracts variables and arrays from a .mat file.

    Parameters:
    mat_file (Path): The path to the .mat file.

    Keyword Arguments:
    var_list (list): A list of variables to extract from the .mat file. (default: None, extract all variables)

    Returns:
    dict: A dictionary of variables.
    """

    # check that the mat_file exists
    if not mat_file.exists():
        raise FileNotFoundError(f"File '{mat_file}' does not exist.")

    # read the .mat file and supress warnings
    return_dict = {}
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        sio.loadmat(
            file_name=mat_file,
            mdict=return_dict,
            squeeze_me=True,
            variable_names=var_list)

    # extract variables
    return return_dict


# define a function to convert a simulation into multiple data frames
print(extract_mat_data(Path("data/taumap6_2024Jul15-104243.mat"), DEF_VAR_LIST))
