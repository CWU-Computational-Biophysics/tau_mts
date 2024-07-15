# processing.py
# data processing for 'tau_map6.m' simulation '.mat' files

# import packages
import warnings
from os import PathLike
from pathlib import Path

import pandas as pd
import numpy as np
import scipy.io as sio

from configurations import DEF_VAR_LIST, DEF_REM_COL_LIST


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

    # initialize the dictionary
    mat_file_path = Path(mat_file)
    return_dict = {
        "mat_file_path": mat_file_path,
        "sim_name": mat_file_path.stem,
    }

    # define the mat_data keyword arguments
    mat_data_kwargs = {
        "file_name": mat_file,
        "squeeze_me": True,
        "variable_names": var_list
    }

    # read the .mat file to return dict and supress warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        sio.loadmat(mdict=return_dict, **mat_data_kwargs)

    # load with proper dtyping for boolean checks
    mat_dtype_data = sio.loadmat(mat_dtype=True, **mat_data_kwargs)

    # check for booleans in mat_dtype_data and update the return dict value
    for key, value in mat_dtype_data.items():
        if isinstance(value, bool):
            return_dict[key] = value

    # return variables
    return return_dict


# define a function to convert the return dict to a df
def mat_data_df(mat_dict: dict, rem_col_list: list = None) -> pd.DataFrame:
    """ Converts a dictionary of variables and arrays from a .mat file to a pandas DataFrame.

    Arguments:
        mat_dict -- A dictionary of variables and arrays from a .mat file.

    Keyword Arguments:
        rem_col_list -- A list of columns to remove from the DataFrame. (default: None, remove no columns)

    Returns:
        pd.DataFrame: A pandas DataFrame of the variables.
    """

    # seperate the dictionary into an array dictionary and variable dictionary
    array_dict = {}
    var_dict = {}
    for key, value in mat_dict.items():
        # skip keys for faster loading
        if key in rem_col_list: continue

        # keep in the array dict if an array or string
        if isinstance(value, (np.ndarray, str, PathLike)):
            array_dict[key] = value
        else:
            var_dict[key] = [value]

    # convert the var dict into a dataframe
    # the keys should be columns
    var_df = pd.DataFrame(var_dict)

    # set the index to the sim_name from the dictionary
    var_df["sim_name"] = mat_dict["sim_name"]
    var_df.set_index("sim_name", inplace=True)

    # save the df to the array dict under "df"
    array_dict["df"] = var_df

    # return the array dict
    return array_dict


# define a function to convert a simulation into multiple data frames
a = extract_mat_data(Path("data/taumap6_2024Jul15-104243.mat"), DEF_VAR_LIST)
b = mat_data_df(a, DEF_REM_COL_LIST)
print(b["df"])
