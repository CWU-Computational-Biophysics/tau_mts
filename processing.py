# processing.py
# data processing for 'tau_map6.m' simulation '.mat' files

# import packages
import warnings
from os import PathLike
from pathlib import Path

import pandas as pd
import numpy as np
import scipy.io as sio


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


# define a function to load a directory of .mat files
# output should be a dictionary with keys for each simulation
# within those keys are those simulations arrays
# there should also be a key for a combined dataframe for all simulations
def load_mat_dir(mat_dir: PathLike, var_list: list = None, rem_col_list: list = None) -> dict:
    """ Loads a directory of .mat files into a dictionary of variables and arrays.

    Arguments:
        mat_dir -- The path to the directory of .mat files.

    Keyword Arguments:
        var_list -- A list of variables to extract from the .mat files. (default: None, extract all variables)
        rem_col_list -- A list of columns to remove from the DataFrames. (default: None, remove no columns)

    Returns:
        dict: A dictionary of variables and arrays.
    """

    # check that the mat_dir exists
    if not Path(mat_dir).exists():
        raise FileNotFoundError(f"Directory '{mat_dir}' does not exist.")

    # initialize the dictionary
    mat_dir_path = Path(mat_dir)
    return_dict = {
        "mat_dir_path": mat_dir_path,
        "sims": {}
    }

    # initialize the combined DataFrame
    combined_df = pd.DataFrame()

    # iterate over the .mat files in the directory
    for mat_file in mat_dir_path.glob("*.mat"):
        # extract the data from the .mat file
        mat_data = extract_mat_data(mat_file, var_list)

        # convert the data to a DataFrame
        mat_df = mat_data_df(mat_data, rem_col_list)

        # add the DataFrame and array list to the combined DataFrame
        combined_df = pd.concat([combined_df, mat_df["df"]])
        return_dict["sims"][mat_data["sim_name"]] =  {key:val for key, val in mat_df.items() if key != "df"}

    # add the combined DataFrame to the return dictionary
    return_dict["df"] = combined_df

    # return the dictionary
    return return_dict
