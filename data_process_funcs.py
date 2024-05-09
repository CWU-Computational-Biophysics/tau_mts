# reference to data processing functions created in data_process.ipynb
# imports
import warnings
from os import PathLike
from pathlib import Path

import pandas as pd
import numpy as np
import scipy.io as sio

from tqdm import tqdm


# functions
# define a function to extract data from a .mat file
def extract_mat_data(mat_file: Path, value_list: list) -> tuple[dict, dict]:
    """
    Extracts variables and arrays from a .mat file.

    Parameters:
    mat_file (Path): The path to the .mat file.
    value_list (list): A list of strings containing the names of the variables and arrays to extract from the .mat file.

    Returns:
    tuple[dict, dict]: A tuple containing two dictionaries: the first dictionary contains the extracted variables, the second dictionary contains the extracted arrays.
    """

    # check that the mat_file exists
    if not mat_file.exists():
        raise FileNotFoundError(f"File '{mat_file}' does not exist.")

    # check that value_list is not empty
    if not value_list:
        raise ValueError("The value_list cannot be empty.")

    # read the .mat file and supress warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        mat_data = sio.loadmat(mat_file, squeeze_me=True)

    # create dictionaries to store extracted data
    var_dict = {}
    array_dict = {}

    # extract variables, check if the extracted variable is a single item, var, or array
    for value in value_list:
        extracted_value = mat_data.get(value)
        if isinstance(extracted_value, np.ndarray):
            array_dict[value] = extracted_value
        else:
            var_dict[value] = extracted_value

    # return the dictionaries
    return var_dict, array_dict


# define a function to extract the sim dictionary from some .mat files
def extract_sim_dict(mat_files: list[PathLike], value_list: list[str], suppress_warnings: bool = True, show_prog: bool = False) -> dict:
    """
    Extracts the simulation dictionary from a directory.

    Parameters:
    mat_files (list[PathLike]): A list of paths to the .mat files.
    value_list (list[str]): A list of strings containing the names of the variables and arrays to extract from the .mat files.
    suppress_warnings (bool): Whether to suppress warnings.
    show_prog (bool): Whether to show progress bars.

    Returns:
    dict: The extracted simulation dictionary.
    """

    # initialize the sim_dict
    # the key of sim_dict shall be the mat_file save_str or the file name
    # the value of sim_dict shall be a dictionary
    # the value shall have one key, "vars", containing the extracted var_dict
    # the value shall have a key for each extracted array whose key name is the array name
    sim_dict = {}

    for mat_file in tqdm(mat_files, desc="Extracting data", unit="file", disable=~show_prog):
        # extract the data
        var_dict, array_dict = extract_mat_data(mat_file, value_list)

        # get the save_str and check that it is not None
        save_str = var_dict["save_str"]
        if save_str is None:
            save_str = mat_file.stem
            if not suppress_warnings:
                print(f"Warning: save_str is None for '{mat_file}'. Using file name, '{save_str}' as save_str.")

        # check that the save str is not already in the dict
        if save_str in sim_dict:
            raise ValueError(f"Duplicate save_str, '{save_str}', found in sim_dict.")

        # create a dictionary to store the extracted data
        sim_dict[save_str] = {
            "vars": var_dict
        }

        # add the extracted arrays to the dictionary
        for array_name, array_data in array_dict.items():
            # check that the array_name is not already in the sim_dict
            if array_name in sim_dict[save_str]:
                raise ValueError(f"Duplicate array name, '{array_name}', found in sim_dict['{save_str}'].")
            sim_dict[save_str][array_name] = array_data

    # return the sim_dict
    return sim_dict


# define a function to create sim_df from extract_sim_dict
def create_sim_df(sim_dict: dict) -> pd.DataFrame:
    """
    Creates a DataFrame from a simulation dictionary.

    Parameters:
    sim_dict (dict): The simulation dictionary.

    Returns:
    pd.DataFrame: The DataFrame containing the simulation data.
    """

    # create a dataframe of each simulations "vars" whose index is the simulation
    sim_df = pd.DataFrame({save_str: sim_data["vars"] for save_str, sim_data in sim_dict.items()}).T
    # remove the save_str column
    sim_df = sim_df.drop(columns=["save_str"])
    # return the sim_df
    return sim_df


def extract_sim_df(mat_files: list[PathLike], value_list: list[str], **kwargs) -> pd.DataFrame:
    """
    Creates a DataFrame from a simulation dictionary.

    Parameters:
    mat_files (list[PathLike]): A list of paths to the .mat files.
    value_list (list[str]): A list of strings containing the names of the variables and arrays to extract from the .mat files.
    **kwargs: Additional keyword arguments to pass to extract_sim_dict.

    Returns:
    pd.DataFrame: The DataFrame containing the simulation data.
    """

    # extract the variables to the sim_dict
    sim_dict = extract_sim_dict(mat_files, value_list, **kwargs)
    # create sim_df
    sim_df = create_sim_df(sim_dict)
    # return the sim_df
    return sim_df
