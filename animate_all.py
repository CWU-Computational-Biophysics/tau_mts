# animate.py
# animate all simulations from a directory


# import packages
import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import colormaps
from rich import print

# import modules
from configurations import DEF_VAR_LIST, DEF_REM_COL_LIST, CSTYLE_FILE_STR
from processing import load_mat_dir
from distributions import create_protein_animation


# get the data directory from the user
user_input = input("Enter the data directory: ")
data_dir = Path(user_input)
if not os.path.isdir(user_input):
    raise ValueError("Invalid data directory")
print(f"Data directory set to '{data_dir}'.")

# configure the figure save directory
FIGURE_DIR = Path("figures", data_dir.parts[-1])
os.makedirs(FIGURE_DIR, exist_ok=True)

# configure the animation save directory
ANIMATION_DIR = Path("animations", data_dir.parts[-1])
os.makedirs(ANIMATION_DIR, exist_ok=True)

# set plotting style
plt.style.use(["default", CSTYLE_FILE_STR])

# set the colormap
DEFAULT_COLORMAP = list(colormaps["tab10"].colors)


# import the data dictionary
# get the data dictionary from the data directory
data_dict = load_mat_dir(data_dir, DEF_VAR_LIST, DEF_REM_COL_LIST)


# add some calculations
# add the final MT length to the sim_df
data_dict["sims"].keys()
length_vals = {}
for sim_name, sim_dict in data_dict["sims"].items():
    # extract the final length as the last entry in first col of mt_length
    length_vals[sim_name] = sim_dict["mt_length"][-1, 0]

# add the dictionary to the df
data_dict["df"]["final_length_units"] = data_dict["df"].index.map(length_vals)

# add final length column
data_dict["df"]["final_length"] = data_dict["df"]["final_length_units"]*data_dict["df"]["dx"]

# create an averages df grouped by tm_ratio
data_avg_df = data_dict["df"].groupby("tm_ratio").mean()


# generate animations of all simulations
# configure animation properties
frame_rate = 30
anim_time = 10
full_frame_rate = 50

# generate animations
sim_names = list(data_dict["sims"].keys())
for sim_name in sim_names:
    # make a low res animation
    create_protein_animation(
        save_path=ANIMATION_DIR / f"{sim_name}.mp4",
        sim_name=sim_name,
        data_dict=data_dict,
        frame_rate=frame_rate,
        anim_time=anim_time,
        overwrite=True,
        protein_points_size=0,
        binding_ticks=False)

    # make a full res animation that shows each time step
    # get the number of time steps
    # time_steps = data_dict["df"].loc[sim_name]["steps"]
    # full_anim_time = int(time_steps / full_frame_rate)

    # make the high rest animation
    # create_protein_animation(
    #     save_path=ANIMATION_DIR / f"{sim_name}_full.mp4",
    #     sim_name=sim_name,
    #     data_dict=data_dict,
    #     frame_rate=full_frame_rate,
    #     anim_time=full_anim_time,
    #     overwrite=True,
    #     protein_points_size=0,
    #     binding_ticks=False)
