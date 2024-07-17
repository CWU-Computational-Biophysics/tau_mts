# run_me.py

# imports
from pathlib import Path
from plotter import main as plotter_main
from animator import main as animator_main

# define data dir
data_dir = Path("data")

# plotter
plotter_main(data_dir)

# animator
animator_main(data_dir)
