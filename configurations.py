# configurations.py
# global variables for use in tau_mts

# define the default variable list for extraction
DEF_VAR_LIST = ["alpha_m", "alpha_t", "dt", "dx", "f_cat", "f_res", "m_off", "m_on", "map6_frac_length", "map6_frac_tip", "map6_plus_end_asym", "mt_length", "mt_grid", "state_empty", "state_growing", "state_map6", "state_notexist", "state_shrinking", "state_static", "state_tau", "steps", "t_force", "t_off", "t_on", "tau_frac_length", "tau_frac_tip", "tau_plus_end_asym", "tm_ratio", "ttot", "sim_name", "export_file", "time"]

# define the default column remove list
DEF_REM_COL_LIST = ["__header__", "__version__", "__globals__"]

# define the location of the cstyle.mplstyle sheet
CSTYLE_FILE_STR = "cstyle.mplstyle"
