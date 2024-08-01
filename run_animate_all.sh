#!/bin/bash

# Define the data directory
DATA_DIR=$1

# Check if the data directory argument is provided
if [ -z "$DATA_DIR" ]; then
  echo "Usage: $0 /path/to/data_directory"
  exit 1
fi

# Run the Python script with nohup
nohup python3 animate_all.py "$DATA_DIR" > animate_all_output.log 2>&1 &

# Get the PID of the last background process
PID=$!

# Print the PID and log file information
echo "Script is running in the background with PID $PID."
echo "Output is being logged to animate_all_output.log."

# Write the PID to a text file
echo $PID >> animate_all_pid.txt
