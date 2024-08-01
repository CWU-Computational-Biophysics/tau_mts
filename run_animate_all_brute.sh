#!/bin/bash

# Define the data directories
DATA_DIRS=("data/paper_50_true" "data/paper_50_false" "data/paper_100_true" "data/paper_100_false")

# Loop through each data directory and run the Python script
for DATA_DIR in "${DATA_DIRS[@]}"; do
  # Define the log file name
  LOG_FILE="output_${DATA_DIR//\//_}.log"
  
  # Run the Python script and redirect output to the log file
  python3 animate_all.py "$DATA_DIR" > "$LOG_FILE" 2>&1 &
  
  # Get the PID of the last background process
  PID=$!
  
  # Use tail to continuously output new lines from the log file to the terminal
  tail -f "$LOG_FILE" &
  
  # Wait for the Python script to finish
  wait $PID
  
  # Kill the tail process after the script finishes
  kill $!
  
  # Check if the script ran successfully
  if [ $? -eq 0 ]; then
    echo "Script for $DATA_DIR completed successfully."
  else
    echo "Script for $DATA_DIR failed."
  fi
done
