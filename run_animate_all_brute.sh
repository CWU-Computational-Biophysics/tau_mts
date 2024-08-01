#!/bin/bash

# Define the data directories
DATA_DIRS=("data/paper_50_true" "data/paper_50_false" "data/paper_100_true" "data/paper_100_false")

# Loop through each data directory and run the Python script
for DATA_DIR in "${DATA_DIRS[@]}"; do
  # Run the Python script
  python3 animate_all.py "$DATA_DIR" > "output_${DATA_DIR//\//_}.log" 2>&1

  # Check if the script ran successfully
  if [ $? -eq 0 ]; then
    echo "Script for $DATA_DIR completed successfully."
  else
    echo "Script for $DATA_DIR failed."
  fi
done
