#!/bin/bash

# Directory containing the results files
result_dir="results"

# Check if the directory exists
if [ ! -d "$result_dir" ]; then
    echo "Directory $result_dir does not exist."
    exit 1
fi

# Loop through each file in the results directory
for file in "$result_dir"/*; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Count entries where the number of mismatches is <= 1
        count=$(awk '$5 <= 1' "$file" | wc -l)
        echo "File: $file has $count hits with 1 or fewer mismatches."
    fi
done
