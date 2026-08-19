#!/bin/bash

# Function to display help
show_help() {
    echo "Usage: $0 [--process] [--help]"
    echo
    echo "Options:"
    echo "  --process      Look for unprocessed sample_data.h5ad in the data dir. Will output processed_data.h5ad."
    echo "  --help, -h     Display this help message and exit."
    echo
    echo "If --process is not specified, the script will attempt to run the visualization on processed_data.h5ad file."
    echo "The app requires 'processed_data.h5ad' to exist in the ./data/ directory."
    echo
    echo "Note: this script assumes your Python environment (conda/venv) is already"
    echo "activated with dependencies installed. See README.md for setup instructions."
    echo
    exit 0
}

# Check for --help or -h as the first argument
if [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
    show_help
fi

# Check if process flag is provided
if [ "$1" == "--process" ]; then
    python3 single_cell_processing.py
    python3 app.py
else
    if [ ! -f "./data/processed_data.h5ad" ]; then
        echo "Error: 'data/processed_data.h5ad' not found. Please run with --process flag first."
        echo
        exit 1
    fi
    python3 app.py
fi