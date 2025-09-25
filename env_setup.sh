#!/bin/bash

# env_setup.sh
# Defines environment variables and shortcuts to use in project scripts

# This project has a really long name lol, so lets make this easier on us all:
export PROJ_ROOT="$phylogenetics_carnivory_fish_morphology"


# Environment Variables:
export DATA_DIR="$PROJ_ROOT/data"
export RAW_DIR="$DATA_DIR/raw"
export CLEAN_DIR="$DATA_DIR/clean"
export SCRIPTS_DIR="$PROJ_ROOT/scripts"
export RESULTS_DIR="$PROJ_ROOT/results"
export LOGS_DIR="$PROJ_ROOT/logs"

# NOTE: If you close out of your terminal session, you will have to re-run this file in your next session to make these shortcuts useable.
# NOTE: Remember when you are referencing these shortcuts, you have to put the $ before you write the variable name:
	# Example:
		# echo $DATA_DIR
