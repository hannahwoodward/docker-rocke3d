#!/bin/sh

# Load in socrates ENV vars
. /home/app/socrates/set_rad_env

# Run Docker CMD
exec "$@"
