#!/bin/bash

# generate nuding files for TVPDS

# 24-hour interval
intervals_per_day=1
# 3-hour interval (24/3)
#intervals_per_day=8
days_state1=4
days_state2=3
state1_file="wrfinput_d01_state1"
state2_file="wrfinput_d01_state2"

# Generate state 1 links
for ((i=1; i<=days_state1*intervals_per_day+1; i++)); do
  ln -s "$state1_file" "wrfinput_d01_t$(printf '%02d' $i)"
done

# Generate state 2 links
for ((j=1; j<=days_state2*intervals_per_day; j++)); do
  k=$((days_state1*intervals_per_day + j+1))
  ln -s "$state2_file" "wrfinput_d01_t$(printf '%02d' $k)"
done

