#!/bin/bash

oat generate_forces input_pre_relax last_conf.dat
oat forces2pairs forces.txt
oat bond_analysis *input.txt *_trj.dat pairs.txt -p 4
oat mean -p 4 -o mean.dat -d devs.json *_trj.dat
