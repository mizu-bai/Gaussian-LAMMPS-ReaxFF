# Gaussian-LAMMPS-ReaxFF

## Introduction

For some reason, I have to compare the ReaxFF and DFT results, so I have no choice but to write this glue code.

## Usage

In `src`, there is a script for `CHON-2017_w` forcefield. If one wants to switch to other forcefields, modify all settings about atoms.

In `example`, folder `H2O` is an optimization task while `HCN` is a transition state search task.

**NOTE**: `opt=nomicro` & `freq=num`

## Reference

(1) Gaussian, Inc. (n.d.). _External_. Gaussian.com. Retrieved March 16, 2024, from https://gaussian.com/external/

(2) Zhang, W.; van Duin, A. C. T. Improvement of the ReaxFF Description for Functionalized Hydrocarbon/Water Weak Interactions in the Condensed Phase. _J. Phys. Chem. B_ **2018**, _122_ (14), 4083–4092. https://doi.org/10.1021/acs.jpcb.8b01127.
