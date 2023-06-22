# Networks Gradients Sampling Toolbox

This is a MATLAB toolbox for sampling regional timeseries data with either network or gradient constraints. The toolbox is based on nullspace sampling methods. Although the toolbox is primarily designed with analyiss of neuroimaging data in mind, the methods can be applied to other neural timeseries data.

The toolbox has two main algorithms (functions):
1. sampling timeseries with preserved intra- or all-network correlations (`network_sampler.m`).
2. sampling timeseries with preserved _k_ principal components (`gradient_sampler.m`).

# Getting started

Download or clone this repository into your preferred location. To download, simply press the _Code_ button and select _Download ZIP_.  To clone, enter `git clone https://github.com/AdityaNanda/Networks-Gradients-Sampling-Toolbox` on the command line.

To add the toolbox to the MATLAB path, use the command `addpath(genpath(path_to_toolbox))` or use the _Set path_ button in the _Environment_ section of the _Home_ ribbon, and click _add with subfolders_. Now you can directly access the relevant functions.

# Demo

The scripts `demo_networks.m` and `demo_gradients` include demos of the two main functions, and use the enclosed file `hcp_1subj.mat`, which includes:

* A parcellated regional timeseries for one resting-state recording of one subject from the Human Connectome Project.
* A partition of this vector into 17 * 2 = 34 networks (the parcellation and partition are from Schaefer et al. (2018) https://doi.org/10.1093/cercor/bhx179).

# Licence
This software is free to use for all academic and research purposes. See the License file for details. 

# Reference

Nanda A and Rubinov M (2023) Unbiased and Efficient Sampling of Timeseries Reveals Redundancy of Brain Network and Gradient Structure, _OSF Preprints_, https://doi.org/10.31219/osf.io/ce9xv .
