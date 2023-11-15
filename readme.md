This repository contains the ODE model for the paper:

## Sonal et al., PLOS Biology, 2023

Below we provide instructions for reproducing each figure in the file "Supplementary_Text_1". Scripts are intended to be run on the command line. For example, to run the file `fit_affinity_param.py`, enter `python fit_affinity_param.py` at the command line.

### Fig ST-1 and ST-2:

Run the script `fit_affinity_param.py` (runtime: ~ 2 min)
 * Files for Fig ST-1 will appear as `figures/contour_plot_zoomed_in.pdf` and `figures/contour_plot_zoomed_out.pdf`
 * Files for Fig ST-2 will appear in the folder `figures/check_fit/`

### Fig ST-3:

Run the script `leak_vs_lag.py` (runtime: ~ 1 sec)
 * Files for this figure will appear in the folder `figures/leakage_vs_lag/`

### Fig ST-4:

Run the script `pipette_error.py` (runtime: ~ 1 sec)
 * Files for this figure will appear as `figures/pipette_error.pdf` and `figures/pipette_error_zoom.pdf`

### Version notes

Scripts were run using the following libraries and versions:

 * Python version 3.10.9
 * NumPy version 1.23.5
 * SciPy version 1.10.0
 * Pandas version 1.5.3
 * Matplotlib version 3.7.0

Runtimes were determined on a 2023 MacBook Pro.
