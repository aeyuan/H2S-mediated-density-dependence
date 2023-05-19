This repository contains the ODE model for the paper:

## When is an auxotroph not an auxotroph: budding yeast lacking MET17 can collectively overcome their metabolic defect

Below we provide instructions for reproducing each figure in the file "Supplementary_Text_1". Scripts are intended to be run on the command line. For example, to run the file `fit_affinity_param.py`, enter `python fit_affinity_param.py` at the command line.

### Fig ST-1 and ST-2:

Run the script `fit_affinity_param.py`
 * Files for Fig ST-1 will appear as `figures/contour_plot_zoomed_in.pdf` and `figures/contour_plot_zoomed_out.pdf`
 * Files for Fig ST-2 will appear in the folder `figures/check_fit/`

### Fig ST-3:

Run the script `test_model.py`
 * Files for this figure will appear in the folder `figures/leakage_vs_lag/`

### Fig ST-4:

Run the script `pipette_error.py`
 * Files for this figure will appear as `figures/pipette_error.pdf` and `figures/pipette_error_zoom.pdf`
