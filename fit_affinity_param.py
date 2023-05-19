import numpy as np
from scipy.integrate import solve_ivp
from ode_model import prefix, unit, rate_of_change, default_params_with_gas, default_params_without_gas
import pandas as pd
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, colormaps, colors

############################
### DEFINE PLOT SETTINGS ###
############################

# Linear colormap
red = '#cb4154'  # brickred
black = '#000000'  # pure black
cmap_colors = [(0, red), (1, black)]
cmap_name = 'RedBlack'
cmap = colors.LinearSegmentedColormap.from_list(cmap_name, cmap_colors)

# plotting style
rcParams['axes.linewidth'] = 1.5
rcParams['xtick.major.width'] = 1.5
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.width'] = 1.5
rcParams['ytick.major.size'] = 5
font = {'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 12}
rc('font', **font)

########################
### HELPER FUNCTIONS ###
########################

def OD_to_cell_density(OD):
    """
    Converts OD measurements to cellular density in units of cells per ml, using a calibration curve specific to our instrument.

    Args:
      * OD (1d numpy array): optical density measurements.

    Returns: cell density in units of cells per ml.
    """
    OD_corrected = np.copy(OD)
    apply_corr = OD > 0.5
    OD_corrected[apply_corr] = 2.282 * OD[apply_corr] / (2.748 - OD[apply_corr])
    return OD_corrected * 7 * 10 ** 7 / (pfx.milli * u.liter)

def get_sum_sqr_resid(params, time, density, init_H2S, return_predictions=False):
    """Get sum of squared residuals between model and data.
    Args:
      * params (dict): A dicitonary of model parameters
      * time (1d numpy array): Time in units of seconds
      * density (1d numpy array): density in units of cells/liter
      * init_H2S (float): initial aqueous concentration of H2S in units of moles
        per liter.

    Note: initial gas H2S concentration is assumed to be zero.

    Returns: sum((log2(x[i]) - log2(xhat[i]))^2) where x[i] and xhat[i] are the
    observed and model-expected densities at timepoint i
    """
    def f(t, state):
        return rate_of_change(state, params)
    tspan = [0, np.max(time)]
    init_cond = [density[0], init_H2S, 0]
    result = solve_ivp(f, tspan, init_cond, t_eval=time, method='LSODA')
    density_hat = result.y[0,:]
    resid = np.log2(density) - np.log2(density_hat)
    if return_predictions:
        return np.sum(resid**2), density_hat, result.y[1,:]
    return np.sum(resid**2)

#################
### LOAD DATA ###
#################

u = unit()
pfx = prefix()

# read in excel sheet
df = pd.read_excel('20201019_NaHS_rm11.xlsx', sheet_name='setup', header=49, usecols="A:P").iloc[2:26]
# only want NaHS ≤ 50 because higher values are toxic
df = df[df['NaHS concentration'] <= 50]
df = df.drop(df.columns[-3:], axis=1) # get rid of last 3 columns because these are after 100 hours.
# remove trials in which cells failed to grow
df = df[df.values[:,-1] > 1]
# extract relevant variables
init_H2S = df['NaHS concentration'].values * pfx.micro * u.mole / u.liter # Sonal: are the units correct?
density = df.values[:,2:].astype(float)
density = OD_to_cell_density(density) # convert from OD600 to density
time = np.array(df.columns[2:]).astype(float) * u.hour

######################
### FIT PARAMETERS ###
######################

params = deepcopy(default_params_without_gas)
n_trials = density.shape[0]

# define parameter values to try
aff_fine = np.arange(0.1, 15, 0.1)
aff_coarse = np.arange(15,51,1)
r_fine = np.arange(0, 0.7, 0.01)
r_coarse = np.arange(0.7,20.5,0.5)

aff_grid = np.concatenate([aff_fine, aff_coarse] ) * pfx.micro * u.mole / u.liter
r_grid = np.concatenate([r_fine, r_coarse]) * (pfx.femto * u.mole) / (u.cell * u.hour)

# loop over combinations of parameter values and record model error
total_squared_error = np.zeros([aff_grid.size, r_grid.size]) + np.nan
for r_idx, r in enumerate(r_grid):
    for aff_idx, affinity in enumerate(aff_grid):
        params['affinity'] = affinity
        params['release'] = r
        sqr_err = np.sum([get_sum_sqr_resid(params, time, density[i,:], init_H2S[i]) for i in range(n_trials)])
        total_squared_error[aff_idx, r_idx] = sqr_err

# convert to root-mean squared log error (rmsle)
rmsle = np.sqrt(total_squared_error / (time.size * n_trials))

# find the best-fit parameter values
aff_idx, r_idx = np.where(rmsle==np.min(rmsle))
best_aff = aff_grid[aff_idx.item()] / (pfx.micro * u.mole / u.liter)
best_r = r_grid[r_idx.item()] / ((pfx.femto * u.mole) / (u.cell * u.hour))
print(f'model without gas: best-fit k = {np.round(best_aff,3)} µM')
print(f'model without gas: best-fit r = {np.round(best_r,3)} fmole/cell/hr')
print(f'minimum RMSLE: {np.round(np.min(rmsle), 2)}')

################################
### MAKE RMSLE CONTOUR PLOTS ###
################################

# zoomed-in contour plot
small_r_idx = r_grid <= 0.5 * (pfx.femto * u.mole) / (u.cell * u.hour)
small_aff_idx = aff_grid <= 9 * (pfx.micro * u.mole / u.liter)
plt.figure(figsize=[2,2])
X, Y = np.meshgrid(r_grid[small_r_idx] / ((pfx.femto * u.mole) / (u.cell * u.hour)), aff_grid[small_aff_idx] / (pfx.micro * u.mole / u.liter))
contour_plot = plt.contour(X, Y, rmsle[small_aff_idx,:][:,small_r_idx], levels=[0.9, 1, 1.2, 1.4], cmap=cmap)
plt.clabel(contour_plot, inline=True, fontsize=10)
plt.xlabel('release rate\n(fmole/cell/hr)')
plt.xticks([0, 0.25, 0.5])
plt.yticks([0, 4, 8])
plt.xlim([0, 0.5])
plt.ylim([0, 9])
plt.ylabel('affinity (µM)')
plt.title('RMSLE', **font)
plt.plot(best_r, best_aff, 'xr')
plt.savefig('figures/contour_plot_zoomed_in.png', bbox_inches='tight')
plt.savefig('figures/contour_plot_zoomed_in.pdf', bbox_inches='tight')
plt.close()

# zoomed-out contour plot
plt.figure(figsize=[2,2])
X, Y = np.meshgrid(r_grid / ((pfx.femto * u.mole) / (u.cell * u.hour)), aff_grid / (pfx.micro * u.mole / u.liter))
plt.semilogx()
contour_plot = plt.contour(X, Y, rmsle, cmap=cmap, levels=[1.0, 1.4, 2.1, 3.5])
plt.clabel(contour_plot, inline=True, fontsize=10)
contour_plot = plt.contour(X, Y, rmsle, cmap=cmap, levels=[1.0, 1.4, 2.1, 3.5])
plt.xlabel('release rate\n(fmole/cell/hr)')
plt.xlim([0.01, 20])
plt.ylim([0, 50])
plt.ylabel('affinity (µM)')
plt.title('RMSLE', **font)
plt.plot(best_r, best_aff, 'xr')
plt.savefig('figures/contour_plot_zoomed_out.png', bbox_inches='tight')
plt.savefig('figures/contour_plot_zoomed_out.pdf', bbox_inches='tight')
plt.close()

########################################################
### SHOW FIT QUALITY FOR DIFFERENT PARAMETER CHOICES ###
########################################################

for i in range(n_trials):
    plt.figure(figsize=[1.5,1.5])
    plt.semilogy(time / u.hour, density[i,:] / (u.cell/(pfx.milli * u.liter)), '.-k')
    plt.ylim([10**5, 10**9])
    # plot best fit in red
    params['affinity'] = best_aff * pfx.micro * u.mole / u.liter
    params['release'] = best_r * (pfx.femto * u.mole) / (u.cell * u.hour)
    _, density_hat, H2S_hat = get_sum_sqr_resid(params, time, density[i,:], init_H2S[i], return_predictions=True)
    plt.semilogy(time / u.hour, density_hat / (u.cell/(pfx.milli * u.liter)), '.-', color='red')
    # plot lower-affinity in blue
    params['affinity'] = 20 * pfx.micro * u.mole / u.liter
    params['release'] = best_r * (pfx.femto * u.mole) / (u.cell * u.hour)
    _, density_hat, H2S_hat = get_sum_sqr_resid(params, time, density[i,:], init_H2S[i], return_predictions=True)
    plt.semilogy(time / u.hour, density_hat / (u.cell/(pfx.milli * u.liter)), '.-', color='blue')
    plt.xlim([-5, 105])
    plt.xticks([0, 50, 100])
    plt.title(f"{df['NaHS concentration'].values[i]} µM, # {int(df.replicate.values[i])}")
    plt.xlabel('hour')
    plt.ylabel('cells / ml')
    plt.savefig(f'figures/check_fit/{i+1}.png', bbox_inches='tight', dpi=300)
    plt.savefig(f'figures/check_fit/{i+1}.pdf', bbox_inches='tight', dpi=300)
    plt.close()

aff_idx = np.where(aff_grid == best_aff * pfx.micro * u.mole / u.liter)[0][0]
r_idx = np.where(r_grid == best_r * (pfx.femto * u.mole) / (u.cell * u.hour))[0][0]
print(f'For best_fit params: RMSLE = {np.round(rmsle[aff_idx, r_idx],2)}')

aff_idx = np.where(aff_grid == 20 * pfx.micro * u.mole / u.liter)[0][0]
r_idx = np.where(r_grid == best_r * (pfx.femto * u.mole) / (u.cell * u.hour))[0][0]
print(f'For k=20, RMSLE = {np.round(rmsle[aff_idx, r_idx],2)}')
