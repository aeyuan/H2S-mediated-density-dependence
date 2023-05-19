import numpy as np
from ode_model import prefix, unit, run_model, default_params_with_gas, default_params_without_gas
import matplotlib.pyplot as plt
from copy import deepcopy
from matplotlib import rc, rcParams

############################
### DEFINE PLOT SETTINGS ###
############################

rcParams['axes.linewidth'] = 1.5
rcParams['xtick.major.width'] = 1.5
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.width'] = 1.5
rcParams['ytick.major.size'] = 5
font = {'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 12}
rc('font', **font)

###############################################################
### SIMULATE GROWTH AT SLIGHTLY DIFFERENT INITIAL DENSITIES ###
###############################################################

u = unit()
pfx = prefix()

target_dens = 10**5
init_densities = np.array([target_dens*1.05,
                           target_dens,
                           target_dens*0.95]) * u.cell/(pfx.milli * u.liter)

colors = ['red', 'black', 'cyan']

params = deepcopy(default_params_without_gas)
dynamics = {}
for init_dens in init_densities:
    result = run_model(initial_cell_density=init_dens,
                  duration=300 * u.hour,
                  tstep=1 * u.hour,
                  initial_aq_chem_conc=0 * pfx.milli * u.mole / u.liter,
                  initial_gas_chem_conc=0,params=params)
    dynamics[init_dens] = result['cell_density']
    time = result['time']

plt.figure(figsize=[2.5,1.8])
for i, init_dens in enumerate(init_densities):
    plt.semilogy(time, dynamics[init_dens], colors[i])
plt.xlabel('hours')
plt.ylabel('density (cell/ml)')
plt.ylim([5 * 10**4, 5 * 10**8])
plt.yticks([10**5, 10**6, 10**7, 10**8])
plt.xlim([0, 300])
plt.legend(['$1.05 \\times 10^5$', '$10^5$', '$0.95 \\times 10^5$'])
plt.savefig('figures/pipette_error.pdf', bbox_inches='tight')
plt.close()

plt.figure(figsize=[1,1.8])
for i, init_dens in enumerate(init_densities):
    plt.plot(time, dynamics[init_dens], colors[i], linewidth=2)
plt.xlabel('hours')
plt.ylabel('density (cell/ml)')
plt.ylim([10**5 * 0.93, 10**5 * 1.07])
plt.xlim([0, 3])
plt.xticks([0, 1, 2, 3])
plt.yticks([10**5 * 0.95, 10**5, 10**5 * 1.05])
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig('figures/pipette_error_zoom.pdf', bbox_inches='tight')
plt.close()
