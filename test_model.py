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

############################
### SET PARAMETER VALUES ###
############################

u = unit()
pfx = prefix()

init_densities = np.array([10**5, 10**6, 10**7]) * u.cell/(pfx.milli * u.liter)
leakage_rates = np.array([0, 0.02, 0.04, 0.06]) / u.hour
aff_grid = [ 7.1 * pfx.micro * u.mole / u.liter,
             0.71* pfx.micro * u.mole / u.liter]

#####################################
### PERFORM NUMERICAL INTEGRATION ###
#####################################

counter = 1
for k in aff_grid:
    params = deepcopy(default_params_without_gas)
    params['affinity'] = k
    dynamics = {}
    for init_dens in init_densities:
        for leak_rate in leakage_rates:
            params['aq_leak'] = leak_rate
            result = run_model(initial_cell_density=init_dens,
                          duration=300 * u.hour,
                          tstep=1 * u.hour,
                          initial_aq_chem_conc=0 * pfx.milli * u.mole / u.liter,
                          initial_gas_chem_conc=0,params=params)
            dynamics[(init_dens, leak_rate)] = result['cell_density']
            time = result['time']

    plt.figure(figsize=[2.5,1.8])
    for init_dens in init_densities:
        plt.semilogy(time, dynamics[init_dens, leakage_rates[0]], 'b')
        plt.semilogy(time, dynamics[init_dens, leakage_rates[1]], 'r')
        plt.semilogy(time, dynamics[init_dens, leakage_rates[2]], 'g')
        plt.semilogy(time, dynamics[init_dens, leakage_rates[3]], 'm')
    plt.xlabel('hours')
    plt.ylabel('density (cell/ml)')
    plt.ylim([5 * 10**4, 5 * 10**8])
    plt.yticks([10**5, 10**6, 10**7, 10**8])
    plt.xlim([0, 300])
    plt.legend(['0/hr', '0.02/hr', '0.04/hr', '0.06/hr'])
    plt.title(f"gmax: {params['gmax'] * u.hour}/h, k: {params['affinity'] / (pfx.micro * u.mole / u.liter )} µM", **font)
    plt.savefig(f'figures/leakage_vs_lag/{counter}.png', bbox_inches='tight')
    plt.savefig(f'figures/leakage_vs_lag/{counter}.pdf', bbox_inches='tight')
    plt.close()
    counter += 1
