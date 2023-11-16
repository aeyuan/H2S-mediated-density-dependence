import numpy as np
from ode_model import prefix, unit, run_model, default_params_with_gas, default_params_without_gas
import matplotlib.pyplot as plt
from copy import deepcopy
from matplotlib import rc, rcParams

u = unit()
pfx = prefix()

############################
### DEFINE PLOT SETTINGS ###
############################

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "font.weight": "bold",
    "xtick.major.size": 4,
    "ytick.major.size": 5,
    "ytick.minor.size": 3,
    "xtick.major.pad": 2, # gap between ticks and labels
    "ytick.major.pad": 2, # gap between ticks and labels
    "xtick.labelsize": 12,
    "grid.color": "0.8",
    "grid.linestyle": "-",
    "grid.linewidth": 0.5,
    "lines.linewidth": 3,
    "lines.color": "k",
    "figure.figsize": (3.5,2.5),
    "svg.fonttype": "none", # to export svg text as text instead of paths
    "pdf.fonttype": 42, # to export pdf text as text instead of paths
    "figure.autolayout": True
    #"axes.labelpad": 0.2 # gap between tick labels and axis label
})

#########################
### VARIED PARAMETERS ###
#########################

init_densities = np.array([10**4, 10**5, 10**6]) *7* u.cell/(pfx.milli * u.liter)
leakage_rates = np.array([0, 0.02, 0.04, 0.06]) / u.hour

######################################
### PANEL E: Hsu1-dependent growth ###
######################################

dynamics = {}
sulfide = {}
for init_dens in init_densities:
    for leak_rate in leakage_rates:
        params = deepcopy(default_params_without_gas)
        params['aq_leak'] = leak_rate
        result = run_model(initial_cell_density=init_dens,
                      duration=300 * u.hour,
                      tstep=1 * u.hour,
                      initial_aq_chem_conc=0 * pfx.milli * u.mole / u.liter,
                      initial_gas_chem_conc=0,params=params)
        dynamics[('no_gas', init_dens, leak_rate)] = result['cell_density']
        sulfide[('no_gas', init_dens, leak_rate)] = result['aq_chem_conc']
        time = result['time']

plt.figure(figsize=[3.7,2.5])
for init_dens in init_densities:

    plt.semilogy(time, dynamics['no_gas', init_dens, leakage_rates[1]],  '--', color='#ab4500', linewidth=2)
    plt.semilogy(time, dynamics['no_gas', init_dens, leakage_rates[2]], '--', color='#f46200', linewidth=2)
    plt.semilogy(time, dynamics['no_gas', init_dens, leakage_rates[3]], '--', color='#ff964f', linewidth=2)
    plt.semilogy(time, dynamics['no_gas', init_dens, leakage_rates[0]], 'k', linewidth=2)


plt.xlabel('Time (h)')
plt.ylabel('Cell density (/ml)')
#plt.ylim([6 * 10**4, 10**5])
plt.yticks([10**5, 10**6, 10**7, 10**8], fontsize='12')
plt.tick_params(axis ='both', which ='minor', length = 3.5)
plt.xlim([0, 200])
plt.legend(['0.02/h', '0.04/h', '0.06/h', '0/h'])
plt.savefig('figures/Fig5E.svg', bbox_inches='tight')
plt.show()

############################################
### PANEL F: FIFTY-FOLD LOWER K CONSTANT ###
############################################

dynamics = {}
for init_dens in init_densities:
    for leak_rate in leakage_rates:
        params = deepcopy(default_params_without_gas)
        params['aq_leak'] = leak_rate
        params['affinity'] = 7.1 * 0.02 * pfx.micro * u.mole / u.liter # 50-fold lower k constant
        result = run_model(initial_cell_density=init_dens,
                      duration=300 * u.hour,
                      tstep=1 * u.hour,
                      initial_aq_chem_conc=0 * pfx.milli * u.mole / u.liter,
                      initial_gas_chem_conc=0,params=params)
        dynamics[('no_gas', init_dens, leak_rate)] = result['cell_density']
        time = result['time']

plt.figure(figsize=[3.7,2.5])
for init_dens in init_densities:

    plt.semilogy(time, dynamics['no_gas', init_dens, leakage_rates[1]],  '--', color='#ab4500', linewidth=2)
    plt.semilogy(time, dynamics['no_gas', init_dens, leakage_rates[2]], '--', color='#f46200', linewidth=2)
    plt.semilogy(time, dynamics['no_gas', init_dens, leakage_rates[3]], '--', color='#ff964f', linewidth=2)
    plt.semilogy(time, dynamics['no_gas', init_dens, leakage_rates[0]], 'k', linewidth=2)


plt.xlabel('Time (h)')
plt.ylabel('Cell density (/ml)')
#plt.ylim([6 * 10**4, 10**5])
plt.yticks([10**5, 10**6, 10**7, 10**8], fontsize='12')
plt.tick_params(axis ='both', which ='minor', length = 3.5)
plt.xlim([0, 200])
plt.legend(['0.02/h', '0.04/h', '0.06/h', '0/h'])
plt.savefig('figures/Fig5F.svg', bbox_inches='tight')
plt.show()
