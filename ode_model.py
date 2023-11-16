import numpy as np
from scipy.integrate import solve_ivp

#################################
### DEFINE PREFIXES AND UNITS ###
#################################

"""
Base units are:
 * time: seconds
 * volume: liters
 * chemical mass: moles
 * biological mass: cells
"""

class prefix:
    def __init__(self):
        self.peta = 10**15
        self.tera = 10**12
        self.giga = 10**9
        self.mega = 10**6
        self.kilo = 10**3
        self.centi = 10**-2
        self.milli = 10**-3
        self.micro = 10**-6
        self.nano = 10 **-9
        self.pico = 10**-12
        self.femto = 10**-15

class unit:
    def __init__(self):
        self.hour = 60 * 60
        self.minute = 60
        self.second = 1
        self.cell = 1
        self.mole = 1
        self.liter = 1

####################
### DEFINE MODEL ###
####################

def rate_of_change(state, params):
  # unpack state vector
  x, saq, sgas = state
  # unpack parameters
  gmax =   params['gmax']
  x_loss = params['death']
  r =      params['release']
  c =      params['consumption']
  k =      params['affinity']
  n =      params['hill_coef']
  cap =    params['carrying_capacity']
  saq_out =params['aq_leak']
  # these are for liquid-gas phase transition
  v_liq =  params['liquid_volume']
  v_gas =  params['gas_volume']
  k_evap = params['k_evap']
  k_sol  = params['k_dissolve']
  sgas_out=params['gas_leak']
  float_cutoff = params['float_cutoff'] # handle small-value numerical problems
  # calculate rates of change
  if saq < float_cutoff:
    g = 0
  else:
    g = gmax * saq ** n / (k ** n + saq ** n) * (1 - x / cap)
  x_dot = g * x - x_loss * x
  saq_dot = r * x - c * g * x - k_evap * saq + k_sol * sgas * v_gas / v_liq - saq_out * saq
  sgas_dot = k_evap * saq * v_liq / v_gas - k_sol * sgas - sgas_out * sgas
  return x_dot, saq_dot, sgas_dot

##########################
### DEFAULT PARAMETERS ###
##########################

u = unit()
pfx = prefix()

default_params_with_gas = {
    'gmax' : 0.26 / u.hour,
    'release': 0.37 * (pfx.femto * u.mole) / (u.cell * u.hour),
    'consumption' : 3 * pfx.femto * u.mole / u.cell,
    'affinity' : 2.8 * pfx.micro * u.mole / u.liter,
    'hill_coef' : 1,
    'death': 0 / u.hour,
    'carrying_capacity' : 1.6 * 10**8 * u.cell / (pfx.milli * u.liter),
    'aq_leak' : 0 / u.hour,
    'liquid_volume' : 3 * pfx.milli * u.liter,
    'gas_volume' : 7 * pfx.milli * u.liter,
    'k_evap' : 1.3 * 10**2 / u.second,
    'k_dissolve' : 1 * 10**2 / u.second,
    'gas_leak' : 0 / u.hour,
    'float_cutoff' : 0
}

default_params_without_gas = {
    'gmax' : 0.26 / u.hour,
    'release': 0.39 * (pfx.femto * u.mole) / (u.cell * u.hour),
    'consumption' : 3 * pfx.femto * u.mole / u.cell,
    'affinity' : 7.1 * pfx.micro * u.mole / u.liter,
    'hill_coef' : 1,
    'death': 0 / u.hour,
    'carrying_capacity' : 1.6 * 10**8 * u.cell / (pfx.milli * u.liter),
    'aq_leak' : 0 / u.hour,
    'liquid_volume' : 3 * pfx.milli * u.liter,
    'gas_volume' : 7 * pfx.milli * u.liter,
    'k_evap' : 0 / u.second,
    'k_dissolve' : 0 * 10**2 / u.second,
    'gas_leak' : 0 / u.hour,
    'float_cutoff' : 0
}

#################
### RUN MODEL ###
#################

def run_model(duration=300 * u.hour,
              tstep=1 * u.hour,
              initial_cell_density=10**7 * u.cell / (pfx.milli * u.liter),
              initial_aq_chem_conc=0 * u.mole / u.liter,
              initial_gas_chem_conc=0 * u.mole / u.liter,
              params=default_params_without_gas, solver='LSODA'):
    """
    Run the model of biological and chemical dynamics by solving an in initial-value problem.
    All inputs use base units (seconds, liters, moles, cells), but the returned dictionary does not.

    Args:
      * duration (numeric): The duration of the simulation.
      * tstep (numeric): The time step for reported values.
      * initial_cell_density (numeric): The initial cell density.
      * initial_aq_chem_conc (numeric): The initial concentration of an aqueous chemical species.
      * initial_gas_chem_conc (numeric): The initial concentration of a gaseous chemical species.
      * params (dictionary): Dictionary containing model parameters. Default is default_params_without_gas.
      * solver (str): The solver method to use for solving the ordinary differential equation. Default is 'LSODA'.

    Returns: A dictionary containing the following:
      'time': An array of time points in hours.
      'cell_density': An array of cell densities in cells per milliliter.
      'aq_chem_conc': An array of aqueous H2S concentrations in moles per liter.
      'gas_chem_conc': An array of gaseous H2S concentrations in moles per liter.
    """
    def f(t, state):
        return rate_of_change(state, params)
    tspan = [0, duration]
    t_eval = np.arange(0, duration+tstep, tstep)
    init_cond = [initial_cell_density, initial_aq_chem_conc, initial_gas_chem_conc]
    result = solve_ivp(f, tspan, init_cond, t_eval=t_eval, method=solver)
    result_dict = {}
    result_dict['time'] = result.t / u.hour
    result_dict['cell_density'] = result.y[0,:] * (pfx.milli * u.liter)/ u.cell
    result_dict['aq_chem_conc'] = result.y[1,:] * u.liter / u.mole
    result_dict['gas_chem_conc'] = result.y[2,:] * u.liter / u.mole
    return result_dict
