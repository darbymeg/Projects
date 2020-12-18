""" Battery single particle model."""
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
import numpy as np
from space_battery_spm_function import residual, ptr, phi_dl_an_0, phi_dl_ca_0, T_0

from space_battery_spm_init import t_final, pars


#Sub-function to make solution vector for new Temps
def make_SV(phi_dl_an_0,phi_dl_ca_0, T_0):
    # Set solution vector
    return np.array([phi_dl_an_0, T_0, T_0, phi_dl_ca_0, T_0])

#Main function to call residuals
def cycle(C_rate = None, thermal_flags=None):
    solution = []       #empty list for solution
    
    if C_rate:
        pars.i_ext *= C_rate/pars.C_rate
        t_fac = pars.C_rate/C_rate
        pars.C_rate = C_rate
    else:
        t_fac = 1.0

    if not thermal_flags:
         #Class of flags to turn certain thermal source/sink terms on 
         # (flag = 1) and off (flag = 0).  Default is 'on'
         class thermal_flags:
            rxn = 1 # heat due to surface reactions
            ohm_el = 1 # ohmic/Joule heating from electron conduction
            ohm_io = 1 # ohmic/Joule heating from ion conduction
            cond = 1 # Heat transfer via thermal conduction
            conv = 1 # Heat transfer via external convection
            rad = 1 # Heat tranfer via external radiation
        
    time_span = np.array([0,t_final*t_fac])
    
    for i, el in enumerate(T_0):
        SV_0 = make_SV(phi_dl_an_0,phi_dl_ca_0, el)
        
        sol =solve_ivp(lambda t, y: residual(t, y, pars, ptr, thermal_flags), 
        time_span, SV_0,method='BDF', t_eval = np.linspace(0,t_final*t_fac,750),rtol=1e-6, atol=1e-8)
        
        solution.append(sol)   #add to list each loop
                        
    return solution


if __name__ == '__main__':
    
    # Run the model, using the parameters in the `inputs.py` file:
    solution = cycle()
    
    T_convert = 273.15   #convert from K to C
    
    #Initalize battery regions 
    T_an = np.zeros((len(T_0),len(solution[0].y[ptr.T_an,:])))
    T_elyte = np.zeros((len(T_0),len(solution[0].y[ptr.T_elyte,:])))
    T_ca = np.zeros((len(T_0),len(solution[0].y[ptr.T_ca,:])))
    
    #Call plotting functions to make temperature plots
    from plotting_fcn1 import ind_temps
    ind_temps(T_an,T_elyte, T_ca, solution, T_0, T_convert)
    
    
    from plotting_fcn2 import loglog_combined
    loglog_combined(T_an,T_elyte, T_ca, solution, T_0, T_convert)





