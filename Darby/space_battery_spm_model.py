""" Battery single particle model."""
from scipy.integrate import solve_ivp
import numpy as np
from space_battery_spm_function import residual

from space_battery_spm_init import phi_dl_an_0, phi_dl_ca_0, T_0, t_final, pars, ptr

def make_SV_0(phi_dl_an_0,phi_dl_ca_0, T_0):
    # Set solution vector
    return np.array([phi_dl_an_0, T_0, T_0, phi_dl_ca_0, T_0])

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
            cond = 1 # Heat transfer via thermal conduciton
            conv = 1 # Heat transfer via external convection
            rad = 1 # Heat tranfer via external radiation
        
    time_span = np.array([0,t_final*t_fac])
    
    for i, el in enumerate(T_0):
        SV_0 = make_SV_0(phi_dl_an_0,phi_dl_ca_0, el)
        
        sol =solve_ivp(lambda t, y: residual(t, y, pars, ptr, thermal_flags), 
        time_span, SV_0, method='BDF',rtol=1e-6, atol=1e-8)
        
        solution.append(sol)
                        
    return solution


if __name__ == '__main__':
    
    # Run the model, using the parameters in the `inputs.py` file:
    solution = cycle()
    
    T_convert = 273   #convert from K to C
    
    T_an = solution.y[ptr.T_an,:] -  T_convert 
    T_elyte = solution.y[ptr.T_elyte,:]-  T_convert 
    T_ca = solution.y[ptr.T_ca,:]- T_convert
    
    # If running this from the command line, make and print some simple plots:
    from matplotlib import pyplot as plt
    
    plt.plot(solution.t, T_an)
    plt.plot(solution.t, T_elyte)
    plt.plot(solution.t, T_ca)
        
    plt.legend(['Anode temperature','Separator temperature', 'Cathode temperature'])
    plt.show()