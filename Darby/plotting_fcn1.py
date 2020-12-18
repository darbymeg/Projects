#Plotting Function 1 - Individual Temperatures
from matplotlib import pyplot as plt
from space_battery_spm_function import ptr


def ind_temps(T_an,T_elyte, T_ca, solution, T_0, T_convert):
    for j,el in enumerate(T_0):
        T_an[j] = solution[j].y[ptr.T_an,:] -  T_convert 
        T_elyte[j] = solution[j].y[ptr.T_elyte,:]-  T_convert 
        T_ca[j] = solution[j].y[ptr.T_ca,:]- T_convert
        
        
        plt.plot(solution[0].t, T_an[j])
        plt.plot(solution[0].t, T_elyte[j])
        plt.plot(solution[0].t, T_ca[j])
        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (C)')
        plt.legend(['Anode temperature','Separator temperature', 'Cathode temperature'])
        plt.show()