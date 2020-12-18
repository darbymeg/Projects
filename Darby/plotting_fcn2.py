#Plotting Function 2 - loglog 
from matplotlib import pyplot as plt
from space_battery_spm_function import ptr

def loglog_combined(T_an,T_elyte, T_ca, solution, T_0, T_convert):
    for j,el in enumerate(T_0):
        T_an[j] = solution[j].y[ptr.T_an,:] -  T_convert 
        T_elyte[j] = solution[j].y[ptr.T_elyte,:]-  T_convert 
        T_ca[j] = solution[j].y[ptr.T_ca,:]- T_convert
        
        plt.loglog(solution[0].t, abs(T_an[j]))
        # plt.loglog(solution[0].t, abs(T_ca[j]))
        # plt.loglog(solution[0].t, abs(T_elyte[j]))
        plt.xlabel('Time (s)')
        plt.ylabel('Abs(Temperature) (C)')
    plt.legend(['Mars Ambient', '- 40 C', '- 30 C', '- 20 C', '- 10 C', '0 C'])
       
    plt.show()

