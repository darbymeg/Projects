# Inputs:

C_rate = 0.2 # How many charges per hour? Source example from: https://trs.jpl.nasa.gov/bitstream/handle/2014/11058/02-2970.pdf?sequence=1

T_amb = 210.372  # Ambient (surrounding) temperature, K
#https://mars.nasa.gov/all-about-mars/facts/

" Microstructure and geometry "
r_p_an = 4e-6    # Anode particle radius, m
r_p_ca = 0.3e-6  # Cathode particle radius, m

H_an = 66e-6        # Anode thickness, m    Ref: Chang paper
H_elyte = 20e-6     # Electrolyte separator thickness, m  #roughly 1/3rd electrodes
H_ca = 58e-6        # Cathode thickness, m  Ref: Chang paper

eps_graphite = .65     # Volume fraction of graphite in anode
eps_elyte_sep = 0.65   # Volume fraction of elyte in separator
eps_LCO = 0.65         # Volume fraction of LCO phase in cathode

" Electrochemical parameters "
# Initial voltages (used to calculate phi_dl guesses)
phi_an_0 = 0 #V
phi_sep_0 = 1.8  #V
phi_ca_0 = 4.6  #V

##add Gibbs as f(T)
##add U as f(T)

#gibbs_an_elyte = 
#gibbs_ca_elyte= 

#dPhi_eq_an = 
#dPhi_eq_ca = 

# Equilibrium double layer potential (phi_electrode - phi_elyte)
dPhi_eq_an = 0.104   #OG 0.104
dPhi_eq_ca = 3.92    #OG 3.92


# Molar enthalpy of Li species (Li or Li+)
h_Li_an = -10e3      #J/mol
h_Li_elyte = 0       #J/mol
h_Li_ca = -3.9e5     #J/mol

C_dl_an = 1e-2    # Double layer capacitance, F/m2 of dl interface
C_dl_ca = 1e-3    # Double layer capacitance, F/m2 of dl interface

i_o_an = 4.0e-3  # Exchange current density, A/m2
n_an = 1 # Charge equivalents transfered to anode
beta_an = 0.5 # Symmetry parameter


i_o_ca = 1e-3 # Exchange current density, A/m2
n_ca = 1# Charge equivalents transfered to cathode
beta_ca = 0.5# Symmetry parameter


" Material properties "
# Anode - Graphite
density_graphite = 1060 # mass density, kg/m3 2260 ,Ref: Chang paper
capacity_graphite = 370 # Anode charge storage capacity, Ah/kg
#capacity: file:///C:/Users/E9974559/Downloads/energies-12-01074.pdf
Cp_graphite = 688.1 #J/kg-K (for avg. T)
# Cp from https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782425&Mask=2
# Thermal conductivity from: 
# https://tfaws.nasa.gov/wp-content/uploads/TFAWS18-PT-11.pdf
lambda_cond_an = 1.4 #W/m-K
# conductivity taken as an average from:
# https://en.m.wikipedia.org/wiki/Electrical_resistivity_and_conductivity#Resistivity_of_various_materials
sigma_el_graphite = 2e4 #S/m

# Elyte separator
density_elyte = 1132 #kg/m3
Cp_elyte = 1453.43 #J/kg-K 
# Cp from https://webbook.nist.gov/cgi/inchi?ID=C96491&Mask=2#Thermo-Condensed
# Thermal conductivity from: 
# https://tfaws.nasa.gov/wp-content/uploads/TFAWS18-PT-11.pdf
lambda_cond_elyte = 0.1 # W/m-K
# ionic conductivity of elte phase:
sigma_io_elyte = 1.280 #S/m

# Cathode - LCoO2
density_LCO = 1540   # mass density kg/m3 2292 Ref: Chang paper
capacity_LCO = 140  # Cathode charge storage capacity, Ah/kg
#capacity: file:///C:/Users/E9974559/Downloads/energies-12-01074.pdf
# Cp from:
#  https://www.sciencedirect.com/science/article/abs/pii/S0021961414003784
Cp_LCO = 613.7 #J/kg-K (for avg. T)
# Thermal conductivity from: 
# https://tfaws.nasa.gov/wp-content/uploads/TFAWS18-PT-11.pdf
lambda_cond_ca = 0.5 #W/m-K
# ELectronic conductivity:
sigma_el_LCO = 2e6 #S/m

" Heat Transfer Parameters "
# Radiation heat transfer emissivity:
emmissivity = 0.75
# Battery external surface area per unit volume:
A_ext = 1000
# Convection heat transfer coefficient:
h_conv = 6.303 #W/m2-K
# v= 4 m/s on mars,    Ref: Soria

# How deep do we want to charge/discharge?
charge_frac = 0.9