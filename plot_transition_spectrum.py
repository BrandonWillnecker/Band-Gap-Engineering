import numpy as np
import matplotlib.pyplot as plt

material_properties = [["Si", 1.12, 4.05 ,0.98 ,0.49],
                       ["Ge", 0.67, 4.03, 0.57 ,0.043],
                       ["AlN", 6.28 ,0.6 ,0.4, 3.53],
                       ["AlAs", 3.06 ,3.5 ,0.15 ,0.42],
                       ["GaN" ,3.44, 4.1 ,0.22, 1.4],
                       ["GaP" ,2.26, 3.8 ,1.12, 0.79],
                       ["GaAs", 1.519 ,4.57 ,0.067 ,0.54],
                       ["GaSb", 0.73, 4.06, 0.041, 0.4],
                       ["InN" ,1.9, 4.6, 0.11, 1.63],
                       ["InP" ,1.35, 4.48,0.08 ,0.85],
                       ["InAs", 0.36, 4.9 ,0.023 ,0.41],
                       ["InSb", 0.17, 4.59, 0.014, 0.43],
                       ["CdSe", 1.74, 4.8, 0.13 ,0.45],
                       ["CdSe", 2.42, 4.87, 0.21, 0.8],
                       ["ZnO", 3.37, 4.3, 1.88, 2.9],
                       ["ZnSe", 2.7, 4.09, 0.21 ,0.63],
                       ["ZnTe", 2.30, 3.53, 0.18, 0.65],
                       ["PbS", 0.37, 3.5 ,0.22 ,0.29],
                       ["PbTe", 0.32 ,4.6 ,0.1, 0.2]]

def get_stack_details(stack_str):
    return None

#InN:0.620164, AlAs:2.0681, GaAs:5.33157
layer_sizes = [11.0073,5.03757,13.9841]
layer_types = [3,6,3]

def get_transition_energy_matrix(layer_sizes, layer_types):
    mu = 0.03
    total_size = 0.0
    cumulative_thickness = []
    dz = 0.1
    
    for i in range(len(layer_sizes)):
        print(material_properties[layer_types[i]],layer_sizes[i])
        total_size += layer_sizes[i]
        cumulative_thickness.append(total_size)
        
    PSI_SIZE = int(total_size/dz)
    PSI_SIZE = max(PSI_SIZE,500)
    dz = total_size/PSI_SIZE
    
    H_electron = np.matrix(np.zeros((PSI_SIZE,PSI_SIZE)))
    H_hole = np.matrix(np.zeros((PSI_SIZE,PSI_SIZE)))
    
    V_electron_max = -np.inf
    V_hole_min = np.inf
    
    Ve_values = []
    Vh_values = []
    
    layer_index = 0
    for i in range(PSI_SIZE):
        if i*dz >= cumulative_thickness[layer_index]:
            layer_index +=1
        
        mu_e = -mu/(dz*dz*material_properties[layer_types[layer_index]][3])
        mu_h = mu/(dz*dz*material_properties[layer_types[layer_index]][4])
        
        V_e = -material_properties[layer_types[layer_index]][2]
        V_h = V_e - material_properties[layer_types[layer_index]][1]
        
        Ve_values.append(V_e)
        Vh_values.append(V_h)
        
        if V_e > V_electron_max:
            V_electron_max = V_e
        if V_h < V_hole_min:
            V_hole_min = V_h
            
        if i>=1:
            H_electron[i,i-1] = mu_e
        H_electron[i,i] = -2*mu_e + V_e
        if i+1<PSI_SIZE:
            H_electron[i,i+1] = mu_e
            
        if i>=1:
            H_hole[i,i-1] = mu_h
        H_hole[i,i] = -2*mu_h + V_h
        if i+1<PSI_SIZE:
            H_hole[i,i+1] = mu_h
            
    plt.figure()
    z_values = np.linspace(0,total_size,PSI_SIZE)
    plt.plot(z_values,Ve_values,color="Blue")
    plt.plot(z_values,Vh_values,color="Red")
            
    electron_energy_levels, electron_energy_states = np.linalg.eig(H_electron)
    hole_energy_levels, hole_energy_states = np.linalg.eig(H_hole)
    
    electron_energy_levels_sorted = []
    hole_energy_levels_sorted  = []
    
    for i in range(PSI_SIZE):
        if electron_energy_levels[i] < V_electron_max:
            electron_energy_levels_sorted.append(electron_energy_levels[i])
            plt.plot([0,total_size],[electron_energy_levels[i],electron_energy_levels[i]],color="Black")
        if hole_energy_levels[i] > V_hole_min:
            hole_energy_levels_sorted.append(hole_energy_levels[i])
            plt.plot([0,total_size],[hole_energy_levels[i],hole_energy_levels[i]],color="Black")
    
    
    electron_energy_levels_sorted = np.sort(electron_energy_levels_sorted)
    hole_energy_levels_sorted = np.sort(hole_energy_levels_sorted)
    
    print(electron_energy_levels_sorted)
    print(hole_energy_levels_sorted)
    
    n_electron_states = len(electron_energy_levels_sorted)
    n_hole_states = len(hole_energy_levels_sorted)
    
    transition_energy_matrix = np.zeros((n_electron_states,n_hole_states))
    
    for i in range(n_electron_states):
        for j in range(n_hole_states):
            transition_energy_matrix[i,j] = abs(hole_energy_levels_sorted[j]-electron_energy_levels_sorted[i])

    plt.xlabel("Stack Height (nm)")
    plt.ylabel("Energy (eV)")
    plt.title("AlAs - GaAs - AlAs")
    plt.show()
    
    return transition_energy_matrix


print(get_transition_energy_matrix(layer_sizes,layer_types))