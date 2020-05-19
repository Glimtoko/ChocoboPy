def perfect_gas(energy, ρ, pressure, material, gamma):
    ncells = energy.size
    
    for i in range(ncells):
        pressure[i] = (gamma[material[i]] - 1.0)*ρ[i]*energy[i]

