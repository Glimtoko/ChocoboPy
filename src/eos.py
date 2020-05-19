import inputs

def perfect_gas(energy, ρ, pressure, material):
    ncells = energy.size
    
    for i in range(ncells):
        pressure[i] = (inputs.gamma[material[i]] - 1.0)*ρ[i]*energy[i]

