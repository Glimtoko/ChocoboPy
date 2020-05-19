def set_initial_conditions_sod(mesh):
    # Sod problem initial conditions
    pressure_init = 1.0/10.0
    ρ_init = 1.0/8.0
    uvelocity_init = 0.0
    vvelocity_init = 0.0
    
    # Initial pressure
    mesh.pressure.fill(pressure_init)
    
    # Initial density
    mesh.ρ.fill(ρ_init)
    
    # Initial velocities
    mesh.velocities.u.fill(uvelocity_init)
    mesh.velocities.v.fill(vvelocity_init)
    
    # Bubble
    xbubble = 0.0
    ybubble = 0.0
    
    for index in range(mesh.ncells):
        xbubble = sum(mesh.coordinates.x[mesh.nodelist[index, :]])/4.0
        ybubble = sum(mesh.coordinates.y[mesh.nodelist[index, :]])/4.0
    
        if xbubble**2.0 + ybubble**2.0 <= 0.4**2.0 + 0.000001:
            mesh.ρ[index] = 1.0
            mesh.pressure[index] = 1.0
            
    return [1.4, ]

def set_initial_conditions(eos_file, mesh):
    # import inputs
    
    # Initial velocities
    mesh.velocities.u.fill(0.0)
    mesh.velocities.v.fill(0.0)
    
    gamma_dict = {}
    # Read EOS data
    mat_data = {}
    with open(eos_file, "r") as efile:
        for line in efile:
            line = line.strip()
            if "#" not in line:
                data = line.split(",")
                if len(data) == 4:
                    mat = int(data[0])
                    rho = float(data[1])
                    P = float(data[2])
                    gamma = float(data[3])
                    
                    mat_data[mat] = {"rho":rho, "P":P, "gamma":gamma}
                    
    print(mat_data)
    
    for mat in mat_data.keys():
        gamma_dict[mat] = mat_data[mat]["gamma"]
        for index in range(mesh.ncells):
            if mesh.material[index] == mat:
                mesh.ρ[index] = mat_data[mat]["rho"]
                mesh.pressure[index] = mat_data[mat]["P"]
                
    return gamma_dict