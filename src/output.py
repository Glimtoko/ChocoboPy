def output_text(mesh, time):
    output_text.count += 1
    
    # Physical Quantities
    name = "results/mypre_no{:03d}.txt".format(output_text.count)
    with open(name, "w") as io:
        for i in range(mesh.ncells):
            line = "{:6d}  {:17.12f}  {:17.12f}  {:19.12f}\n".format(
                i+1, mesh.pressure[i], mesh.ρ[i], mesh.energy[i]
            )
            io.write(line)
    
    name = "results/myvel_no{:03d}.txt".format(output_text.count)
    with open(name, "w") as io:
        for i in range(mesh.nnodes):
            line = "{:17.12f}  {:17.12f}  {:17.12f}\n".format(
                time, mesh.velocities.u[i], mesh.velocities.v[i]
            )
            io.write(line)   
    
    name = "results/myx_no{:03d}.txt".format(output_text.count)
    with open(name, "w") as io:
        for i in range(mesh.nnodes):
            line = "{:17.12f}  {:17.12f}\n".format(
                mesh.coordinates.x[i], mesh.coordinates.y[i]
            )
            io.write(line)
    
    name = "results/nodelist_no{:03d}.txt".format(output_text.count)
    with open(name, "w") as io:
        for i in range(mesh.ncells):
            line = "{:6d}  {:6d}  {:6d}  {:6d}\n".format(
                mesh.nodelist[i,0]+1,
                mesh.nodelist[i,1]+1,
                mesh.nodelist[i,2]+1,
                mesh.nodelist[i,3]+1
            )
            io.write(line)

output_text.count = -1