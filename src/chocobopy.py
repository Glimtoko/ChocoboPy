import copy
import time as timer

#import inputs
import user_input
import data_structures as ds
import read_mesh as rm
import geometry
import initial_conditions as ic
import lagrangian_hydro as hydro
import eos
import output
 

def chocobopy(meshfile, eosfile, controlfile, usehcmesh, outputloc):  
    inputs = user_input.Inputs()
    inputs.read_input(controlfile)
    
    if usehcmesh:
        nnodes = (inputs.meshx+1)*(inputs.meshy+1)
        ncells = inputs.meshx*inputs.meshy
        
        # Create mesh
        mesh = ds.Mesh(nnodes, ncells, inputs.meshx, inputs.meshy)
        
        # Setup geometry
        geometry.set_geometry_onereg(inputs.x0, inputs.x1, inputs.y0, inputs.y1, mesh)
        
        # Set initial conditions
        gamma = ic.set_initial_conditions_sod(mesh)
    else:
        nnodes, ncells = rm.get_mesh_size(meshfile)
        
        # Create mesh
        mesh = ds.Mesh(nnodes, ncells, 0, 0)
        
        # Read geometry
        rm.read_mesh(meshfile, inputs.material_list, mesh)
        
        # Set initial conditions
        gamma = ic.set_initial_conditions(eosfile, mesh)
        
    # Calculate initial volume
    hydro.calculate_area_volume(
        mesh.coordinates.x, mesh.coordinates.y, mesh.nodelist,
        mesh.volume, mesh.area
    )
    
    # Set initial mass
    hydro.calculate_mass(mesh.volume, mesh.ρ, mesh.mass)
    
    # Set initial pressure from ideal gas EoS
    for i in range(mesh.ncells):
        m = mesh.material[i]
        mesh.energy[i] = mesh.pressure[i]/((gamma[m] - 1.0)*mesh.ρ[i])
    
    # Calculate total energy (initially, KE is zero)
    total_energy, total_ke, total_ie = hydro.calculate_total_energy(
        mesh.energy, mesh.ρ,
        mesh.velocities.u, mesh.velocities.v,
        mesh.nodelist, mesh.fem.element_weights,
        mesh.mass
    )
    
    print("Initial total energy = {}".format(total_energy))
    
    # Main loop - For testing, we'll insert a break
    time = inputs.t0
    dt = inputs.dtinit
    step = 1
    output_dumped = False
    elapsed1 = timer.perf_counter()
    while time <= inputs.tend:
        # Get FEM elements
        hydro.get_finite_elements(
            mesh.coordinates.x, mesh.coordinates.y, mesh.nodelist,
            mesh.fem.iN, mesh.fem.idNdx, mesh.fem.idNdy, mesh.fem.dNdx, mesh.fem.dNdy,
            mesh.fem.element_weights
        )
        
        # Get divergence
        hydro.get_divv(
            mesh.velocities.u, mesh.velocities.v, mesh.fem.dNdx, mesh.fem.dNdy,
            mesh.nodelist, mesh.fem.divv
        )
    
        # Set soundspeed
        hydro.set_soundspeed(mesh.pressure, mesh.ρ, mesh.soundspeed, mesh.material, gamma)
    
        # Calculate artificial viscosity
        hydro.get_q(mesh.ρ, mesh.soundspeed, mesh.fem.divv, mesh.area, mesh.q, inputs.cq, inputs.cl)
    
        # Calculate timestep
        dt, control = hydro.get_dt(
            mesh.area, mesh.soundspeed, mesh.ρ, mesh.q, dt, time, inputs.t0,
            inputs.dtmax, inputs.dtinit, inputs.growth
        )
        time += dt
        print("{:4d}  {:10.7f}  {:15.12f}  {:4d}".format(step, time, dt, control))
        
        # Half timestep positions
        dt05 = dt/2
        hydro.move_nodes(
            dt05, mesh.coordinates.x, mesh.coordinates.y,
            mesh.velocities.u, mesh.velocities.v,
            mesh.coordinates05.x, mesh.coordinates05.y
        )
        
        # Get half-timestep FEM elements
        hydro.get_finite_elements(
            mesh.coordinates05.x, mesh.coordinates05.y, mesh.nodelist,
            mesh.fem.iN, mesh.fem.idNdx, mesh.fem.idNdy, mesh.fem.dNdx, mesh.fem.dNdy,
            mesh.fem.element_weights
        )
        
        # Half-timestep volume
        # volumeold = copy.copy(mesh.volume)
        hydro.calculate_area_volume(
            mesh.coordinates05.x, mesh.coordinates05.y, mesh.nodelist,
            mesh.volume05, mesh.area
        )
        
        # Half-timestep density
        hydro.get_rho(mesh.mass, mesh.volume05, mesh.ρ05)
        
        # Integrate divergence of v
        hydro.get_idivv(
            mesh.velocities.u, mesh.velocities.v, mesh.nodelist,
            mesh.fem.idNdx, mesh.fem.idNdy,
            mesh.fem.idivv
        )
        
        # Half-timestep energy
        hydro.get_energy(
            dt05, mesh.pressure, mesh.q, mesh.mass, mesh.energy,
            mesh.fem.idivv, mesh.energy05
        )
        
        
        # Use EoS to calculate half-timestep pressure
        eos.perfect_gas(mesh.energy05, mesh.ρ05, mesh.pressure05, mesh.material, gamma)
    
        # Save velocities
        uold = copy.copy(mesh.velocities.u)
        vold = copy.copy(mesh.velocities.v)
        
        # Momentum update
        hydro.perform_momentum_update(
            dt, uold, vold, mesh.ρ05, mesh.pressure05, mesh.q,
            mesh.nodelist, mesh.fem.iN, mesh.fem.idNdx, mesh.fem.idNdy,
            mesh.bc.type,
            mesh.velocities.u, mesh.velocities.v
        )
    
        # Time-averaged velocities
        ubar = (uold + mesh.velocities.u)/2
        vbar = (vold + mesh.velocities.v)/2
        
        # Full timestep nodal positions
        hydro.move_nodes(
            dt, mesh.coordinates.x, mesh.coordinates.y,
            ubar, vbar,
            mesh.coordinates.x, mesh.coordinates.y
        )
        
        # Get FEM elements, again
        hydro.get_finite_elements(
            mesh.coordinates.x, mesh.coordinates.y, mesh.nodelist,
            mesh.fem.iN, mesh.fem.idNdx, mesh.fem.idNdy, mesh.fem.dNdx, mesh.fem.dNdy,
            mesh.fem.element_weights
        )
        
        # Full timestep volume
        hydro.calculate_area_volume(
            mesh.coordinates.x, mesh.coordinates.y, mesh.nodelist,
            mesh.volume, mesh.area
        )
    
        # Full-timestep density
        hydro.get_rho(mesh.mass, mesh.volume, mesh.ρ)
        
        # Integrate divergence of v
        hydro.get_idivv(
            ubar, vbar, mesh.nodelist,
            mesh.fem.idNdx, mesh.fem.idNdy,
            mesh.fem.idivv
        )
        
        # Full-timestep energy
        hydro.get_energy(
            dt, mesh.pressure05, mesh.q, mesh.mass, mesh.energy,
            mesh.fem.idivv, mesh.energy
        )
        
        # Use EoS to calculate full-timestep pressure
        eos.perfect_gas(mesh.energy, mesh.ρ, mesh.pressure, mesh.material, gamma)
        
        if (time >= 0.2 and not output_dumped) or inputs.debug_step_count > 0:
            output.output_text(mesh, time, outputloc)
            output_dumped = True
        
        if step == inputs.debug_step_count:
            break
    
        step += 1
        
    elapsed2 = timer.perf_counter()
    
    print("Time taken: {}".format(elapsed2-elapsed1))


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="ChocoboPy CFD Code")
    
    parser.add_argument(
        "--mesh", "-m",
        dest="meshfile",
        default="sod.key", 
        help="File containing mesh"
    )
    
    parser.add_argument(
        "--eos", "-e",
        dest="eosfile",
        default="eos_data.in", 
        help="File containing EoS data"
    )
    
    parser.add_argument(
        "--control", "-c",
        dest="controlfile",
        default="inputs.in", 
        help="File containing input (control) data"
    )

    parser.add_argument(
        "--output", "-o",
        dest="outputloc",
        default="results", 
        help="Directory to put results in. Will be created if necessary"
    )
    
    parser.add_argument(
        "--usehc",
        dest="usehcmesh",
        action="store_true",
        default=False, 
        help="Use hard-coded spherical Sod setup"
    )
    
    
    
    args = parser.parse_args()

    usehcmesh = False
    
    chocobopy(
        args.meshfile, args.eosfile, args.controlfile, args.usehcmesh,
        args.outputloc
    )