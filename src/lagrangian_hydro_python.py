import math
import numpy

import inputs

def calculate_total_energy(
    energy,
    rho,
    u,
    v,
    nodelist,
    element_weightc,
    mass
):
    f = 1.0    # "TWOPI in non axisymmetric"`

    ncells = energy.size

    total_energy = 0.0
    total_ke = 0.0
    total_ie = 0.0

    temp_ke = 0.0
    element_energy = 0.0

    for i in range(ncells):
        temp_ke = 0.0

        for j in range(4):
            node = nodelist[i,j]
            temp_ke += 0.5*rho[i]*element_weightc[i,j]*(u[node]**2 + v[node]**2)*f

        element_energy = mass[i]*energy[i]*f + temp_ke

        total_ke += temp_ke
        total_ie += mass[i]*energy[i]*f
        total_energy += element_energy

    return total_energy, total_ke, total_ie


def get_finite_elements(
        x,
        y,
        nodelist,
        iN,
        idNdx,
        idNdy,
        dNdx,
        dNdy,
        element_weightc
):
    ncells = iN.shape[0]

    for i in range(ncells):
        n = nodelist[i,:]

        # Calculate coefficients for shape function integeral
        a1 = (-x[n[0]] + x[n[1]] + x[n[2]] - x[n[3]])/4.0
        a2 = (x[n[0]] - x[n[1]] + x[n[2]] - x[n[3]])/4.0
        a3 = (-x[n[0]] - x[n[1]] + x[n[2]] + x[n[3]])/4.0

        b1 = (-y[n[0]] + y[n[1]] + y[n[2]] - y[n[3]])/4.0
        b2 = (y[n[0]] - y[n[1]] + y[n[2]] - y[n[3]])/4.0
        b3 = (-y[n[0]] - y[n[1]] + y[n[2]] + y[n[3]])/4.0

        # Shape function integrals
        iN[i,0] = ((3.0*b3-b2)*(3.0*a1-a2)-(3.0*a3-a2)*(3.0*b1-b2))/9.0
        iN[i,1] = ((3.0*b3+b2)*(3.0*a1-a2)-(3.0*a3+a2)*(3.0*b1-b2))/9.0
        iN[i,2] = ((3.0*b3+b2)*(3.0*a1+a2)-(3.0*a3+a2)*(3.0*b1+b2))/9.0
        iN[i,3] = ((3.0*b3-b2)*(3.0*a1+a2)-(3.0*a3-a2)*(3.0*b1+b2))/9.0

        # Partial derivative integral terms
        idNdx[i,0] = -b3+b1
        idNdx[i,1] = b3+b1
        idNdx[i,2] = b3-b1
        idNdx[i,3] = -b3-b1

        idNdy[i,0] = a3-a1
        idNdy[i,1] = -a3-a1
        idNdy[i,2] = -a3+a1
        idNdy[i,3] = a3+a1

        # Partial derivatives
        J = a1*b3 - a3*b1
        for node in range(4):
            dNdx[i,node] = 0.25*idNdx[i,node]/J
            dNdy[i,node] = 0.25*idNdy[i,node]/J

        # Element weights
        for node in range(4):
            element_weightc[i,node] = iN[i,node]


def get_divv(
    u,
    v,
    dNdx,
    dNdy,
    nodelist,
    divv
):
    ncells = dNdx.shape[0]

    divv.fill(0.0)

    for i in range(ncells):
        for j in range(4):
            node = nodelist[i,j]
            divv[i] += u[node]*dNdx[i,j] + v[node]*dNdy[i,j]


def set_soundspeed(
    pressure,
    rho,
    soundspeed,
    material
):
    ncells = pressure.size

    for i in range(ncells):
        if rho[i] < 0.0 or pressure[i] < 0.0:
            soundspeed[i] = 0.0
        else:
            gamma = inputs.gamma[material[i]]
            soundspeed[i] =(gamma*pressure[i]/rho[i])**0.5


def get_q(
    rho,
    soundspeed,
    divv,
    area,
    q
):
    ncells = rho.size


    for i in range(ncells):
        if divv[i] < 0.0:
            dudx = math.sqrt(area[i])*divv[i]
            q[i] = inputs.CQ*rho[i]*dudx**2 + inputs.CL*rho[i]*soundspeed[i]*abs(dudx)
        else:
            q[i] = 0.0


def get_dt(
    area,
    soundspeed,
    rho,
    q,
    dtold,
    time
):
    ncells = area.size
    delta_t = numpy.zeros(ncells)

    dtmin = 1.0
    dt = 0.0
    control = 0

    for i in range(ncells):
        if area[i] < 0.0:
            print("Negative area in cell {}. Area = {}".format(i, area[i]))
            delta_t[i] = 9999.9
        else:
            delta_t[i] = math.sqrt(
                area[i]/max(inputs.rhocutoff, soundspeed[i]**2 + 2.0*q[i]/rho[i])
            )/2.0

        if delta_t[i] < dtmin:
            dtmin = delta_t[i]
            control = i

    if time <= inputs.t0:
        dt = min(dtmin, inputs.dtmax, inputs.dtinit)
    else:
        dt = min(dtmin, inputs.dtmax, dtold*inputs.growth)
        if dt == dtold*inputs.growth:
            control = -1

    return dt, control


def move_nodes(
    dt,
    x,
    y,
    u,
    v,
    xnew,
    ynew
):
    nnodes = x.size

    for i in range(nnodes):
        xnew[i] = x[i] + dt*u[i]
        ynew[i] = y[i] + dt*v[i]


def get_rho(
    mass,
    volume,
    rho
):
    ncells = mass.size

    for i in range(ncells):
        rho[i] = mass[i]/volume[i]


def get_idivv(
    u,
    v,
    nodelist,
    idNdx,
    idNdy,
    idivv
):
    """
    Calculates the summed quantity in equation 3.51 of Andy's thesis
    """
    ncells = nodelist.shape[0]

    idivv.fill(0.0)

    for i in range(ncells):
        for j in range(4):
            node = nodelist[i,j]
            idivv[i] += u[node]*idNdx[i,j] + v[node]*idNdy[i,j]


def get_energy(
    dt,
    pressure,
    q,
    mass,
    energy0,
    idivv,
    energy
):
    ncells = pressure.size

    for i in range(ncells):
        energy[i] = energy0[i] - dt*(pressure[i] + q[i])*idivv[i]/mass[i]


def perform_momentum_update(
    dt,
    u,
    v,
    rho,
    pressure,
    q,
    nodelist,
    iN,
    idNdx,
    idNdy,
    boundary,
    uout,
    vout
):

    nnodes = u.size
    ncells = rho.size

    mstn = numpy.zeros(nnodes, float)
    fstnx = numpy.zeros(nnodes, float)
    fstny = numpy.zeros(nnodes, float)

    #fstnx = force_scatter_to_nodes_x
    #mstn = mass_scatter_to_nodes
    #fstny = force_scatter_to_nodes_y

    # Scatter masses and forces to nodes
    for i in range(ncells):
        for j in range(4):
            node = nodelist[i,j]
            mstn[node] += rho[i]*iN[i,j]
            fstnx[node] += (pressure[i] + q[i])*idNdx[i,j]
            fstny[node] += (pressure[i] + q[i])*idNdy[i,j]

    # F = ma => v = Fdelta_t/m
    for node in range(nnodes):
        uout[node] = u[node] + dt*fstnx[node]/mstn[node]
        vout[node] = v[node] + dt*fstny[node]/mstn[node]

        # Apply boundary conditions
        if boundary[node] == -1 or boundary[node] == -3:
            uout[node] = u[node]

        if boundary[node] == -2 or boundary[node] == -3:
            vout[node] = v[node]


def calculate_area_volume(
    x,
    y,
    nodelist,
    volume,
    area
):
    ncells = volume.size

    for i in range(ncells):
        n = nodelist[i,:]

        a1 = (-x[n[0]] + x[n[1]] + x[n[2]] - x[n[3]])/4.0
        a3 = (-x[n[0]] - x[n[1]] + x[n[2]] + x[n[3]])/4.0

        b1 = (-y[n[0]] + y[n[1]] + y[n[2]] - y[n[3]])/4.0
        b3 = (-y[n[0]] - y[n[1]] + y[n[2]] + y[n[3]])/4.0

        volume[i] = 4.0*(a1*b3 - a3*b1)
        area[i] = volume[i]


def calculate_mass(volume, rho, mass):
    ncells = volume.size

    for i in range(ncells):
        mass[i] = volume[i]*rho[i]
