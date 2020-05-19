import math
import numpy

#import inputs

cimport numpy
cimport cython

from cython.parallel import prange
from libc.math cimport fabs

# Data types
DTYPE = numpy.double
ITYPE = numpy.int

ctypedef numpy.double_t DTYPE_t
ctypedef numpy.int_t ITYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def calculate_total_energy(
    numpy.ndarray[DTYPE_t, ndim=1] energy,
    numpy.ndarray[DTYPE_t, ndim=1] rho,
    numpy.ndarray[DTYPE_t, ndim=1] u,
    numpy.ndarray[DTYPE_t, ndim=1] v,
    numpy.ndarray[ITYPE_t, ndim=2] nodelist,
    numpy.ndarray[DTYPE_t, ndim=2] element_weightc,
    numpy.ndarray[DTYPE_t, ndim=1] mass
):
    cdef double f = 1.0    # "TWOPI in non axisymmetric"`
    
    cdef int ncells = energy.size
    
    cdef double total_energy = 0.0
    cdef double total_ke = 0.0
    cdef double total_ie = 0.0
    
    cdef double temp_ke = 0.0
    cdef double element_energy = 0.0
    
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


@cython.boundscheck(False)
@cython.wraparound(False)
def get_finite_elements(
        numpy.ndarray[DTYPE_t, ndim=1] x,
        numpy.ndarray[DTYPE_t, ndim=1] y,
        numpy.ndarray[ITYPE_t, ndim=2] nodelist,
        numpy.ndarray[DTYPE_t, ndim=2] iN,
        numpy.ndarray[DTYPE_t, ndim=2] idNdx,
        numpy.ndarray[DTYPE_t, ndim=2] idNdy,
        numpy.ndarray[DTYPE_t, ndim=2] dNdx,
        numpy.ndarray[DTYPE_t, ndim=2] dNdy,
        numpy.ndarray[DTYPE_t, ndim=2] element_weightc
):
    cdef int ncells = iN.shape[0]
    
    cdef int node
    cdef numpy.ndarray[ITYPE_t, ndim=1] n
    cdef double a1, a2, a3, b1, b2, b3
    cdef double J
    
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


@cython.boundscheck(False)
@cython.wraparound(False)
def get_divv(
    numpy.ndarray[DTYPE_t, ndim=1] u,
    numpy.ndarray[DTYPE_t, ndim=1] v,
    numpy.ndarray[DTYPE_t, ndim=2] dNdx,
    numpy.ndarray[DTYPE_t, ndim=2] dNdy,
    numpy.ndarray[ITYPE_t, ndim=2] nodelist,
    numpy.ndarray[DTYPE_t, ndim=1] divv
):
    cdef int ncells = dNdx.shape[0]
    cdef int node
    
    divv.fill(0.0)
    
    cdef int i, j
    for i in prange(ncells, nogil=True):
        for j in prange(4):
            node = nodelist[i,j]
            divv[i] += u[node]*dNdx[i,j] + v[node]*dNdy[i,j]


@cython.boundscheck(False)
@cython.wraparound(False)
def set_soundspeed(
    numpy.ndarray[DTYPE_t, ndim=1] pressure,
    numpy.ndarray[DTYPE_t, ndim=1] rho,
    numpy.ndarray[DTYPE_t, ndim=1] soundspeed,
    numpy.ndarray[ITYPE_t, ndim=1] material,
    gamma_dict
):
    cdef int ncells = pressure.size
    cdef double gamma = 0
    
    cdef int i
    for i in range(ncells):
        if rho[i] < 0.0 or pressure[i] < 0.0:
            soundspeed[i] = 0.0
        else:
            gamma = gamma_dict[material[i]]
            soundspeed[i] =(gamma*pressure[i]/rho[i])**0.5


@cython.boundscheck(False)
@cython.wraparound(False)
def get_q(
    numpy.ndarray[DTYPE_t, ndim=1] rho,
    numpy.ndarray[DTYPE_t, ndim=1] soundspeed,
    numpy.ndarray[DTYPE_t, ndim=1] divv,
    numpy.ndarray[DTYPE_t, ndim=1] area,
    numpy.ndarray[DTYPE_t, ndim=1] q,
    double cq,
    double cl
):
    cdef int ncells = rho.size
    cdef double dudx
    
    cdef double CQ = cq
    cdef double CL = cl
    
    cdef int i
    for i in prange(ncells, nogil=True):
        if divv[i] < 0.0:
            dudx = (area[i]**0.5)*divv[i]
            q[i] = CQ*rho[i]*dudx**2 + CL*rho[i]*soundspeed[i]*fabs(dudx)
        else:
            q[i] = 0.0


@cython.boundscheck(False)
@cython.wraparound(False)
def get_dt(
    numpy.ndarray[DTYPE_t, ndim=1] area,
    numpy.ndarray[DTYPE_t, ndim=1] soundspeed,
    numpy.ndarray[DTYPE_t, ndim=1] rho,
    numpy.ndarray[DTYPE_t, ndim=1] q, 
    double dtold, 
    double time,
    double t0,
    double dtmax,
    double dtinit,
    double growth
):
    cdef int ncells = area.size
    cdef numpy.ndarray[DTYPE_t, ndim=1] delta_t = numpy.zeros(ncells)
    
    cdef double dtmin = 1.0
    cdef double dt = 0.0
    cdef int control = 0
    cdef double rhocutoff = 1.0e-6
    
    for i in range(ncells):
        if area[i] < 0.0:
            print("Negative area in cell {}. Area = {}".format(i, area[i]))
            delta_t[i] = 9999.9
        else:
            delta_t[i] = math.sqrt(
                area[i]/max(rhocutoff, soundspeed[i]**2 + 2.0*q[i]/rho[i])
            )/2.0
        
        if delta_t[i] < dtmin:
            dtmin = delta_t[i]
            control = i
    
    if time <= t0:
        dt = min(dtmin, dtmax, dtinit)
    else:
        dt = min(dtmin, dtmax, dtold*growth)
        if dt == dtold*growth:
            control = -1
    
    return dt, control


@cython.boundscheck(False)
@cython.wraparound(False)
def move_nodes(
    double dt,
    numpy.ndarray[DTYPE_t, ndim=1] x,
    numpy.ndarray[DTYPE_t, ndim=1] y,
    numpy.ndarray[DTYPE_t, ndim=1] u,
    numpy.ndarray[DTYPE_t, ndim=1] v,
    numpy.ndarray[DTYPE_t, ndim=1] xnew,
    numpy.ndarray[DTYPE_t, ndim=1] ynew
):
    cdef int nnodes = x.size
    
    cdef int i
    for i in prange(nnodes, nogil=True):
        xnew[i] = x[i] + dt*u[i]
        ynew[i] = y[i] + dt*v[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def get_rho(
    numpy.ndarray[DTYPE_t, ndim=1] mass, 
    numpy.ndarray[DTYPE_t, ndim=1] volume, 
    numpy.ndarray[DTYPE_t, ndim=1] rho
):
    cdef int ncells = mass.size
    
    cdef int i
    for i in prange(ncells, nogil=True):
        rho[i] = mass[i]/volume[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def get_idivv(
    numpy.ndarray[DTYPE_t, ndim=1] u, 
    numpy.ndarray[DTYPE_t, ndim=1] v, 
    numpy.ndarray[ITYPE_t, ndim=2] nodelist, 
    numpy.ndarray[DTYPE_t, ndim=2] idNdx, 
    numpy.ndarray[DTYPE_t, ndim=2] idNdy, 
    numpy.ndarray[DTYPE_t, ndim=1] idivv
):
    """
    Calculates the summed quantity in equation 3.51 of Andy's thesis
    """
    cdef int ncells = nodelist.shape[0]
    cdef int node
    
    idivv.fill(0.0)
    
    cdef int i, j
    for i in prange(ncells, nogil=True):
        for j in prange(4):
            node = nodelist[i,j]
            idivv[i] += u[node]*idNdx[i,j] + v[node]*idNdy[i,j]
            

@cython.boundscheck(False)
@cython.wraparound(False)            
def get_energy(
    double dt,
    numpy.ndarray[DTYPE_t, ndim=1] pressure, 
    numpy.ndarray[DTYPE_t, ndim=1] q, 
    numpy.ndarray[DTYPE_t, ndim=1] mass, 
    numpy.ndarray[DTYPE_t, ndim=1] energy0, 
    numpy.ndarray[DTYPE_t, ndim=1] idivv, 
    numpy.ndarray[DTYPE_t, ndim=1] energy
):
    cdef int ncells = pressure.size
    
    cdef int i
    for i in prange(ncells, nogil=True):
        energy[i] = energy0[i] - dt*(pressure[i] + q[i])*idivv[i]/mass[i]


@cython.boundscheck(False)
@cython.wraparound(False) 
def perform_momentum_update(
    double dt,
    numpy.ndarray[DTYPE_t, ndim=1] u, 
    numpy.ndarray[DTYPE_t, ndim=1] v, 
    numpy.ndarray[DTYPE_t, ndim=1] rho, 
    numpy.ndarray[DTYPE_t, ndim=1] pressure, 
    numpy.ndarray[DTYPE_t, ndim=1] q, 
    numpy.ndarray[ITYPE_t, ndim=2] nodelist, 
    numpy.ndarray[DTYPE_t, ndim=2] iN, 
    numpy.ndarray[DTYPE_t, ndim=2] idNdx, 
    numpy.ndarray[DTYPE_t, ndim=2] idNdy, 
    numpy.ndarray[ITYPE_t, ndim=1] boundary, 
    numpy.ndarray[DTYPE_t, ndim=1] uout, 
    numpy.ndarray[DTYPE_t, ndim=1] vout
):

    cdef int nnodes = u.size
    cdef int ncells = rho.size
    
    cdef numpy.ndarray[DTYPE_t, ndim=1] mass_scatter_to_nodes = numpy.zeros(nnodes, float)
    cdef numpy.ndarray[DTYPE_t, ndim=1] force_scatter_to_nodes_x = numpy.zeros(nnodes, float)
    cdef numpy.ndarray[DTYPE_t, ndim=1] force_scatter_to_nodes_y = numpy.zeros(nnodes, float)
    
    cdef double [:] mstn = mass_scatter_to_nodes
    cdef double [:] fstnx = force_scatter_to_nodes_x
    cdef double [:] fstny = force_scatter_to_nodes_y
    
    cdef int node, i, j
    
    # Scatter masses and forces to nodes
    with nogil:
        for i in prange(ncells):
            for j in prange(4):
                node = nodelist[i,j]
                mstn[node] += rho[i]*iN[i,j]
                fstnx[node] += (pressure[i] + q[i])*idNdx[i,j]
                fstny[node] += (pressure[i] + q[i])*idNdy[i,j]
    
    # F = ma => v = Fdelta_t/m
    for node in prange(nnodes, nogil=True):
        uout[node] = u[node] + dt*fstnx[node]/mstn[node]
        vout[node] = v[node] + dt*fstny[node]/mstn[node]
    
        # Apply boundary conditions
        if boundary[node] == -1 or boundary[node] == -3:
            uout[node] = u[node]

        if boundary[node] == -2 or boundary[node] == -3:
            vout[node] = v[node]


@cython.boundscheck(False)
@cython.wraparound(False) 
def calculate_area_volume(
    numpy.ndarray[DTYPE_t, ndim=1] x,
    numpy.ndarray[DTYPE_t, ndim=1] y,
    numpy.ndarray[ITYPE_t, ndim=2] nodelist,
    numpy.ndarray[DTYPE_t, ndim=1] volume,
    numpy.ndarray[DTYPE_t, ndim=1] area
):
    cdef int ncells = volume.size
    cdef numpy.ndarray[ITYPE_t, ndim=1] n
    cdef double a1, a3, b1, b3
    
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
