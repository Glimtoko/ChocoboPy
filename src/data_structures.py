import numpy as np

class BoundaryConditions():
    def __init__(self, x, y, n):
        self.low = np.ndarray(x, int)
        self.up = np.ndarray(x, int)
        self.left = np.ndarray(y, int)
        self.right = np.ndarray(y, int)
        self.type =  np.ndarray(n, int)
        
        
class Coordinates():
    def __init__(self, n):
        self.x = np.ndarray(n, np.double)
        self.y = np.ndarray(n, np.double)
        self.n = n
        
class Velocities():
    def __init__(self, n):
        self.u = np.ndarray(n, np.double)
        self.v = np.ndarray(n, np.double)
        self.n = n
        
class FEM():
    def __init__(self, n):
        self.iN = np.ndarray((n,4), np.double)
        self.idNdx = np.ndarray((n,4), np.double)
        self.idNdy = np.ndarray((n,4), np.double)
        self.dNdx = np.ndarray((n,4), np.double)
        self.dNdy = np.ndarray((n,4), np.double)
        self.dNdy = np.ndarray((n,4), np.double)
        self.element_weights = np.ndarray((n,4), np.double)
        
        self.divv = np.ndarray(n, np.double)
        self.idivv = np.ndarray(n, np.double)
        
class Mesh():
    def __init__(self, nn, nc, nx, ny):
        self.xcells = nx
        self.ycells = ny
        self.ncells = nc
        self.nnodes = nn
        
        self.coordinates = Coordinates(nn)
        self.coordinates05 = Coordinates(nn)
        self.velocities = Velocities(nn)
        self.fem = FEM(nc)
        
        self.bc = BoundaryConditions(nx+1, ny+1, nn)
        
        self.ρ = np.ndarray(nc, np.double)
        self.ρ05 = np.ndarray(nc, np.double)
        self.pressure = np.ndarray(nc, np.double)
        self.pressure05 = np.ndarray(nc, np.double)
        self.volume = np.ndarray(nc, np.double)
        self.volume05 = np.ndarray(nc, np.double)
        self.soundspeed = np.ndarray(nc, np.double)
        self.area = np.ndarray(nc, np.double)
        self.mass = np.ndarray(nc, np.double)
        self.energy = np.ndarray(nc, np.double)
        self.energy05 = np.ndarray(nc, np.double)
        self.q = np.ndarray(nc, np.double)
        
        self.nodelist = np.ndarray((nc,4), int)
        
        self.regioncelliterator = {}
        self.region = np.ndarray(nc, int)
        self.material = np.ndarray(nc, int)
