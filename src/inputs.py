# EoS
gamma = 1.4

# Artificial Viscosity
CL = 0.5
CQ = 0.75

# Mesh Extents
x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 1.0

# Timestep, etc.
t0 = 0.0
tend = 0.205
dtinit = 0.0001
dtmax = 0.0001
growth = 1.02

# Mesh resolution
nregions = 1
meshx = 50
meshy = 50

# Debug settings
debug_step_count = 0   # Set to zero to disable debugging

# Cutoffs
rhocutoff = 1.0e-6

boundary_values = {
    "XLOW": -1,
    "XHIGH": -1,
    "YLOW": -2,
    "YHIGH": -2,
}

material_list = {
    1: 1,
    2: 2
}
