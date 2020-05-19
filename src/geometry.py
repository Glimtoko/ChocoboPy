def set_geometry_onereg(x0, x1, y0, y1, mesh):
    deltax = (x1 - x0)/mesh.xcells
    deltay = (y1 - y0)/mesh.ycells
    
    low = 0
    up = 0
    left = 0
    right = 0
    
    nodeindex = 0
    for j in range(mesh.ycells+1):
       for i in range(mesh.xcells+1):
           mesh.coordinates.x[nodeindex] = x0 + i*deltax
           mesh.coordinates.y[nodeindex] = y0 + j*deltay
    
           # Stuff
           node1 = nodeindex - j
           node2 = nodeindex - j - 1
           node3 = nodeindex - j - 1 - mesh.xcells
           node4 = nodeindex - j - mesh.xcells
    
           # Internal node
           if i != 0 and i != mesh.xcells:
               if j != 0 and j != mesh.ycells:
                   mesh.nodelist[node1, 0] = nodeindex
                   mesh.nodelist[node2, 1] = nodeindex
                   mesh.nodelist[node3, 2] = nodeindex
                   mesh.nodelist[node4, 3] = nodeindex
    
           # Boundary/Corner nodes
           # Lower boundary
           if j == 0:
               if i != 0 and i != mesh.xcells:
                   mesh.nodelist[node2, 1] = nodeindex
                   mesh.nodelist[node1, 0] = nodeindex
                   mesh.bc.low[low] = nodeindex
                   low += 1
                   mesh.bc.type[nodeindex] = -2
    
           # Upper boundary
           if j == mesh.ycells:
               if i != 0 and i != mesh.xcells:
                   mesh.nodelist[node3, 2] = nodeindex
                   mesh.nodelist[node4, 3] = nodeindex
                   mesh.bc.up[up] = nodeindex
                   up += 1
                   mesh.bc.type[nodeindex] = -2
    
           # Left-hand boundary and corners
           if i == 0:
               if j == 0:
                   # Bottom-left corner
                   mesh.nodelist[node1, 0] = nodeindex
                   mesh.bc.left[left] = nodeindex
                   mesh.bc.low[low] = nodeindex
                   low += 1
                   left += 1
                   mesh.bc.type[nodeindex] = -3
               elif j == mesh.ycells:
                   # Top-left corner
                   mesh.nodelist[node4, 3] = nodeindex
                   mesh.bc.left[left] = nodeindex
                   mesh.bc.up[up] = nodeindex
                   up += 1
                   left += 1
                   mesh.bc.type[nodeindex] = -3
               else:
                   # Left-side boundary
                   mesh.nodelist[node4, 3] = nodeindex
                   mesh.nodelist[node1, 0] = nodeindex
                   mesh.bc.left[left] = nodeindex
                   left += 1
                   mesh.bc.type[nodeindex] = -1
    
           # Right-hand boundary and corners
           if i == mesh.xcells:
               if j == 0:
                   # Bottom-right corner
                   mesh.nodelist[node2, 1] = nodeindex
                   mesh.bc.right[right] = nodeindex
                   mesh.bc.low[low] = nodeindex
                   right += 1
                   low += 1
                   mesh.bc.type[nodeindex] = -3
               elif j == mesh.ycells:
                   # Top-right corner
                   mesh.nodelist[node3, 2] = nodeindex
                   mesh.bc.right[right] = nodeindex
                   mesh.bc.up[up] = nodeindex
                   right += 1
                   up += 1
                   mesh.bc.type[nodeindex] = -3
               else:
                   # Right-side boundary
                   mesh.nodelist[node3, 2] = nodeindex
                   mesh.nodelist[node2, 1] = nodeindex
                   mesh.bc.right[right] = nodeindex
                   right += 1
                   mesh.bc.type[nodeindex] = -1

           nodeindex += 1


#def calculate_area_volume(x, y, nodelist, volume, area):
    #ncells = volume.size
    
    #for i in range(ncells):
        #n = nodelist[i,:]
    
        #a1 = (-x[n[0]] + x[n[1]] + x[n[2]] - x[n[3]])/4.0
        #a3 = (-x[n[0]] - x[n[1]] + x[n[2]] + x[n[3]])/4.0
    
        #b1 = (-y[n[0]] + y[n[1]] + y[n[2]] - y[n[3]])/4.0
        #b3 = (-y[n[0]] - y[n[1]] + y[n[2]] + y[n[3]])/4.0
    
        #volume[i] = 4.0*(a1*b3 - a3*b1)
        #area[i] = volume[i]


#def calculate_mass(volume, ρ, mass):
    #ncells = volume.size
    
    #for i in range(ncells):
        #mass[i] = volume[i]*ρ[i]
