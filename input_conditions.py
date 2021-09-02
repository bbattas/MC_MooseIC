# Various Constants and Parameters for the MC Code
#
import numpy as np
from MC_functions import *

dim = 2                                           # Dimensions
min_radius = 0.3                                  # Minimum Particle Radius
b_upper = 10 * np.ones(dim)                       # Upper bound (x,y,z)
b_lower = 0 * np.ones(dim)                        # Lower bound (x,y,z)
N = 15                                            # Number of particles
max_ic_its = 120                                  # Maximum number of tries to generate initial condition
xmin = [0, 0]                                     # Lowest energy location (x,y,z)
overlapWeight = 500                               # Weight for particle overlap penalty

p = init_rad_and_loc(N,dim,min_radius,b_lower,b_upper,max_ic_its)
E_old = 0
for n in range(len(p)):
    p[n].E = Energy(n,p,xmin,overlapWeight)
    E_old = E_old + p[n].E


plotCircles(p,b_lower,b_upper,"black")

plt.show()
