# Various Constants and Parameters for the MC Code
#
import numpy as np
from MC_functions import *


dim = 2                                           # Dimensions
min_radius = 5                                  # Minimum Particle Radius
b_upper = 100 * np.ones(dim)                      # Upper bound (x,y,z)
b_lower = -100 * np.ones(dim)                        # Lower bound (x,y,z)
N = 15                                            # Number of particles
max_ic_its = 120                                  # Maximum number of tries to generate initial condition
xmin = [0, 0]                                   # Lowest energy location (x,y,z)
overlapWeight = 500                               # Weight for particle overlap penalty
n_steps = 100                                     # Number of MC iterations
it_perParticle = 3                                # Number of iterations to try per particle with Ian's approach
disp_max = 2                                      # Maximum particle displacement distance
pusherTF = False                                  # Whether or not to use the pusher function (W.I.P. atm)

# Text output details
txtName = "mp4Test"
header = True

# Graphing Details
showGraph = True
pltTime = 0.02
# Animation Details
saveAnimation = True
aniType = "mp4"
aniName = "mp4Test"
aniDPI = 200


p = init_rad_and_loc(N,dim,min_radius,b_lower,b_upper,max_ic_its)

for n in range(len(p)):
    p[n].E = Energy(n,p,xmin,overlapWeight)

p = MC_Main(n_steps, p, b_lower, b_upper, it_perParticle, disp_max, xmin, overlapWeight,pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI)


writeText(txtName,p,dim,header,b_lower,b_upper)