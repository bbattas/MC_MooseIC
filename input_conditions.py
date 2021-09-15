# Various Constants and Parameters for the MC Code
#
import numpy as np
from MC_functions import *
import ffmpeg


dim = 2                                           # Dimensions
min_radius = 3                                  # Minimum Particle Radius
b_upper = [100,200]#100 * np.ones(dim)                      # Upper bound (x,y,z)
b_lower = 0 * np.ones(dim)                        # Lower bound (x,y,z)
N = 60                                            # Number of particles
max_ic_its = 120                                  # Maximum number of tries to generate initial condition
xmin = 0#[0, 0]                                   # Lowest energy location (x,y,z), drop uses 1 coordinate
energyType = "Point"                              # Point, Drop
overlapWeight = 500                               # Weight for particle overlap penalty
n_steps = 150                                     # Number of MC iterations
it_perParticle = 3                                # Number of iterations to try per particle with Ian's approach
disp_max = 2                                      # Maximum particle displacement distance
pusherTF = False                                  # Whether or not to use the pusher function (W.I.P. atm)

# Text output details
txtName = "60_Drop"
header = False

# Graphing Details
showGraph = True
pltTime = 0.02
# Animation Details
saveAnimation = True
aniType = "gif"
aniName = "60_Drop"
aniDPI = 200


p = init_rad_and_loc(N,dim,min_radius,b_lower,b_upper,max_ic_its)

# for n in range(len(p)):
#     p[n].E = Energy(n,p,xmin,overlapWeight)

# Particles all converge to a low energy point (xmin)
# p = MC_Main_Point(n_steps, p, b_lower, b_upper, it_perParticle, disp_max, xmin, overlapWeight,pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI)

# Particles fall to y = xmin
p = MC_Main_Drop(n_steps, p, b_lower, b_upper, it_perParticle, disp_max, xmin, overlapWeight,pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI)


# Convert Gif output to MP4
# stream = ffmpeg.input('mp4Test.gif')
# stream = ffmpeg.output(stream, 'mp4Test.mp4')
# ffmpeg.run(stream)


writeText(txtName,p,dim,header,b_lower,b_upper)