# Various Constants and Parameters for the MC Code
#
import numpy as np
from MC_functions import *
import ffmpeg
from random import randint
import time

dim = 3                                           # Dimensions
min_diameter = 400                                  # Minimum Particle Radius (nm)
# b_upper = [2000000,2500000]#[4000,8000,4000]#100 * np.ones(dim)                      # Upper bound [x,y,z]
# min_diameter = 5
max_diameter = 5000
b_upper = [5000,5000,8000]
b_lower = 0 * np.ones(dim)                        # Lower bound [x,y,z]
periodic = False
N = 50                                        # Number of particles
max_ic_its = 120                                  # Maximum number of tries to generate initial condition
xmin = [1000,1000,1000]#0#[1000000,1250000]#0#[0,0,0]#[0, 0]                                   # Lowest energy location [x,y,z], drop uses 1 coordinate
dropAxis = 2                                      # Axis (0,1,2 = x,y,z) for particles to drop (if using drop)
# energyType = "Point"                              # Point, Drop
overlapWeight = 100000000                               # Weight for particle overlap penalty
n_steps = 40                                     # Number of MC iterations
it_perParticle = 3                                # Number of iterations to try per particle with Ian's approach
disp_max = 1000                                      # Maximum particle displacement distance
pusherTF = False                                  # Whether or not to use the pusher function (W.I.P. atm)

# volumeFractionSampling Details
sampleNum = 10000                                 # Number of points to sample in convex hull for vol frac
samplePlot = True                                # Plot the sampled points on the particles
targetVF = 0.65
vfTol = 0.05
red_steps = 5                                      # Number of iterations for MC Main in the particle deletion vol frac
maxloop = 20
max_part_del = 20

# Text output details
txtName = "converge_1Large"
header = False

# Graphing Details
showGraph = True
pltTime = 0.02
# Animation Details
saveAnimation = True
aniType = "gif"
aniName = "3D_corner_1Large"
aniDPI = 400

# Outside bounds move
# diff = p[p1].x + p[p2].x
# print(diff)
# print(((diff<b_lower)|(diff>b_upper)))
#
# pD = b_upper-b_lower # periodic domain lengths
# print(p[p1].x + pD)
# perArray = np.array([[-1,-1]*pD,[-1,0]*pD,[-1,1]*pD,
#                      [0,-1]*pD,[0,0]*pD,[0,1]*pD,
#                      [1,-1]*pD,[1,0]*pD,[1,1]*pD])
# perDist = (list(np.linalg.norm(p[p1].x - (p[p2].x + perArray[n])) for n in range(len(perArray))))
# print(perDist)
# print(min(perDist))
# print(perDist.index(min(perDist)))
# perShift = perArray[perDist.index(min(perDist))]










p = init_rad_and_loc(N,dim,min_diameter,max_diameter,b_lower,b_upper,max_ic_its,periodic)
# plots(p,b_lower,b_upper,dim)
# plt.show()

#
#
# # # Particles all converge to a low energy point (xmin)
p = MC_Main_Point(n_steps, p, dim, b_lower, b_upper, it_perParticle, disp_max, xmin, periodic, overlapWeight,pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI)

# # Particles fall to y = xmin
# p = MC_Main_Drop(n_steps, p, dim, b_lower, b_upper, it_perParticle, disp_max, xmin, periodic, overlapWeight, dropAxis,
#                  pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI)

print("Done Main")
writeText(txtName,p,dim,header,b_lower,b_upper)


#
# # # Turn off interactive plotting
plt.ioff()
# plots(p,b_lower,b_upper,dim)
print("Wrote Predelete")
# #
# # # Convert Gif output to MP4
# # # stream = ffmpeg.input('mp4Test.gif')
# # # stream = ffmpeg.output(stream, 'mp4Test.mp4')
# # # ffmpeg.run(stream)
# #



# p = densityReduction(p,dim,b_lower,b_upper,red_steps,periodic,sampleNum,samplePlot,it_perParticle, disp_max, xmin,
#                      overlapWeight,pusherTF,pltTime,aniName,aniType,aniDPI,targetVF,vfTol,maxloop,max_part_del)
#
# print("reduced")
# writeText(txtName,p,dim,header,b_lower,b_upper)




plt.show()
