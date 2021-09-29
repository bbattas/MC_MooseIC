# Various Constants and Parameters for the MC Code
#
import numpy as np
from MC_functions import *
import ffmpeg
from random import randint
import time

dim = 2                                           # Dimensions
min_radius = 200                                  # Minimum Particle Radius
b_upper = [4000,5000]#[4000,8000,4000]#100 * np.ones(dim)                      # Upper bound [x,y,z]
b_lower = 0 * np.ones(dim)                        # Lower bound [x,y,z]
periodic = True
N = 25                                            # Number of particles
max_ic_its = 120                                  # Maximum number of tries to generate initial condition
xmin = [1000,1000]#0#[0,0,0]#[0, 0]                                   # Lowest energy location [x,y,z], drop uses 1 coordinate
dropAxis = 1                                      # Axis (0,1,2 = x,y,z) for particles to drop (if using drop)
# energyType = "Point"                              # Point, Drop
overlapWeight = 10000                               # Weight for particle overlap penalty
n_steps = 100                                     # Number of MC iterations
it_perParticle = 3                                # Number of iterations to try per particle with Ian's approach
disp_max = 50                                      # Maximum particle displacement distance
pusherTF = False                                  # Whether or not to use the pusher function (W.I.P. atm)

# volumeFractionSampling Details
sampleNum = 10000                                 # Number of points to sample in convex hull for vol frac
samplePlot = True                                # Plot the sampled points on the particles
targetVF = 0.65
vfTol = 0.05

# Text output details
txtName = "p30_drop_2D"
header = False

# Graphing Details
showGraph = True
pltTime = 0.02
# Animation Details
saveAnimation = True
aniType = "gif"
aniName = "periodicTEST_conv2D_domainProblem"
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










p = init_rad_and_loc(N,dim,min_radius,b_lower,b_upper,max_ic_its,periodic)


#
#
# # # Particles all converge to a low energy point (xmin)
p = MC_Main_Point(n_steps, p, dim, b_lower, b_upper, it_perParticle, disp_max, xmin, periodic, overlapWeight,pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI)

# # # Particles fall to y = xmin
# p = MC_Main_Drop(n_steps, p, dim, b_lower, b_upper, it_perParticle, disp_max, xmin, periodic, overlapWeight, dropAxis,
#                  pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI)
#
# # # Turn off interactive plotting
plt.ioff()
# #
# # # Convert Gif output to MP4
# # # stream = ffmpeg.input('mp4Test.gif')
# # # stream = ffmpeg.output(stream, 'mp4Test.mp4')
# # # ffmpeg.run(stream)
# #
# #
# #
# # # Density reduction
# # #   Delete a pore, run a couple steps, repeat
# # #   OR: Delete several and run a couple steps
# # #   OR IN the beginning add small grains that are pores?
# # # Density/volume fraction measurement:
# # #   Either dirty area of hull - area of circles
# # #   OR random points in convex hull? maybe grid but random in hull might be better
# # #
# # # CONSIDER A CONCAVE HULL? (alphashape)
# #
# # # vf = volumeFractionSampling(p,dim,b_lower,b_upper,sampleNum,samplePlot)
# # # vf = vf[0]
# # # maxloop = 20
# # # max_part_del = 10
# # # lp = 0
# # # part_del = 0
# # # while lp < maxloop and part_del < max_part_del:
# # #     lp += 1
# # #     if vf > targetVF + vfTol:
# # #         p_new = np.delete(p, randint(0,len(p)-1),axis=0)
# # #         p_new = MC_Main_Drop(25, p_new, b_lower, b_upper, it_perParticle, disp_max, xmin, overlapWeight,
# # #                              pusherTF,False,pltTime,False,aniName,aniType,aniDPI)
# # #         vf_new = volumeFractionSampling(p_new,dim,b_lower,b_upper,sampleNum,samplePlot=False)
# # #         # If everything is done
# # #         if vf_new <= targetVF + vfTol and vf_new >= targetVF - vfTol:
# # #             vf = vf_new
# # #             p = p_new
# # #             # lp = maxloop + 1
# # #             print("Success- Volume Fraction = ",vf)
# # #             print("* Check the graph to make sure it is one cohesive block of particles")
# # #             break
# # #         elif vf_new > targetVF + vfTol:
# # #             p = p_new
# # #             part_del += 1
# # #             vf = vf_new
# # #             # lp += 1
# # #         # elif vf_new < targetVF - vfTol:
# # #         #     lp += 1
# # #     print("loop = ",lp)
# # #     print("particles deleted = ",part_del)
# # # if part_del == max_part_del:
# # #     print("FAILED: Maximum number of particles deleted")
# vf = volumeFractionSampling(p,dim,b_lower,b_upper,sampleNum,samplePlot)
# print(vf[0])
# #
# #
# #
# # # writeText(txtName,p,dim,header,b_lower,b_upper)
# #
#
# plt.show()
