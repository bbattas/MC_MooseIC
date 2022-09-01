import random
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from MC_functions import *
from numba import jit
from datetime import datetime



# def gauss(x,mean,stdev):
#     return np.exp(- (x - mean)**2 / (2*stdev**2))

# Particle size distribution function, x = random number 0-1
def part_size_distrib(min_radius):
    x = random.gauss(0.5,0.3)
    out = 40000 * x + 10000   #OUTPUT VALUE
    if out < min_radius:
        return min_radius
    elif out > 55000:
        return 55000
    else:
        return out

x = np.linspace(0,1,5000)
dist = []
testdist = []
for n in x:
    dist.append(part_size_distrib(8000))
    testdist.append(random.gauss(0.5,0.3))
plt.figure(1)
# plt.plot(x,x*20000*random.gauss(0.5,0.3))
plt.hist(testdist,bins=40)
plt.figure(2)
plt.hist(dist,bins = 10 ** np.linspace(np.log10(0.1E3), np.log10(100E3), 100))
plt.semilogx()
plt.xlabel("Particle Size (nm)")
plt.show()

# print(datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
# @jit(nopython=True)
# class particle:
#     def __init__(self):
#         self.r = 0
#         self.x = []
#         self.E = 0
#
# # Particle size distribution function, x = random number 0-1
# def part_size_dist(min_radius):
#     x = random.random()
#     out = 1   #OUTPUT VALUE
#     if out < min_radius:
#         return min_radius
#     else:
#         return out
#
# def init_rad_and_loc(N,dim,min_radius,b_lower,b_upper,max_ic_its):
#     p = []
#     it = 0
#     for n in range(N):
#         loop = 1
#         p.append(particle())
#         # RADIUS
#         # Define initial particle radius
#         p[n].r = part_size_dist(min_radius)
#         # LOCATION
#         # Define array for particle centers
#         p[n].x = np.zeros((1, dim))
#         # loop is the cutoff parameter for when all the points are not overlapping
#         while loop == 1:
#             # Define the particle's location
#             for i in range(dim):
#                 p[n].x[0, i] = (b_upper[i] - b_lower[i]) * random.random() + b_lower[i]
#             # Check that the most recent particle isn't overlapping previous particles
#             if n != 0:
#                 for i in range(0, n):
#                     dist = np.linalg.norm(p[n].x - p[i].x)
#                     if dist < (p[n].r + p[i].r):
#                         if it < max_ic_its:
#                             it = it + 1
#                             loop = 1
#                             break
#                         else:
#                             print("ERROR: MAX ITERATIONS REACHED")
#                     else:
#                         loop = 0
#             else:
#                 loop = 0
#     return p


# class ar:
#     def __init__(self):
#         self.r = 0
#         self.b = 0
#         self.c = 0
#
# lis = []
# for n in range(4):
#     print(n)
#     lis.append(ar())
#     lis[n].r = n
#     lis[n].b = random.random()
#     print(lis[n])
#     print(vars(lis[n]))
#     print(lis[n].r,lis[n].b,lis[n].c)

# TEMP VALUES
# min_radius = 0.3
# dim = 2
# b_upper = 100 * np.ones(dim)
# b_lower = 0 * np.ones(dim)
# N = 15 # number of particles
# max_ic_its = 120
#
# xmin = [50, 50]
# overlapWeight = 100
#
# n_steps = int(100)
# disp_max = 1
# T = 5000 #arbitrary for NVT accept/reject in MC
#
# it_perParticle = 3
# pusher = False
# p = []
# it = 0
# for n in range(N):
#     loop = 1
#     print("n = ",n)
#     p.append(particle())
#     # RADIUS
#     # Define initial particle radius
#     p[n].r = part_size_dist(min_radius)
#     # LOCATION
#     # Define array for particle centers
#     p[n].x = np.zeros((1,dim))
#     # loop is the cutoff parameter for when all the points are not overlapping
#     while loop == 1:
#         # Define the particle's location
#         for i in range(dim):
#             p[n].x[0,i] = (b_upper[i] - b_lower[i]) * random.random() + b_lower[i]
#         # Check that the most recent particle isn't overlapping previous particles
#         if n != 0:
#             for i in range(0,n):
#                 print("checking ",i)
#                 dist = np.linalg.norm(p[n].x - p[i].x)
#                 print(dist)
#                 if dist < (p[n].r + p[i].r):
#                     if it < max_ic_its:
#                         print("   Particle IS overlapping")
#                         it = it + 1
#                         loop = 1
#                         break
#                     else:
#                         print("ERROR: MAX ITERATIONS REACHED")
#                 else:
#                     loop = 0
#         else:
#             loop = 0
#
#
#
#     print(vars(p[n]))
# print("Iterations used = ",it)
# ax = plt.gca()
# ax.set_xlim((b_lower[0], b_upper[0]))
# ax.set_ylim((b_lower[1], b_upper[1]))
# for n in range(N):
#     circle = plt.Circle(*p[n].x,p[n].r,fill=False)
#     ax.add_artist(circle)

# Particle r and location
# p = init_rad_and_loc(N,dim,min_radius,b_lower,b_upper,max_ic_its)
# plt.figure(1)
# plotCircles(p,b_lower,b_upper,"black")


# Initial energy calculation
# E_old = 0
# for n in range(len(p)):
#     p[n].E = Energy(n,p,xmin,overlapWeight)
#     E_old = E_old + p[n].E

# CAN ONLY CALC ENERGY ON THE MOVED PARTICLE FOR COMPARISON LIKE IAN DID!!!!!

# WE TIME STEPPIN NOW BITCHES
# This is using a metropolis algorithm, instead of Ian's point-by-point brute force approach
#   Ian's took each point in order and moved it repeatedly until it reached a lover energy state or reached a
#   maximum number of iterations
# t = [0]
# PE = [E_old]
# reject = 0
# acceptUp = 0
# acceptDown = 0
# kB = 8.617E-5
#
#
# for i in range(n_steps):
#     E_new = 0
#     # Select random atom to move
#     rand_n = random.randint(0,len(p)-1)
#     x_last = p[rand_n].x
#     # E_last = p[rand_n].E
#     # Displace the random particle a random distance/direction within box of disp_max
#     disp_x = np.ones((1, dim))*random.random()*2*disp_max - np.ones((1, dim))*disp_max
#     # New position and Energy for the displaced particle
#     p[rand_n].x = p[rand_n].x + disp_x
#     tempE = []
#     E_new = 0
#     for n in range(len(p)):
#         tempE.append(Energy(n, p, xmin, overlapWeight))
#         E_new = E_new + p[n].E
#     # Difference between the new and old energy state for the system
#     delE = E_new - PE[i]
#     # If E increased accept with NVT probability (for arbitrary T)
#     # If E decreased accept automatically
#     if delE > 0:
#         if random.random() > np.exp(-delE/(kB*T)): # Reject
#             # print(np.exp(-delE/(kB*T)))
#             E_new = PE[i]
#             p[rand_n].x = x_last
#             reject = reject + 1
#         else:
#             acceptUp = acceptUp + 1
#             for n in range(len(p)):
#                 p[n].E = tempE[n]
#     else:
#         acceptDown = acceptDown + 1
#         for n in range(len(p)):
#             p[n].E = tempE[n]
#     PE.append(E_new)
#     t.append(i)
#
#
# print("Increases Rejected = ",reject)
# print("Increases Accepted = ",acceptUp)
# print("Decreases Accepted = ",acceptDown)


# TAKE 2: ARRAYS EVERYWHERE
# r = np.zeros(len(p))
# x = np.zeros(len(p))
# y = np.zeros(len(p))
# z = np.zeros(len(p))
# # xy = np.zeros(len(p), dtype='i,i')
# # xmin = np.zeros(len(p), dtype='i,i')
# E = np.zeros(len(p))
# for n in range(len(p)):
#     r[n] = p[n].r
#     x[n] = p[n].x[0, 0]
#     y[n] = p[n].x[0, 1]
#     # xy[n] = tuple([x[n],y[n]])
#     # xmin[n] = tuple([0,0])
#     if dim == 3:
#         z[n] = p[n].x[0, 2]
#     E[n] = p[n].E
# xy = [x,y]
# # xmin = np.zeros_like(xy) #NEEDS CHANGED IF I CHANGE FREE ENERGY
# print(EnergyNP(n,r,x,y,xmin,overlapWeight))

# print(r)
# print(E)
# print(E[1])
# print(E[8])
# print(xy)
# print(xmin)

# print(np.subtract(xy[0],xmin[0]))
# print(np.linalg.norm(xy[0]-xmin[0]))
# print(tuple(np.subtract(xy,xmin)))
# print(tuple(np.linalg.norm(np.subtract(xy,xmin))))
# print(EnergyNP(5,r,xy,xmin,overlapWeight))
# for i in range(n_steps):


# Take 3: Ian's Approach
# r in dr is reference to distance bewteen particle centers
# def dE_dr(n,p,xmin,overlapWeight):
#     # FE_dist = np.linalg.norm(p[n].x - xmin)  # distance from center of particle to minimum energy location
#     dPE_dr = 2 * (p[n].x - xmin) #** 2 # derivative of the
#     dOE_dr = 0  # overlap energy penalty
#     for m in range(len(p)):
#         if m != n:
#             dist = np.linalg.norm(p[n].x - p[m].x)
#             if dist < (p[n].r + p[m].r):
#                 dOE_dr = dOE_dr + (p[m].x - p[n].x) / dist
#     return dPE_dr + overlapWeight * dOE_dr
#
# def pusher(n, p, it_perParticle, disp_x, xmin, overlapWeight):
#     p_new = p
#     stop = False
#     i = 0
#     m = -1
#     while stop == False:
#         if i != n:
#             dr = p[n].r + p[i].r
#             dx = np.linalg.norm(p[n].x - p[i].x)
#             if dx < dr:
#                 stop = True
#                 m = i
#         if i == len(p):
#             stop = True
#         i = i + 1
#     # Push the overlapping particle
#     if m != -1:
#         p_new[m].x = p[m].x - disp_x
#         p_new[n].E = Energy(n,p_new,xmin,overlapWeight) #IAN used p not pnew in (), but i think this is right
#         p_new[m].E = Energy(m,p_new,xmin,overlapWeight)
#         dE = p_new[n].E + p_new[m].E - p[n].E - p[m].E
#         if dE > 0 and it_perParticle>1:
#             p_new = pusher(m, p_new, it_perParticle-1, disp_x, xmin, overlapWeight)
#             dE = p_new[n].E + p_new[m].E - p[n].E - p[m].E
#         if dE > 0:
#             p_new = p
#     return p_new


# Create updating Figure

#
# plt.ion()
# fig = plt.figure(1)
#
# moviewriter = ani.ImageMagickWriter()
# # moviewriter.setup(fig, 'my_movie.gif', dpi=100)
# with moviewriter.saving(fig, "writer_test.gif", 100):
#     plotCircles(p,b_lower,b_upper,"red")
#     fig.canvas.draw()
#     moviewriter.grab_frame()
#     plt.pause(0.2)
#     t = [0]
#     PE = [E_old]
# # with moviewriter.saving(fig, "writer_test.gif", 100):
#     for i in range(n_steps):
#         E_new = 0
#         for n in range(len(p)):
#             x_last = p[n].x
#             E_last = p[n].E
#             it = 0
#             stop = False
#             while stop == False:
#                 dE = dE_dr(n,p,xmin,overlapWeight)
#                 disp_x = random.random() * disp_max * dE/np.linalg.norm(dE)
#                 p[n].x = p[n].x - disp_x
#                 p[n].E = Energy(n,p,xmin,overlapWeight)
#                 if p[n].E < E_last:
#                     stop = True
#                 elif it >= it_perParticle:
#                     stop = True
#                     p[n].x = x_last
#                     p[n].E = E_last
#                 elif pusher == True:
#                     p = pusher(n, p, it_perParticle, disp_x, xmin, overlapWeight)
#                 it = it + 1
#         fig.clear(True)
#         plt.title(i+1)
#         plotCircles(p, b_lower, b_upper, "red")
#         fig.canvas.draw()
#         moviewriter.grab_frame()
#         plt.pause(0.02)

# plt.ioff()
# moviewriter.finish()


    # # E_last = p[rand_n].E
    # # Displace the random particle a random distance/direction within box of disp_max
    # disp_x = np.ones((1, dim))*random.random()*2*disp_max - np.ones((1, dim))*disp_max
    # # New position and Energy for the displaced particle
    # p[rand_n].x = p[rand_n].x + disp_x
    # tempE = []
    # E_new = 0
    # for n in range(len(p)):
    #     tempE.append(Energy(n, p, xmin, overlapWeight))
    #     E_new = E_new + p[n].E
    # # Difference between the new and old energy state for the system
    # delE = E_new - PE[i]
    # # If E increased accept with NVT probability (for arbitrary T)
    # # If E decreased accept automatically
    # if delE > 0:
    #     if random.random() > np.exp(-delE/(kB*T)): # Reject
    #         # print(np.exp(-delE/(kB*T)))
    #         E_new = PE[i]
    #         p[rand_n].x = x_last
    #         reject = reject + 1
    #     else:
    #         acceptUp = acceptUp + 1
    #         for n in range(len(p)):
    #             p[n].E = tempE[n]
    # else:
    #     acceptDown = acceptDown + 1
    #     for n in range(len(p)):
    #         p[n].E = tempE[n]
    # PE.append(E_new)
    # t.append(i)






# plotCircles(p,b_lower,b_upper,"red")


# plt.figure(3)
# plt.plot(t,PE)





# plt.show()



#
# plt.pause(2)


# # TExT FILE
# file = open("text_test.txt","w+")
# file.write("Header1")
# file.write("Header 2")
# for i in range(5):
#     arr = str([i, i+1, i+2])
#     file.write(arr + '\n')






# VOLUME FRACTION CUTTING PARTICLES OUT
# vf = volumeFractionSampling(p,dim,b_lower,b_upper,sampleNum,samplePlot)
# vf = vf[0]
# maxloop = 20
# max_part_del = 10
# lp = 0
# part_del = 0
# while lp < maxloop and part_del < max_part_del:
#     lp += 1
#     if vf > targetVF + vfTol:
#         p_new = np.delete(p, randint(0,len(p)-1),axis=0)
#         p_new = MC_Main_Drop(25, p_new, b_lower, b_upper, it_perParticle, disp_max, xmin, overlapWeight,
#                              pusherTF,False,pltTime,False,aniName,aniType,aniDPI)
#         vf_new = volumeFractionSampling(p_new,dim,b_lower,b_upper,sampleNum,samplePlot=False)
#         # If everything is done
#         if vf_new <= targetVF + vfTol and vf_new >= targetVF - vfTol:
#             vf = vf_new
#             p = p_new
#             # lp = maxloop + 1
#             print("Success- Volume Fraction = ",vf)
#             print("* Check the graph to make sure it is one cohesive block of particles")
#             break
#         elif vf_new > targetVF + vfTol:
#             p = p_new
#             part_del += 1
#             vf = vf_new
#             # lp += 1
#         # elif vf_new < targetVF - vfTol:
#         #     lp += 1
#     print("loop = ",lp)
#     print("particles deleted = ",part_del)
# if part_del == max_part_del:
#     print("FAILED: Maximum number of particles deleted")
# vf = volumeFractionSampling(p,dim,b_lower,b_upper,sampleNum,samplePlot)
#
