import random
import numpy as np
import matplotlib.pyplot as plt
import random
from MC_functions import *

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
min_radius = 0.3
dim = 2
b_upper = 10 * np.ones(dim)
b_lower = 0 * np.ones(dim)
N = 15 # number of particles
max_ic_its = 120

xmin = [0, 0]
overlapWeight = 500

n_steps = int(1e6)
disp_max = 0.01
T = 5000    #arbitrary for NVT accept/reject in MC
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
p = init_rad_and_loc(N,dim,min_radius,b_lower,b_upper,max_ic_its)
plt.figure(1)
plotCircles(p,b_lower,b_upper,"black")


# Initial energy calculation
E_old = 0
for n in range(len(p)):
    p[n].E = Energy(n,p,xmin,overlapWeight)
    E_old = E_old + p[n].E



# WE TIME STEPPIN NOW BITCHES
# This is using a metropolis algorithm, instead of Ian's point-by-point brute force approach
#   Took each point in order and moved it repeatedly until it reached a lover energy state or reached a
#   maximum number of iterations
t = [0]
PE = [E_old]
reject = 0
acceptUp = 0
acceptDown = 0
kB = 8.617E-5

for i in range(n_steps):
    E_new = 0
    # Select random atom to move
    rand_n = random.randint(0,len(p)-1)
    x_last = p[rand_n].x
    # E_last = p[rand_n].E
    # Displace the random particle a random distance/direction within box of disp_max
    disp_x = np.ones((1, dim))*random.random()*2*disp_max - np.ones((1, dim))*disp_max
    # New position and Energy for the displaced particle
    p[rand_n].x = p[rand_n].x + disp_x
    tempE = []
    E_new = 0
    for n in range(len(p)):
        tempE.append(Energy(n, p, xmin, overlapWeight))
        E_new = E_new + p[n].E
    # Difference between the new and old energy state for the system
    delE = E_new - PE[i]
    # If E increased accept with NVT probability (for arbitrary T)
    # If E decreased accept automatically
    if delE > 0:
        if random.random() > np.exp(-delE/(kB*T)): # Reject
            # print(np.exp(-delE/(kB*T)))
            E_new = PE[i]
            p[rand_n].x = x_last
            reject = reject + 1
        else:
            acceptUp = acceptUp + 1
            for n in range(len(p)):
                p[n].E = tempE[n]
    else:
        acceptDown = acceptDown + 1
        for n in range(len(p)):
            p[n].E = tempE[n]
    PE.append(E_new)
    t.append(i)

print("Increases Rejected = ",reject)
print("Increases Accepted = ",acceptUp)
print("Decreases Accepted = ",acceptDown)

plotCircles(p,b_lower,b_upper,"red")


plt.figure(3)
plt.plot(t,PE)





plt.show()






# plotCircles(p,b_lower,b_upper)
#
#
# plt.show()