import random
import numpy as np
import matplotlib.pyplot as plt
from MC_functions import *

min_radius = 0.3
dim = 2
b_upper = 10 * np.ones(dim)
b_lower = 0 * np.ones(dim)
N = 10 # number of particles
max_ic_its = 120

xmin = np.array([0, 0])
overlapWeight = 500

n_steps = int(1e5)
disp_max = 0.01
# T = 5000    #arbitrary for NVT accept/reject in MC
# def f_dist(x,y,n,refPt):
#     p1 = [x[n],y[n]]
#     return np.linalg.norm([x[n],y[n]]-refPt)
#

# Circle Plotting Function
def plotC(r,x,y,b_lower,b_upper,col):
    ax = plt.gca()
    ax.set_xlim((b_lower[0], b_upper[0]))
    ax.set_ylim((b_lower[1], b_upper[1]))
    for n in range(len(r)):
        circle = plt.Circle((x[n],y[n]), r[n], fill=False,color=col)
        ax.add_artist(circle)
    return ax

# Initialize arrays for properties
r = np.zeros(N)
x = np.zeros(N)
y = np.zeros(N)
if dim == 3:
    z = np.zeros(N)
E = np.zeros(N)

# Initial Particle Radius and Location
it = 0
for n in range(N):
    loop = 1
    # RADIUS
    r[n] = part_size_dist(min_radius)
    # LOCATION
    while loop == 1:
        # Define the particle's location
        x[n] = (b_upper[0] - b_lower[0]) * random.random() + b_lower[0]
        y[n] = (b_upper[1] - b_lower[1]) * random.random() + b_lower[1]
        if dim == 3:
            z[n] = (b_upper[2] - b_lower[2]) * random.random() + b_lower[2]
        # Check that the most recent particle isn't overlapping previous particles
        if n != 0:
            for i in range(n):
                # Probably a better wayt to do this part
                if (np.linalg.norm(np.array([x[n],y[n]]) - np.array([x[i],y[i]]))) < (r[n] + r[i]):
                    if it < max_ic_its:
                        it = it + 1
                        loop = 1
                        break
                    else:
                        raise ValueError('Maximum IC Iterations Reached.')
                else:
                    loop = 0
        else:
            loop = 0

print("Iterations for IC generation = ",it)

plt.figure(1)
plotC(r,x,y,b_lower,b_upper,"black")

# INITIAL E CALCULATION
def ijDistXMIN(i,j):
    # return np.linalg.norm([x[i],y[j]] - xmin)
    return ((x[i] - xmin[0])**2 + (y[j] + xmin[1])**2)**0.5

        # np.where(arr==1,0,)
# print(np.linalg.norm([x[2],y[2]] - xmin))

print(np.fromfunction(lambda i, j: ijDistXMIN(i,j), (N, N), dtype=int))
# PE = np.linalg.norm([x[2],y[2]] - xmin) ** 2
# OE = 0







# plt.show()
