# All the functions to be called in the MC code
#
import random
import numpy as np
import matplotlib.pyplot as plt


class particle:
    def __init__(self):
        self.r = 0
        self.x = []
        self.E = 0

# Particle size distribution function, x = random number 0-1
def part_size_dist(min_radius):
    x = random.random()
    out = 1   #OUTPUT VALUE
    if out < min_radius:
        return min_radius
    else:
        return out

# Initial particle radius and location setup
def init_rad_and_loc(N,dim,min_radius,b_lower,b_upper,max_ic_its):
    p = []
    it = 0
    for n in range(N):
        loop = 1
        p.append(particle())
        # RADIUS
        # Define initial particle radius
        p[n].r = part_size_dist(min_radius)
        # LOCATION
        # Define array for particle centers
        p[n].x = np.zeros((1, dim))
        # loop is the cutoff parameter for when all the points are not overlapping
        while loop == 1:
            # Define the particle's location
            for i in range(dim):
                p[n].x[0, i] = (b_upper[i] - b_lower[i]) * random.random() + b_lower[i]
            # Check that the most recent particle isn't overlapping previous particles
            if n != 0:
                for i in range(0, n):
                    dist = np.linalg.norm(p[n].x - p[i].x)
                    if dist < (p[n].r + p[i].r):
                        if it < max_ic_its:
                            it = it + 1
                            loop = 1
                            break
                        else:
                            print("ERROR: MAX ITERATIONS REACHED")
                    else:
                        loop = 0
            else:
                loop = 0
    print("Iterations used = ", it)
    return p


# Circle Plotting Function
def plotCircles(p,b_lower,b_upper,col):
    ax = plt.gca()
    ax.set_xlim((b_lower[0], b_upper[0]))
    ax.set_ylim((b_lower[1], b_upper[1]))
    for n in range(len(p)):
        circle = plt.Circle(*p[n].x, p[n].r, fill=False,color=col)
        ax.add_artist(circle)

    return ax

# Potential Energy Function
def PotentialE(distFromIdeal):
    # distFromIdeal is the norm between the particle and the lowest energy position
    return distFromIdeal**2

# Calculates the energy of the given particle n in p[n]
def Energy(n,p,xmin,overlapWeight):
    FE_dist = np.linalg.norm(p[n].x - xmin)    # distance from center of particle to minimum energy location
    OE = 0  # overlap energy penalty
    for m in range(len(p)):
        if m != n:
            dist = np.linalg.norm(p[n].x - p[m].x)
            if dist < (p[n].r + p[m].r):
                OE = OE + (p[n].r + p[m].r) - dist
    return PotentialE(FE_dist) + overlapWeight * OE

# # Initial timestep energy calculation for each particle
# def initE(p,xmin,overlapWeight):
#     E_old = 0
#     for n in range(len(p)):
#         p[n].E = Energy(n, p, xmin, overlapWeight)
#         E_old = E_old + p[n].E

