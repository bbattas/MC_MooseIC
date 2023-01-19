# All the functions to be called in the MC code
#
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from datetime import datetime
from scipy.spatial import ConvexHull, Delaunay
from numpy.linalg import det
from scipy.stats import dirichlet
from collections import namedtuple
from random import randint
# from particle_distribution import particle_size_dist
import math
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits import mplot3d



class particle:
    def __init__(self):
        self.r = 0
        self.x = []
        self.E = 0

# def gauss(x,mean,stdev):
#     return np.exp(- (x - mean)**2 / (2*stdev**2))

# Particle size distribution function, x = random number 0-1
def part_size_dist(min_diameter,max_diameter):
    x = random.gauss(0.5,0.3)
    out = 40000 * x + 10000   #OUTPUT VALUE
    if out < min_diameter:
        return min_diameter/2
    elif out > max_diameter:
        return max_diameter/2
    else:
        return out/2


# PARTICLE SIZE DIST BUILDING SCRIPTS
# they were running the script particle_distribution.py for no reason when i just imported the function normally
# with from particle_distribution import particle_size_dist
def test_tot_vol(d):
    radii = d/2
    vol = 4 * math.pi * (radii ** 3) / 3
    tot_weight = np.sum(vol)
    return tot_weight

def pdf(x,mu,sigma):
    return (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2)) / (x * sigma * np.sqrt(2 * np.pi)))

# Lognormal mu from average value and standard deviation
def mu(avg,stdev):
    return np.log(avg) - (stdev**2)/2

def part_nums(total_particles,avg1,stdev1,volper1,avg2,stdev2,volper2):
    # Particle size (diameter) arrays for lognormal distributions
    d1 = np.random.lognormal(mu(avg1, stdev1), stdev1, 10000)
    d2 = np.random.lognormal(mu(avg2, stdev2), stdev2, 10000)
    # distribution total volumes for 10k particles each
    vol1 = test_tot_vol(d1)
    vol2 = test_tot_vol(d2)
    num2 = (total_particles * vol1 / volper1) / ((vol2 / volper2) + (vol1 / volper1))
    num1 = total_particles - num2
    return num1, num2

def particle_size_dist(number_particles,min_d,max_d,avg1,stdev1,volper1,avg2,stdev2,volper2):
    small_percent = part_nums(1,avg1,stdev1,volper1,avg2,stdev2,volper2)[0]
    # TEMPORARY: force all to be large
    #small_percent = 0.0
    part_d = np.where(np.random.rand(number_particles) < small_percent, np.random.lognormal(mu(avg1, stdev1), stdev1, number_particles),
                      np.random.lognormal(mu(avg2, stdev2), stdev2, number_particles))
    if min_d != 0.0:
        part_d = np.where(part_d < min_d, np.full_like(part_d, min_d, dtype=np.double),
                          np.where(part_d > max_d, np.full_like(part_d, max_d, dtype=np.double), part_d))
    return part_d





# Periodic shift if x outside domain (set p[n].x = domainShift)
def domainShift(p,n,b_lower,b_upper,periodic):
    if periodic == False:
        return p[n].x
    domainTest = np.any((p[n].x>b_upper)|(p[n].x<b_lower))
    if domainTest == False:
        return p[n].x
    elif domainTest == True:
        overTest = p[n].x > b_upper
        underTest = p[n].x < b_lower
        return p[n].x - overTest*(b_upper-b_lower) + underTest*(b_upper-b_lower)

# Distance function to implement periodic boundary conditions if needed
# p1 and p2 will mainly be the particle x cooridnates (p[n].x)
def distance(p1,p2,dim,b_upper,b_lower,periodic):
    if periodic == False:
        return np.linalg.norm(p1 - p2)
    if periodic == True:
        diff = p1 - p2
        perCheck = np.any(abs(diff) > 0.5*(b_upper-b_lower))
        if perCheck == False:
            return np.linalg.norm(p1 - p2)
        elif perCheck == True:
            if dim == 2:
                pD = b_upper-b_lower # periodic domain lengths
                # make an array of the image cell shifts possible in 2D
                perArray = np.array([[-1,-1]*pD,[-1,0]*pD,[-1,1]*pD,
                                     [0,-1]*pD,[0,0]*pD,[0,1]*pD,
                                     [1,-1]*pD,[1,0]*pD,[1,1]*pD])
                # shift the second particle into all possible image cells
                perDist = (list(np.linalg.norm(p1 - (p2 + perArray[n])) for n in range(len(perArray))))
                return min(perDist)
            elif dim == 3:
                pD = b_upper-b_lower # periodic domain lengths
                # make an array of the image cell shifts possible in 3D
                perArray = np.array([[-1,-1, -1]*pD,[-1,0, -1]*pD,[-1,1, -1]*pD,
                                     [0,-1, -1]*pD,[0,0, -1]*pD,[0,1, -1]*pD,
                                     [1,-1, -1]*pD,[1,0, -1]*pD,[1,1, -1]*pD,
                                     [-1, -1, 0] * pD, [-1, 0, 0] * pD, [-1, 1, 0] * pD,
                                     [0, -1, 0] * pD, [0, 0, 0] * pD, [0, 1, 0] * pD,
                                     [1, -1, 0] * pD, [1, 0, 0] * pD, [1, 1, 0] * pD,
                                     [-1, -1, 1] * pD, [-1, 0, 1] * pD, [-1, 1, 1] * pD,
                                     [0, -1, 1] * pD, [0, 0, 1] * pD, [0, 1, 1] * pD,
                                     [1, -1, 1] * pD, [1, 0, 1] * pD, [1, 1, 1] * pD])
                # shift the second particle into all possible image cells
                perDist = (list(np.linalg.norm(p1 - (p2 + perArray[n])) for n in range(len(perArray))))
                return min(perDist)

# for the free energy gradient to be periodic
def d_distance(p1,p2,dim,b_upper,b_lower,periodic):
    if periodic == False:
        return p1 - p2
    if periodic == True:
        diff = p1 - p2
        perCheck = np.any(abs(diff) > 0.5*(b_upper-b_lower))
        if perCheck == False:
            return p1 - p2
        elif perCheck == True:
            if dim == 2:
                pD = b_upper-b_lower # periodic domain lengths
                # make an array of the image cell shifts possible in 2D
                perArray = np.array([[-1,-1]*pD,[-1,0]*pD,[-1,1]*pD,
                                     [0,-1]*pD,[0,0]*pD,[0,1]*pD,
                                     [1,-1]*pD,[1,0]*pD,[1,1]*pD])
                # shift the second particle into all possible image cells
                perDist = (list(np.linalg.norm(p1 - (p2 + perArray[n])) for n in range(len(perArray))))
                return p1 - (p2 + perArray[perDist.index(min(perDist))])
            elif dim == 3:
                pD = b_upper-b_lower # periodic domain lengths
                # make an array of the image cell shifts possible in 3D
                perArray = np.array([[-1,-1, -1]*pD,[-1,0, -1]*pD,[-1,1, -1]*pD,
                                     [0,-1, -1]*pD,[0,0, -1]*pD,[0,1, -1]*pD,
                                     [1,-1, -1]*pD,[1,0, -1]*pD,[1,1, -1]*pD,
                                     [-1, -1, 0] * pD, [-1, 0, 0] * pD, [-1, 1, 0] * pD,
                                     [0, -1, 0] * pD, [0, 0, 0] * pD, [0, 1, 0] * pD,
                                     [1, -1, 0] * pD, [1, 0, 0] * pD, [1, 1, 0] * pD,
                                     [-1, -1, 1] * pD, [-1, 0, 1] * pD, [-1, 1, 1] * pD,
                                     [0, -1, 1] * pD, [0, 0, 1] * pD, [0, 1, 1] * pD,
                                     [1, -1, 1] * pD, [1, 0, 1] * pD, [1, 1, 1] * pD])
                # shift the second particle into all possible image cells
                perDist = (list(np.linalg.norm(p1 - (p2 + perArray[n])) for n in range(len(perArray))))
                return p1 - (p2 + perArray[perDist.index(min(perDist))])

# Initial particle radius and location setup
def init_rad_and_loc(N,dim,min_diameter,max_diameter,b_lower,b_upper,max_ic_its,periodic):
    p = []
    it = 0
    # Units should all be in nm now
    # particle_size_dist(part_num, min_d, max_d, avg1, stdev1, volper1, avg2, stdev2, volper2)
    part_rad = particle_size_dist(N,min_diameter,max_diameter,700,0.5,0.593,5890,0.5,0.407)
    for n in range(N):
        loop = 1
        p.append(particle())
        # RADIUS
        # Define initial particle radius
        # p[n].r = part_size_dist(min_diameter,max_diameter)
        p[n].r = part_rad[n] / 2
        print(p[n].r)
        # TEMPORARY: To manually force 1 particle to be this size
        if n == 0:
            p[n].r = 2500#5000#2500 #particle size = 2x radius
        # LOCATION
        # Define array for particle centers
        p[n].x = np.zeros((1, dim))
        # loop is the cutoff parameter for when all the points are not overlapping
        while loop == 1:
            # Define the particle's location
            for i in range(dim):
                p[n].x[0, i] = (b_upper[i] - b_lower[i]) * random.random() + b_lower[i]
            # TEMPORARY: To manually force largest particel to be on the minimum energy coord
            if n == 0:
                p[n].x[0] = [0, 0]#, 0]
            # Check that the most recent particle isn't overlapping previous particles
            if n != 0:
                for i in range(0, n):
                    # dist = np.linalg.norm(p[n].x - p[i].x)
                    dist = distance(p[n].x,p[i].x,dim,b_upper,b_lower,periodic)
                    if dist < (p[n].r + p[i].r):
                        if it < max_ic_its:
                            it = it + 1
                            loop = 1
                            break
                        else:
                            print("ERROR: MAX ITERATIONS REACHED - there will be particle overlap")
                    else:
                        loop = 0
            else:
                loop = 0
    print("Iterations used for initial condition = ", it)
    print(" ")
    return p

def plots(p,b_lower,b_upper,dim):
    if dim == 2:
        return plotCircles(p,b_lower,b_upper)
    elif dim == 3:
        return plotSpheres(p, b_lower, b_upper)

# Circle Plotting Function
def plotCircles(p,b_lower,b_upper):
    ax = plt.gca()
    ax.set_xlim((b_lower[0], b_upper[0]))
    ax.set_ylim((b_lower[1], b_upper[1]))
    ax.set_aspect('equal')
    for n in range(len(p)):
        circle = plt.Circle(*p[n].x, p[n].r, fill=False)#,color=col
        ax.add_artist(circle)
    return ax

# Sphere plotting function (for dim = 3)
def plotSpheres(p,b_lower,b_upper):
    ax = plt.subplot(1,1,1,projection='3d')#axes(projection='3d')
    # fig = plt.figure()#figsize=(8, 6)
    # ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.set_xlim3d((b_lower[0], b_upper[0]))
    ax.set_ylim3d((b_lower[1], b_upper[1]))
    ax.set_zlim3d((b_lower[2], b_upper[2]))
    ax.set_box_aspect((b_upper[0]-b_lower[0],b_upper[1]-b_lower[1],b_upper[2]-b_lower[2]))
    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:20j] # changing number before j smooths spheres
    for n in range(len(p)):
        x = p[n].r * np.cos(u) * np.sin(v)
        y = p[n].r * np.sin(u) * np.sin(v)
        z = p[n].r * np.cos(v)
        ax.plot_surface(x + p[n].x[0,0], y + p[n].x[0,1], z + p[n].x[0,2],alpha=0.5) # alpha is transparency
    return ax

# Potential Energy Function
def PotentialE(distFromIdeal):
    # distFromIdeal is the norm between the particle and the lowest energy position
    return distFromIdeal**2

# Calculates the energy of the given particle n in p[n]
def Energy(n,p,xmin,overlapWeight,dim,b_upper,b_lower,periodic):
    # FE_dist = np.linalg.norm(p[n].x - xmin)    # distance from center of particle to minimum energy location
    FE_dist = distance(p[n].x,xmin,dim,b_upper,b_lower,periodic) # for periodic distance
    PE = FE_dist**2
    OE = 0  # overlap energy penalty
    for m in range(len(p)):
        if m != n:
            dist = distance(p[n].x,p[m].x,dim,b_upper,b_lower,periodic)
            if dist < (p[n].r + p[m].r):
                OE = OE + (p[n].r + p[m].r) - dist
    return PE + overlapWeight * OE#PotentialE(FE_dist) + overlapWeight * OE

# Calculates the energy gradient at particle n in p[n]'s location
#   Used to determine the direction of motion for the particle
def dE_dr(n,p,xmin,overlapWeight,dim,b_upper,b_lower,periodic):
    # FE_dist = np.linalg.norm(p[n].x - xmin)  # distance from center of particle to minimum energy location
    # dPE_dr = 2 * (p[n].x - xmin) # derivative of the PE in vector form not as a norm
    dPE_dr = 2 * d_distance(p[n].x,xmin,dim,b_upper,b_lower,periodic) # for periodic gradient??????
    dOE_dr = 0  # overlap energy penalty
    for m in range(len(p)):
        if m != n:
            dist = distance(p[n].x,p[m].x,dim,b_upper,b_lower,periodic)
            if dist < (p[n].r + p[m].r):
                dOE_dr = dOE_dr + (p[m].x - p[n].x) / dist
    return dPE_dr + overlapWeight * dOE_dr

# Pushes overlapping atom away by the same displacement as the last moved atom
# NOT BEING USED- makes everything spaz out
def pusher(n, p, it_perParticle, disp_x, xmin, overlapWeight):
    p_new = p
    stop = False
    i = 0
    m = -1
    while stop == False:
        if i != n:
            dr = p[n].r + p[i].r
            dx = np.linalg.norm(p[n].x - p[i].x)
            if dx < dr:
                stop = True
                m = i
        if i == len(p):
            stop = True
        i = i + 1
    # Push the overlapping particle
    if m != -1:
        p_new[m].x = p[m].x - disp_x
        p_new[n].E = Energy(n,p_new,xmin,overlapWeight) #IAN used p not pnew in (), but i think this is right
        p_new[m].E = Energy(m,p_new,xmin,overlapWeight)
        dE = p_new[n].E + p_new[m].E - p[n].E - p[m].E
        if dE > 0 and it_perParticle>1:
            p_new = pusher(m, p_new, it_perParticle-1, disp_x, xmin, overlapWeight)
            dE = p_new[n].E + p_new[m].E - p[n].E - p[m].E
        if dE > 0:
            p_new = p
    return p_new



#IANS APPROACH:
def MC_Main_Point(n_steps, p, dim, b_lower, b_upper, it_perParticle, disp_max, xmin, periodic, overlapWeight,pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI):
    # Initial Energy Calculation
    for n in range(len(p)):
        p[n].E = Energy(n, p, xmin, overlapWeight,dim,b_upper,b_lower,periodic)
    # Don't show the graph or make an animation
    if showGraph == False and saveAnimation == False:
        for i in range(n_steps):
            print(i)
            for n in range(len(p)):
                x_last = p[n].x
                E_last = p[n].E
                it = 0
                stop = False
                while stop == False:
                    dE = dE_dr(n, p, xmin, overlapWeight,dim,b_upper,b_lower,periodic)
                    disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                    # Set this so that i can manually set a large grain at the minimum point and it wont move
                    if n != 0: #np.array_equal(p[n].x[0],xmin) == False:
                        p[n].x = p[n].x - disp_x
                    # Periodic shift if new position is outside the bounds
                    p[n].x = domainShift(p,n,b_lower,b_upper,periodic)
                    p[n].E = Energy(n, p, xmin, overlapWeight,dim,b_upper,b_lower,periodic)
                    if p[n].E < E_last:
                        stop = True
                    elif it >= it_perParticle:
                        stop = True
                        p[n].x = x_last
                        p[n].E = E_last
                    elif pusherTF == True:
                        p = pusher(n, p, it_perParticle, disp_x, xmin, overlapWeight)
                    it = it + 1
        return p
    # Show the graph but no animation
    elif showGraph == True and saveAnimation == False:
        plt.ion()
        fig = plt.figure()
        # plotCircles(p, b_lower, b_upper, "red")
        plots(p, b_lower, b_upper, dim)
        fig.canvas.draw()
        plt.pause(pltTime)
        for i in range(n_steps):
            for n in range(len(p)):
                x_last = p[n].x
                E_last = p[n].E
                it = 0
                stop = False
                while stop == False:
                    dE = dE_dr(n, p, xmin, overlapWeight,dim,b_upper,b_lower,periodic)
                    disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                    # Set this so that i can manually set a large grain at the minimum point and it wont move
                    if n != 0: #np.array_equal(p[n].x[0],xmin) == False:
                        p[n].x = p[n].x - disp_x
                    # Periodic shift if new position is outside the bounds
                    p[n].x = domainShift(p, n, b_lower, b_upper,periodic)
                    p[n].E = Energy(n, p, xmin, overlapWeight,dim,b_upper,b_lower,periodic)
                    if p[n].E < E_last:
                        stop = True
                    elif it >= it_perParticle:
                        stop = True
                        p[n].x = x_last
                        p[n].E = E_last
                    elif pusherTF == True:
                        p = pusher(n, p, it_perParticle, disp_x, xmin, overlapWeight)
                    it = it + 1
            fig.clear(True)
            # plt.title(i + 1)
            # plotCircles(p, b_lower, b_upper, "red")
            plots(p, b_lower, b_upper, dim)
            plt.title(i + 1)
            fig.canvas.draw()
            plt.pause(pltTime)
        return p
    elif showGraph == True and saveAnimation == True:
        moviewriter = ani.ImageMagickWriter()
        plt.ion()
        fig = plt.figure()
        name = aniName + "." + aniType
        with moviewriter.saving(fig, name, aniDPI):
            # plotCircles(p, b_lower, b_upper, "red")
            plots(p, b_lower, b_upper, dim)
            fig.canvas.draw()
            moviewriter.grab_frame()
            plt.pause(pltTime)
            for i in range(n_steps):
                for n in range(len(p)):
                    x_last = p[n].x
                    E_last = p[n].E
                    it = 0
                    stop = False
                    while stop == False:
                        dE = dE_dr(n, p, xmin, overlapWeight,dim,b_upper,b_lower,periodic)
                        disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                        # Set this so that i can manually set a large grain at the minimum point and it wont move
                        if n != 0: #np.array_equal(p[n].x[0],xmin) == False:
                            p[n].x = p[n].x - disp_x
                        # Periodic shift if new position is outside the bounds
                        p[n].x = domainShift(p, n, b_lower, b_upper,periodic)
                        p[n].E = Energy(n, p, xmin, overlapWeight,dim,b_upper,b_lower,periodic)
                        if p[n].E < E_last:
                            stop = True
                        elif it >= it_perParticle:
                            stop = True
                            p[n].x = x_last
                            p[n].E = E_last
                        elif pusherTF == True:
                            p = pusher(n, p, it_perParticle, disp_x, xmin, overlapWeight)
                        it = it + 1
                fig.clear(True)
                # plt.title(i+1)
                # plotCircles(p, b_lower, b_upper, "red")
                plots(p, b_lower, b_upper, dim)
                plt.title(i + 1)
                fig.canvas.draw()
                moviewriter.grab_frame()
                plt.pause(pltTime)
        moviewriter.finish()
        return p
    elif showGraph == False and saveAnimation == True:
        raise ValueError("Can't save the animation without showing the graph.")



def writeText(fileName,p,dim,header,b_lower,b_upper):
    name = fileName + ".txt"
    file = open(name, "w+")
    if header == True:
        file.write('Generated at: ' + datetime.now().strftime("%m/%d/%Y %H:%M:%S") + '\n')
        file.write('Dimensions: ' + str(dim) + '\n')
        file.write('Number of Particles: ' + str(len(p)) + '\n')
        file.write('Domain: ' + str(b_lower) + ' to ' + str(b_upper) + '\n')
    if dim == 2:
        file.write('%8s %12s %12s %8s\n' % ('x', 'y', 'z', 'r'))
        for n in range(len(p)):
            file.write('%12.6f %12.6f %12.6f %8.4f\n' % (p[n].x[0, 0], p[n].x[0, 1], 0, p[n].r))
    elif dim == 3:
        file.write('%8s %12s %12s %8s\n' % ('x', 'y', 'z', 'r'))
        for n in range(len(p)):
            file.write('%12.6f %12.6f %12.6f %8.4f\n' % (p[n].x[0, 0], p[n].x[0, 1], p[n].x[0, 2], p[n].r))
    else:
        raise ValueError('Dimension Error in writeText()')
    # return file

# def gifToMP4()


# script for volumeFractionSampling
def dist_in_hull(points, n):
    dims = points.shape[-1]
    hull = points[ConvexHull(points).vertices]
    deln = hull[Delaunay(hull).simplices]

    vols = np.abs(det(deln[:, :dims, :] - deln[:, dims:, :])) / np.math.factorial(dims)
    sample = np.random.choice(len(vols), size=n, p=vols / vols.sum())

    return np.einsum('ijk, ij -> ik', deln[sample], dirichlet.rvs([1] * (dims + 1), size=n))


# Samples sampleNum random(?) points in a convex hull of particle centers to calculate volume fraction solid
def volumeFractionSampling(p,dim,b_lower,b_upper,sampleNum,samplePlot):
    x = np.zeros((len(p), dim))
    if dim == 2:
        Circle = namedtuple("Circle", "x y r")
        circles = []
        # Restructure the data
        for n in range(len(p)):
            x[n, 0] = p[n].x[0, 0]
            x[n, 1] = p[n].x[0, 1]
            circles.append(Circle(p[n].x[0, 0], p[n].x[0, 1], p[n].r))
        # Calculate the convex hull and the array of points in the hull
        hull = ConvexHull(x)
        sample = dist_in_hull(x, sampleNum)
        count = 0
        colorCode = np.zeros(len(sample))
        for n in range(len(sample)):
            if any((sample[n, 0] - circle.x) ** 2 + (sample[n, 1] - circle.y) ** 2 <= (circle.r ** 2) for circle in
                   circles):
                count += 1
                colorCode[n] = 1
        if samplePlot == True:
            fig = plt.figure()
            plots(p, b_lower, b_upper, dim)
            for simplex in hull.simplices:
                plt.plot(x[simplex, 0], x[simplex, 1], 'b-')
            plt.scatter(sample[:, 0], sample[:, 1], c=np.where(colorCode == 1, 'red', 'gray'), s=2)
            plt.autoscale()
            # print("Of the ", sampleNum, " points, ", count, " were in particles.")
            # print("-> Approximated volume fraction = ", count / sampleNum)
            return count/sampleNum, fig

    elif dim == 3:
        Circle = namedtuple("Circle", "x y z r")
        circles = []
        # Restructure the data
        for n in range(len(p)):
            x[n, 0] = p[n].x[0, 0]
            x[n, 1] = p[n].x[0, 1]
            x[n, 2] = p[n].x[0, 2]
            circles.append(Circle(p[n].x[0, 0], p[n].x[0, 1], p[n].x[0, 2], p[n].r))

        # Calculate the convex hull and the array of points in the hull
        hull = ConvexHull(x)
        sample = dist_in_hull(x, sampleNum)

        count = 0
        colorCode = np.zeros(len(sample))
        for n in range(len(sample)):
            if any((sample[n, 0] - circle.x) ** 2 + (sample[n, 1] - circle.y) ** 2 + (sample[n, 2] - circle.z) ** 2 <=
                   (circle.r ** 2) for circle in circles):
                count += 1
                colorCode[n] = 1

        if samplePlot == True:
            fig = plt.figure()
            ax = plots(p, b_lower, b_upper, dim)
            for simplex in hull.simplices:
                plt.plot(x[simplex, 0], x[simplex, 1], x[simplex,2], 'b-')
            # Now plot the random points!
            ax.scatter3D(sample[:, 0], sample[:, 1], sample[:, 2], c=np.where(colorCode == 1, 'red', 'gray'),s=2,alpha=0.3)
            plt.autoscale()
            # print("Of the ", sampleNum, " points, ", count, " were in particles.")
            # print("-> Approximated volume fraction = ", count / sampleNum)
            return count / sampleNum, fig

    return count/sampleNum



# Density reduction by particle deletion
def densityReduction(p,dim,b_lower,b_upper,red_steps,periodic,sampleNum,samplePlot,it_perParticle, disp_max, xmin, overlapWeight,pusherTF,pltTime,aniName,aniType,aniDPI,targetVF,vfTol,maxloop,max_part_del): #samplePlot
    vf = volumeFractionSampling(p,dim,b_lower,b_upper,sampleNum,False)
    lp = 0
    part_del = 0
    while lp < maxloop and part_del < max_part_del:
        lp += 1
        if vf > targetVF + vfTol:
            p_new = np.delete(p, randint(0,len(p)-1),axis=0)
            p_new = MC_Main_Drop(red_steps, p, dim, b_lower, b_upper, it_perParticle, disp_max, xmin, periodic,
                                 overlapWeight, 2, pusherTF,False,pltTime,False,aniName,aniType,aniDPI)
            # p_new = MC_Main_Point(red_steps, p_new, dim, b_lower, b_upper, it_perParticle, disp_max, xmin, periodic,
            #                       overlapWeight,pusherTF,False,pltTime,False,aniName,aniType,aniDPI)

            vf_new = volumeFractionSampling(p_new,dim,b_lower,b_upper,sampleNum,samplePlot=False)
            # If everything is done
            if vf_new <= targetVF + vfTol and vf_new >= targetVF - vfTol:
                vf = vf_new
                p = p_new
                # lp = maxloop + 1
                print("Success- Volume Fraction = ",vf)
                print("* Check the graph to make sure it is one cohesive block of particles")
                break
            elif vf_new > targetVF + vfTol:
                p = p_new
                part_del += 1
                vf = vf_new
        #         # lp += 1
        #     # elif vf_new < targetVF - vfTol:
        #     #     lp += 1
        # print("loop = ",lp)
        # print("particles deleted = ",part_del)
    if part_del == max_part_del:
        print("FAILED: Maximum number of particles deleted")
    vf = volumeFractionSampling(p,dim,b_lower,b_upper,sampleNum,samplePlot=True)
    print("Volume Fraction = ",vf[0])
    return p


############################################################################################
############################################################################################
############################################################################################
# PARTICLE DROP
############################################################################################

# Distance for periodic not in the drop axis - hard coded for drop axis = 1
# Distance function to implement periodic boundary conditions if needed
# p1 and p2 will mainly be the particle x cooridnates (p[n].x)
def distanceDrop(p1,p2,dim,b_upper,b_lower,periodic):
    if periodic == False:
        return np.linalg.norm(p1 - p2)
    if periodic == True:
        diff = p1 - p2
        perCheck = np.any(abs(diff) > 0.5*(b_upper-b_lower))
        if perCheck == False:
            return np.linalg.norm(p1 - p2)
        elif perCheck == True:
            if dim == 2:
                pD = b_upper-b_lower # periodic domain lengths
                # make an array of the image cell shifts possible in 2D
                perArray = np.array([[-1,0]*pD,[0,0]*pD,[1,0]*pD])
                # shift the second particle into all possible image cells
                perDist = (list(np.linalg.norm(p1 - (p2 + perArray[n])) for n in range(len(perArray))))
                return min(perDist)
            elif dim == 3:
                pD = b_upper-b_lower # periodic domain lengths
                # make an array of the image cell shifts possible in 3D
                perArray = np.array([[-1,0, -1]*pD,
                                     [0,0, -1]*pD,
                                     [1,0, -1]*pD,
                                     [-1, 0, 0] * pD,
                                     [0, 0, 0] * pD,
                                     [1, 0, 0] * pD,
                                     [-1, 0, 1] * pD,
                                     [0, 0, 1] * pD,
                                     [1, 0, 1] * pD])
                # shift the second particle into all possible image cells
                perDist = (list(np.linalg.norm(p1 - (p2 + perArray[n])) for n in range(len(perArray))))
                return min(perDist)


# Gradient distance for drop (so only on the specified axis
# for the free energy gradient to be periodic
# PERIODIC REWRITE means the drop direction periodicity isnt needed anymore
def d_distanceDrop(p1,p2,dim,b_upper,b_lower,periodic,dropAxis):
    return p1 - p2
    # if periodic == False:
    #     return p1 - p2
    # if periodic == True:
    #     diff = p1 - p2
    #     perCheck = np.any(abs(diff) > 0.5*(b_upper[dropAxis]-b_lower[dropAxis]))
    #     if perCheck == False:
    #         return p1 - p2
    #     elif perCheck == True:
    #         pD = b_upper[dropAxis]-b_lower[dropAxis] # periodic domain lengths
    #         # make an array of the image cell shifts possible in 2D
    #         perArray = np.array([-1*pD,0*pD,1*pD])
    #         # shift the second particle into all possible image cells
    #         perDist = (list(np.linalg.norm(p1 - (p2 + perArray[n])) for n in range(len(perArray))))
    #         return p1 - (p2 + perArray[perDist.index(min(perDist))])


# Calculates the energy of the given particle n in p[n]
def EnergyDrop(n,p,xmin,overlapWeight,dropAxis,dim,b_upper,b_lower,periodic):
    # FE_dist = p[n].x[0,dropAxis] - xmin    # distance from center of particle to minimum energy location
    FE_dist = distanceDrop(p[n].x[0,dropAxis],xmin,dim,b_upper,b_lower,periodic)
    PE = FE_dist**2
    OE = 0  # overlap energy penalty
    for m in range(len(p)):
        if m != n:
            dist = distanceDrop(p[n].x,p[m].x,dim,b_upper,b_lower,periodic)
            if dist < (p[n].r + p[m].r):
                OE = OE + (p[n].r + p[m].r) - dist
    return PE + overlapWeight * OE

# Calculates the energy gradient at particle n in p[n]'s location
#   Used to determine the direction of motion for the particle
def dE_drDrop(n,p,dim,xmin,overlapWeight,dropAxis,b_upper,b_lower,periodic):
    # FE_dist = np.linalg.norm(p[n].x - xmin)  # distance from center of particle to minimum energy location
    if dim == 2:
        dPE_dr = np.array([0,0])
    elif dim == 3:
        dPE_dr = np.array([0, 0, 0])
    # dPE_dr[dropAxis] = 2 * (p[n].x[0, dropAxis] - xmin)  # ** 2 # derivative of the
    dPE_dr[dropAxis] = 2 * d_distanceDrop(p[n].x[0,dropAxis],xmin,dim,b_upper,b_lower,periodic,dropAxis)
    dOE_dr = 0  # overlap energy penalty
    for m in range(len(p)):
        if m != n:
            dist = distanceDrop(p[n].x,p[m].x,dim,b_upper,b_lower,periodic)
            if dist < (p[n].r + p[m].r):
                dOE_dr = dOE_dr + (p[m].x - p[n].x) / dist
    return dPE_dr + overlapWeight * dOE_dr

# Particles drop in the x direction (x=0 is the lowest energy)
def MC_Main_Drop(n_steps, p, dim, b_lower, b_upper, it_perParticle, disp_max, xmin, periodic, overlapWeight, dropAxis, pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI):
    # Initial Energy Calculation
    for n in range(len(p)):
        p[n].E = EnergyDrop(n, p, xmin, overlapWeight,dropAxis,dim,b_upper,b_lower,periodic)
    # Don't show the graph or make an animation
    if showGraph == False and saveAnimation == False:
        for i in range(n_steps):
            print(i)
            for n in range(len(p)):
                x_last = p[n].x
                E_last = p[n].E
                it = 0
                stop = False
                while stop == False:
                    dE = dE_drDrop(n, p, dim, xmin, overlapWeight, dropAxis,b_upper,b_lower,periodic)
                    if np.linalg.norm(dE) == 0.0:
                        disp_x = np.zeros(dim)
                    else:
                        disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                    p[n].x = p[n].x - disp_x
                    # Periodic shift if new position is outside the bounds
                    p[n].x = domainShift(p, n, b_lower, b_upper, periodic)
                    p[n].E = EnergyDrop(n, p, xmin, overlapWeight, dropAxis,dim,b_upper,b_lower,periodic)
                    if p[n].E < E_last:
                        stop = True
                    elif it >= it_perParticle:
                        stop = True
                        p[n].x = x_last
                        p[n].E = E_last
                    elif pusherTF == True:
                        p = pusher(n, p, it_perParticle, disp_x, xmin, overlapWeight)
                    it = it + 1
        return p
    # Show the graph but no animation
    elif showGraph == True and saveAnimation == False:
        plt.ion()
        fig = plt.figure()
        # plotCircles(p, b_lower, b_upper, "red")
        plots(p, b_lower, b_upper, dim)
        fig.canvas.draw()
        plt.pause(pltTime)
        for i in range(n_steps):
            for n in range(len(p)):
                x_last = p[n].x
                E_last = p[n].E
                it = 0
                stop = False
                while stop == False:
                    dE = dE_drDrop(n, p, dim, xmin, overlapWeight, dropAxis,b_upper,b_lower,periodic)
                    if np.linalg.norm(dE) == 0.0:
                        disp_x = np.zeros(dim)
                    else:
                        disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                    p[n].x = p[n].x - disp_x
                    # Periodic shift if new position is outside the bounds
                    p[n].x = domainShift(p, n, b_lower, b_upper, periodic)
                    p[n].E = EnergyDrop(n, p, xmin, overlapWeight, dropAxis,dim,b_upper,b_lower,periodic)
                    if p[n].E < E_last:
                        stop = True
                    elif it >= it_perParticle:
                        stop = True
                        p[n].x = x_last
                        p[n].E = E_last
                    elif pusherTF == True:
                        p = pusher(n, p, it_perParticle, disp_x, xmin, overlapWeight)
                    it = it + 1
            fig.clear(True)
            # plt.title(i+1)
            # plotCircles(p, b_lower, b_upper, "red")
            plots(p, b_lower, b_upper, dim)
            plt.title(i + 1)
            fig.canvas.draw()
            plt.pause(pltTime)
        return p
    elif showGraph == True and saveAnimation == True:
        moviewriter = ani.ImageMagickWriter()
        plt.ion()
        fig = plt.figure()
        name = aniName + "." + aniType
        with moviewriter.saving(fig, name, aniDPI):
            # plotCircles(p, b_lower, b_upper, "red")
            plots(p, b_lower, b_upper, dim)
            fig.canvas.draw()
            moviewriter.grab_frame()
            plt.pause(pltTime)
            for i in range(n_steps):
                for n in range(len(p)):
                    x_last = p[n].x
                    E_last = p[n].E
                    it = 0
                    stop = False
                    while stop == False:
                        dE = dE_drDrop(n, p, dim, xmin, overlapWeight, dropAxis,b_upper,b_lower,periodic)
                        if np.linalg.norm(dE) == 0.0:
                            disp_x = np.zeros(dim)
                        else:
                            disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                        p[n].x = p[n].x - disp_x
                        # Periodic shift if new position is outside the bounds
                        p[n].x = domainShift(p, n, b_lower, b_upper, periodic)
                        p[n].E = EnergyDrop(n, p, xmin, overlapWeight, dropAxis,dim,b_upper,b_lower,periodic)
                        if p[n].E < E_last:
                            stop = True
                        elif it >= it_perParticle:
                            stop = True
                            p[n].x = x_last
                            p[n].E = E_last
                        elif pusherTF == True:
                            p = pusher(n, p, it_perParticle, disp_x, xmin, overlapWeight)
                        it = it + 1
                fig.clear(True)
                # plt.title(i+1)
                # plotCircles(p, b_lower, b_upper, "red")
                plots(p, b_lower, b_upper, dim)
                plt.title(i+1)
                fig.canvas.draw()
                moviewriter.grab_frame()
                plt.pause(pltTime)
        moviewriter.finish()
        return p
    elif showGraph == False and saveAnimation == True:
        raise ValueError("Can't save the animation without showing the graph.")

