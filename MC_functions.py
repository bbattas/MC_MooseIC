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
import math



class particle:
    def __init__(self):
        self.r = 0
        self.x = []
        self.E = 0

def gauss(x,mean,stdev):
    return np.exp(- (x - mean)**2 / (2*stdev**2))
# Particle size distribution function, x = random number 0-1
def part_size_dist(min_radius):
    x = random.random()
    out = 300 * gauss(x,0.5,0.15)   #OUTPUT VALUE
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
                            print("ERROR: MAX ITERATIONS REACHED - there will be particle overlap")
                    else:
                        loop = 0
            else:
                loop = 0
    print("Iterations used for initial condition = ", it)
    print(" ")
    return p


# Circle Plotting Function
def plotCircles(p,b_lower,b_upper,col):
    ax = plt.gca()
    ax.set_xlim((b_lower[0], b_upper[0]))
    ax.set_ylim((b_lower[1], b_upper[1]))
    ax.set_aspect('equal')
    for n in range(len(p)):
        circle = plt.Circle(*p[n].x, p[n].r, fill=False)#,color=col
        ax.add_artist(circle)
    return ax

# Potential Energy Function
def PotentialE(distFromIdeal):
    # distFromIdeal is the norm between the particle and the lowest energy position
    return distFromIdeal**2

# Calculates the energy of the given particle n in p[n]
def Energy(n,p,xmin,overlapWeight):
    FE_dist = np.linalg.norm(p[n].x - xmin)    # distance from center of particle to minimum energy location
    PE = FE_dist**2
    OE = 0  # overlap energy penalty
    for m in range(len(p)):
        if m != n:
            dist = np.linalg.norm(p[n].x - p[m].x)
            if dist < (p[n].r + p[m].r):
                OE = OE + (p[n].r + p[m].r) - dist
    return PE + overlapWeight * OE#PotentialE(FE_dist) + overlapWeight * OE

# Calculates the energy gradient at particle n in p[n]'s location
#   Used to determine the direction of motion for the particle
def dE_dr(n,p,xmin,overlapWeight):
    # FE_dist = np.linalg.norm(p[n].x - xmin)  # distance from center of particle to minimum energy location
    dPE_dr = 2 * (p[n].x - xmin) #** 2 # derivative of the
    dOE_dr = 0  # overlap energy penalty
    for m in range(len(p)):
        if m != n:
            dist = np.linalg.norm(p[n].x - p[m].x)
            if dist < (p[n].r + p[m].r):
                dOE_dr = dOE_dr + (p[m].x - p[n].x) / dist
    return dPE_dr + overlapWeight * dOE_dr

# Pushes overlapping atom away by the same displacement as the last moved atom
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
def MC_Main_Point(n_steps, p, b_lower, b_upper, it_perParticle, disp_max, xmin, overlapWeight,pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI):
    # Initial Energy Calculation
    for n in range(len(p)):
        p[n].E = Energy(n, p, xmin, overlapWeight)
    # Don't show the graph or make an animation
    if showGraph == False and saveAnimation == False:
        for i in range(n_steps):
            for n in range(len(p)):
                x_last = p[n].x
                E_last = p[n].E
                it = 0
                stop = False
                while stop == False:
                    dE = dE_dr(n, p, xmin, overlapWeight)
                    disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                    p[n].x = p[n].x - disp_x
                    p[n].E = Energy(n, p, xmin, overlapWeight)
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
        plotCircles(p, b_lower, b_upper, "red")
        fig.canvas.draw()
        plt.pause(pltTime)
        for i in range(n_steps):
            for n in range(len(p)):
                x_last = p[n].x
                E_last = p[n].E
                it = 0
                stop = False
                while stop == False:
                    dE = dE_dr(n, p, xmin, overlapWeight)
                    disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                    p[n].x = p[n].x - disp_x
                    p[n].E = Energy(n, p, xmin, overlapWeight)
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
            plt.title(i+1)
            plotCircles(p, b_lower, b_upper, "red")
            fig.canvas.draw()
            plt.pause(pltTime)
        return p
    elif showGraph == True and saveAnimation == True:
        moviewriter = ani.ImageMagickWriter()
        plt.ion()
        fig = plt.figure()
        name = aniName + "." + aniType
        with moviewriter.saving(fig, name, aniDPI):
            plotCircles(p, b_lower, b_upper, "red")
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
                        dE = dE_dr(n, p, xmin, overlapWeight)
                        disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                        p[n].x = p[n].x - disp_x
                        p[n].E = Energy(n, p, xmin, overlapWeight)
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
                plt.title(i+1)
                plotCircles(p, b_lower, b_upper, "red")
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
        if any((sample[n, 0] - circle.x) ** 2 + (sample[n, 1] - circle.y) ** 2 <= (circle.r ** 2) for circle in circles):
            count += 1
            colorCode[n] = 1

    if samplePlot == True:
        fig = plt.figure()
        plotCircles(p, b_lower, b_upper, 'red')
        for simplex in hull.simplices:
            plt.plot(x[simplex, 0], x[simplex, 1], 'b-')
        plt.scatter(sample[:, 0], sample[:, 1], c=np.where(colorCode == 1, 'red', 'gray'), s=2)
        plt.autoscale()
        # print("Of the ", sampleNum, " points, ", count, " were in particles.")
        # print("-> Approximated volume fraction = ", count / sampleNum)
        return count/sampleNum, fig

    # print("Of the ", sampleNum, " points, ", count, " were in particles.")
    # print("-> Approximated volume fraction = ", count / sampleNum)
    return count/sampleNum


############################################################################################
############################################################################################
############################################################################################
# PARTICLE DROP
############################################################################################

# Calculates the energy of the given particle n in p[n]
def EnergyDrop(n,p,xmin,overlapWeight):
    FE_dist = p[n].x[0,1] - xmin    # distance from center of particle to minimum energy location
    PE = FE_dist**2
    OE = 0  # overlap energy penalty
    for m in range(len(p)):
        if m != n:
            dist = np.linalg.norm(p[n].x - p[m].x)
            if dist < (p[n].r + p[m].r):
                OE = OE + (p[n].r + p[m].r) - dist
    return PE + overlapWeight * OE

# Calculates the energy gradient at particle n in p[n]'s location
#   Used to determine the direction of motion for the particle
def dE_drDrop(n,p,xmin,overlapWeight):
    # FE_dist = np.linalg.norm(p[n].x - xmin)  # distance from center of particle to minimum energy location
    dPE_dr = np.array([0,2 * (p[n].x[0,1] - xmin)]) #** 2 # derivative of the
    dOE_dr = 0  # overlap energy penalty
    for m in range(len(p)):
        if m != n:
            dist = np.linalg.norm(p[n].x - p[m].x)
            if dist < (p[n].r + p[m].r):
                dOE_dr = dOE_dr + (p[m].x - p[n].x) / dist
    return dPE_dr + overlapWeight * dOE_dr

# Particles drop in the x direction (x=0 is the lowest energy)
def MC_Main_Drop(n_steps, p, b_lower, b_upper, it_perParticle, disp_max, xmin, overlapWeight,pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI):
    # Initial Energy Calculation
    for n in range(len(p)):
        p[n].E = EnergyDrop(n, p, xmin, overlapWeight)
    # Don't show the graph or make an animation
    if showGraph == False and saveAnimation == False:
        for i in range(n_steps):
            for n in range(len(p)):
                x_last = p[n].x
                E_last = p[n].E
                it = 0
                stop = False
                while stop == False:
                    dE = dE_drDrop(n, p, xmin, overlapWeight)
                    disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                    p[n].x = p[n].x - disp_x
                    p[n].E = EnergyDrop(n, p, xmin, overlapWeight)
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
        plotCircles(p, b_lower, b_upper, "red")
        fig.canvas.draw()
        plt.pause(pltTime)
        for i in range(n_steps):
            for n in range(len(p)):
                x_last = p[n].x
                E_last = p[n].E
                it = 0
                stop = False
                while stop == False:
                    dE = dE_drDrop(n, p, xmin, overlapWeight)
                    disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                    p[n].x = p[n].x - disp_x
                    p[n].E = EnergyDrop(n, p, xmin, overlapWeight)
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
            plt.title(i+1)
            plotCircles(p, b_lower, b_upper, "red")
            fig.canvas.draw()
            plt.pause(pltTime)
        return p
    elif showGraph == True and saveAnimation == True:
        moviewriter = ani.ImageMagickWriter()
        plt.ion()
        fig = plt.figure()
        name = aniName + "." + aniType
        with moviewriter.saving(fig, name, aniDPI):
            plotCircles(p, b_lower, b_upper, "red")
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
                        dE = dE_drDrop(n, p, xmin, overlapWeight)
                        disp_x = random.random() * disp_max * dE / np.linalg.norm(dE)
                        p[n].x = p[n].x - disp_x
                        p[n].E = EnergyDrop(n, p, xmin, overlapWeight)
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
                plt.title(i+1)
                plotCircles(p, b_lower, b_upper, "red")
                fig.canvas.draw()
                moviewriter.grab_frame()
                plt.pause(pltTime)
        moviewriter.finish()
        return p
    elif showGraph == False and saveAnimation == True:
        raise ValueError("Can't save the animation without showing the graph.")

