import random
import numpy as np
import matplotlib.pyplot as plt
import math
from datetime import datetime

part_num = 200
text_write = False
# fileName = "particle_dist"
# Small particle curve
avg1 = 700 #0.7 #um
stdev1 = 0.5
volper1 = 0.593
# Larger Particle Curve
avg2 = 5890 #5.89 #um
stdev2 = 0.5
volper2 = 0.407

min_d = 100 #0.1 #um
max_d = 5000 #5.0 #um
manual_bins = np.logspace(-1,3,200)

def weight_percent(d):
    radii = np.asarray(d) / 2#/ 1e3 #nm to um
    vol = 4 * math.pi * (radii ** 3) / 3
    tot_weight = np.sum(vol)
    wt_per = (vol/tot_weight) * 100
    return wt_per

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
    part_d = np.where(np.random.rand(number_particles) < small_percent, np.random.lognormal(mu(avg1, stdev1), stdev1, number_particles),
                      np.random.lognormal(mu(avg2, stdev2), stdev2, number_particles))
    if min_d != 0.0:
        part_d = np.where(part_d < min_d, np.full_like(part_d, min_d, dtype=np.double),
                          np.where(part_d > max_d, np.full_like(part_d, max_d, dtype=np.double), part_d))
    return part_d

# d1 = np.random.lognormal(mu(0.7,stdev1), stdev1,10000)
# d2 = np.random.lognormal(mu(5.89,stdev2), stdev2,10000)

pdf1 = pdf(manual_bins,mu(avg1,stdev1),stdev1)
pdf2 = pdf(manual_bins,mu(avg2,stdev2),stdev2)

particle_sizes = particle_size_dist(part_num,min_d,max_d,avg1,stdev1,volper1,avg2,stdev2,volper2)

plt.figure(1)
# plt.hist(d1,bins=manual_bins,density=True)
# plt.hist(d2,bins=manual_bins,density=True)
plt.hist(particle_sizes,bins=manual_bins,density=True)
plt.plot(manual_bins,pdf1)
plt.plot(manual_bins,pdf2)
plt.semilogx()
plt.xlabel('Particle Size (um)')
plt.ylabel('Count')

plt.figure(2)
# particle_sizes = particle_size_dist(1000000,0,0,avg1,stdev1,volper1,avg2,stdev2,volper2)
manual_bins = np.logspace(np.log10(0.2),np.log10(60),200)
plt.hist(particle_sizes,bins=manual_bins,weights=weight_percent(particle_sizes),cumulative=True)
plt.semilogx()
plt.xlim([0.2, 60])
x = [0.2, 0.4, 0.6, 1, 2, 4, 6, 10, 20, 40, 60]
plt.xticks(x, x)
plt.gca().invert_xaxis()

plt.xlabel('Particle Size (um)')
plt.ylabel('Cumulative Weight (%)')

if text_write == True:
    np.save('particle_dist.txt', particle_sizes)

plt.show()