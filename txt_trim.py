from MC_functions import *


# dim = 2                                           # Dimensions
# min_diameter = 8000                                  # Minimum Particle Radius
# # b_upper = [2000000,2500000]#[4000,8000,4000]#100 * np.ones(dim)                      # Upper bound [x,y,z]
# # min_diameter = 5
# max_diameter = 55000
# b_upper = [800000,1250000]
# b_lower = 0 * np.ones(dim)                        # Lower bound [x,y,z]
# periodic = True
# N = 250                                        # Number of particles
# max_ic_its = 120                                  # Maximum number of tries to generate initial condition
# xmin = 0#[1000000,1250000]#0#[0,0,0]#[0, 0]                                   # Lowest energy location [x,y,z], drop uses 1 coordinate
# dropAxis = 1                                      # Axis (0,1,2 = x,y,z) for particles to drop (if using drop)
# # energyType = "Point"                              # Point, Drop
# overlapWeight = 100000000                               # Weight for particle overlap penalty
# n_steps = 120                                     # Number of MC iterations
# it_perParticle = 3                                # Number of iterations to try per particle with Ian's approach
# disp_max = 20000                                      # Maximum particle displacement distance
# pusherTF = False                                  # Whether or not to use the pusher function (W.I.P. atm)
#
# # volumeFractionSampling Details
# sampleNum = 10000                                 # Number of points to sample in convex hull for vol frac
# samplePlot = True                                # Plot the sampled points on the particles
# targetVF = 0.65
# vfTol = 0.05
# red_steps = 3                                      # Number of iterations for MC Main in the particle deletion vol frac
# maxloop = 20
# max_part_del = 20
#
# # Text output details
# txtName = "PeriodicDrop_VF_"
# header = False
#
# # Graphing Details
# showGraph = False
# pltTime = 0.02
# # Animation Details
# saveAnimation = False
# aniType = "gif"
# aniName = "periodicTEST_conv2D_domainProblem"
# aniDPI = 400

dim = 3                                           # Dimensions
min_diameter = 2000#400                                  # Minimum Particle Radius (nm)
# b_upper = [2000000,2500000]#[4000,8000,4000]#100 * np.ones(dim)                      # Upper bound [x,y,z]
# min_diameter = 5
max_diameter = 10000#5000
b_upper = [4000,4000,5000]
b_lower = 0 * np.ones(dim)                        # Lower bound [x,y,z]
periodic = True
N = 50                                        # Number of particles
max_ic_its = 120                                  # Maximum number of tries to generate initial condition
xmin = [6000,6000,6000]#0#[0,0,0]#[0, 0]                                   # Lowest energy location [x,y,z], drop uses 1 coordinate
dropAxis = 2                                      # Axis (0,1,2 = x,y,z) for particles to drop (if using drop)
# energyType = "Point"                              # Point, Drop
overlapWeight = 100000000                               # Weight for particle overlap penalty
n_steps = 200                                     # Number of MC iterations
it_perParticle = 3                                # Number of iterations to try per particle with Ian's approach
disp_max = 10#1000#20000                                      # Maximum particle displacement distance
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
txtName = "PeriodicDrop_VF_"
header = False

# Graphing Details
showGraph = True
pltTime = 0.02
# Animation Details
saveAnimation = True
aniType = "gif"
aniName = "3D_converge_corner_cont"
aniDPI = 400






def txtReconstruct(fileName):
    name = fileName + ".txt"
    p = []
    count = 0
    with open(name, "r") as file:
        next(file)
        for n in file:
            p.append(particle())
            p[count].x = np.zeros((1, 3))
            # print(n.split()[0])
            # print(p[count].x[0])
            p[count].x[0, 0] = float(n.split()[0])
            p[count].x[0, 1] = float(n.split()[1])
            p[count].x[0, 2] = float(n.split()[2])
            p[count].r = float(n.split()[3])
            count = count+1
    return p


p = txtReconstruct("converge_1Large_corner")

# p[n].x = np.zeros((1, dim))

# plotCircles(p,b_lower,b_upper)

plotSpheres(p,b_lower,b_upper)#[300000,300000,150000])#b_upper
plt.show()

# vf = volumeFractionSampling(p,dim,b_lower,b_upper,sampleNum,samplePlot=True)
# plt.show()
# print("Volume Fraction = ",vf[0])
# quit()
# Run more steps
# p = MC_Main_Drop(10, p, dim, b_lower, b_upper, it_perParticle, disp_max, xmin, periodic, overlapWeight, dropAxis,
#                  pusherTF,showGraph,pltTime,saveAnimation,aniName,aniType,aniDPI)
p = MC_Main_Point(200, p, dim, b_lower, b_upper, it_perParticle, disp_max, xmin, False, overlapWeight,
                  pusherTF,showGraph,pltTime,saveAnimation,"converge_Allcorner_cont2",aniType,aniDPI)

writeText("converge_AllLarge_corner_cont2",p,dim,header,b_lower,b_upper)
quit()
p = MC_Main_Drop(15, p, dim, b_lower, b_upper, it_perParticle, disp_max, 0, periodic, overlapWeight, dropAxis,
                 pusherTF,showGraph,pltTime,saveAnimation,"drop",aniType,aniDPI)
writeText("PeriodicDrop_predelete_temp2",p,dim,header,b_lower,b_upper)
# p_2 = densityReduction(p,dim,b_lower,b_upper,red_steps,periodic,sampleNum,samplePlot,it_perParticle, disp_max, xmin,
#                      overlapWeight,pusherTF,pltTime,aniName,aniType,aniDPI,targetVF,vfTol,maxloop,max_part_del)
#
# print("reduced")
# writeText("PeriodicDrop_VF_temp",p,dim,header,b_lower,b_upper)

plt.show()
