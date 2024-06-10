from __future__ import division
import os
import random
import math
import numpy as np
from numpy import zeros, sqrt, log, where, pi, mean, arange, sin, cos
import matplotlib.pyplot as plt
import csv
from gauss import random, gaussian

# Set initial positions
x = -2.5
y = 0.0

# Set parameters
t = 0.0
dt = 0.005
k = 1.0         # for k_y
kBT = 0.5915
diffusion = 1.0
nstep = 100000
dx = 0.01
dy = 0.01
xmin = -5.0
xmax = 5.0
ymin = -5.0
ymax = 5.0

# Plot the 2D Potential Energy
def energy_2D():
    global x, y
    energy_x = []
    energy_y = []
    energy_y2 = []    
    energy_y3 = []
    range_x = np.linspace(xmin, xmax, 101)
    range_y = np.linspace(ymin, ymax, 101) 
    Ex = 0.1 * x ** 4 - x ** 2
    Ey = 0.5 * k * y ** 2
    for x in range_x:
        Ex = 0.1 * x ** 4 - x ** 2
        energy_x.append(Ex)
    for y in range_y:
        Ey = 0.5 * 0.5 * y ** 2
        energy_y.append(Ey)
        Ey2 = 0.1 * y ** 2
        energy_y2.append(Ey2)
        Ey3 = 0.075 * y ** 2
        energy_y3.append(Ey3)
    
    X, Y = np.meshgrid(range_x, range_y)
    Z = (0.1 * X ** 4 - X ** 2) + (0.5 * k * Y ** 2)    
    
    # plt.figure()
    # CP = plt.contour(X, Y, Z, linewidths=2) 
    # plt.title('Contour Plot')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # # plt.savefig('potential.png', dpi=1080)
    # plt.figure()
    # plt.plot(range_x, energy_x, linewidth=2)
    # plt.xlabel('x positions')
    # plt.ylabel('Energy')
    # # plt.savefig('potential_x.png', dpi=1080)
    # plt.figure()
    # plt.plot(range_y, energy_y, linewidth=2, label='k = 0.5')
    # plt.plot(range_y, energy_y2, linewidth=2, label='k = 2')
    # plt.plot(range_y, energy_y3, linewidth=2, label='k = 0.15')
    # plt.xlabel('y positions')
    # plt.ylabel('Energy')
    # plt.legend()
    # # plt.savefig('potential_y.png', dpi=1080)
    # plt.show()

    # Plot positions vs time
    t, x, y = np.loadtxt('trajectories/xy_traj_100k.trj', usecols=(0,1,2), unpack=True)
    t2, x2, y2 = np.loadtxt('trajectories/xy_traj_100k_D4.trj', usecols=(0,1,2), unpack=True)
    plt.figure()
    plt.plot(t, x, 'r')
    plt.figure()
    plt.plot(t2, x2, 'b')
    # plt.plot(t, y)
    plt.show()

# Run Brownian Dynamics (Langevin)
def bd_2D():
    global x, y, t
    xTraj = np.array([x])
    yTraj = np.array([y])
    time = np.array([t])
    for i in np.linspace(0, dt*nstep, nstep):
        Fx = - 0.1 * 4 * x ** 3 + 2 * x
        Fy = - k * y
        fact1 = diffusion * dt / kBT
        delx = Fx * fact1
        dely = Fy * fact1
        fact2 = sqrt(2 * diffusion * dt)
        x = x + delx + fact2 * gaussian()
        y = y + dely + fact2 * gaussian()
        t = t + dt
        xTraj = np.append(xTraj, x)
        yTraj = np.append(yTraj, y)
        time = np.append(time, t)
    # xyTraj = zip(time, xTraj, yTraj)
    # np.savetxt('trajectories/xy_traj_100k.trj', np.c_[time, xTraj, yTraj])

energy_2D()
# bd_2D()
