from __future__ import division
import os
import math
import numpy as np
from numpy import zeros, sqrt, square, log, where, pi, mean, average, arange, sin, cos
import matplotlib.pyplot as plt
import csv
from scipy import linalg


# Set Parameters
images = 1
tau_min = 8                     # iterates from tau = tau_min to tau_max-1.
tau_max = 1400
nstep_min = 0
nstep_max = 10                  # last timestep = nstep_max + 1
# N = nstep_max - nstep_min + 1   # N number of data points, N = nstep_max - nstep_min

# Load Trajectory into Arrays
time, x, y = np.loadtxt('trajectories/xy_traj_100k.trj', usecols=(0,1,2), unpack=True)
time = time[nstep_min:nstep_max+1]
x = x[nstep_min:nstep_max+1]
y = y[nstep_min:nstep_max+1]
N = len(time)                   # N number of data points, N = nstep_max - nstep_min

def moo():
    # global x, y, time
    X_mf = np.array([])
    Y_mf = np.array([])

    # Compute Average of x, y Trajectories
    x_av = np.sum(x) / N
    y_av = np.sum(y) / N

    # Compute Mean-Free Coordinates
    n = 0
    while n < N:
        X_mf = np.append(X_mf, x[n] - x_av)
        Y_mf = np.append(Y_mf, y[n] - y_av)   
        n += 1

    # Calculating the Covariances
    C_xx = np.array([])
    C_xy = np.array([])
    C_yx = np.array([])
    C_yy = np.array([])

    for t in xrange(tau_min,tau_max+1):     # Loop over different lag times
        C_a = np.array([])
        C_b = np.array([])
        C_c = np.array([])
        C_d = np.array([])
        if t > N:
            raise RuntimeError ("Invalid lag time. Lower tau")
        for i in xrange(N-t):        # Loop over the trajectory
            C_a = np.append(C_a, X_mf[i]*X_mf[i+t])
            C_b = np.append(C_b, X_mf[i]*Y_mf[i+t])
            C_c = np.append(C_c, Y_mf[i]*X_mf[i+t])
            C_d = np.append(C_d, Y_mf[i]*Y_mf[i+t])
        C_xx = np.append(C_xx, np.sum(C_a) / (N-t)) 
        C_xy = np.append(C_xy, np.sum(C_b) / (N-t))
        C_yx = np.append(C_yx, np.sum(C_c) / (N-t))
        C_yy = np.append(C_yy, np.sum(C_d) / (N-t))

        cov = zip(time[tau_min:tau_max+1], C_xx, C_xy, C_yx, C_yy)
        with open('temp.txt', 'w') as f:
            np.savetxt(f, cov)
            f.flush()
            os.fsync(f.fileno())
    
    # # Compute Variance (for Normalization)
    # C_x0 = sqrt(np.sum(square(X_mf)) / N)
    # C_y0 = sqrt(np.sum(square(Y_mf)) / N)

    # # Compute Correlations: Normalized Covariances
    # C_xx = C_xx / square(C_x0)
    # C_xy = C_xy / (C_x0 * C_y0)
    # C_yx = C_yx / (C_y0 * C_x0)
    # C_yy = C_yy / square(C_y0)
    # return C_xx, C_xy, C_yx, C_yy, time

# moo()
