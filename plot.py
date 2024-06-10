from __future__ import division
import os
import math
import numpy as np
from numpy import zeros, sqrt, square, log, where, pi, mean, arange, sin, cos
import matplotlib.pyplot as plt
import csv
from scipy import linalg

# Extract Correlations: Normalized Covariances
def norm():
    t, x, y = np.loadtxt('trajectories/xy_traj_100k.trj', usecols=(0,1,2), unpack=True)
    images = 1
    nstep_min = 1
    nstep_max = 100000               # last timestep = nstep_max + 1
    N = nstep_max - nstep_min + 1    # N number of data points, N = nstep_max - nstep_min
    t = t[nstep_min:nstep_max+1]
    x = x[nstep_min:nstep_max+1]
    y = y[nstep_min:nstep_max+1]

    # Compute Mean-Free Coordinates
    x_av = np.sum(x) / N
    y_av = np.sum(y) / N
    X_mf = np.array([])
    Y_mf = np.array([])

    n = 0
    while n < N:
        X_mf = np.append(X_mf, x[n] - x_av)
        Y_mf = np.append(Y_mf, y[n] - y_av)   
        n += 1

    C_x0 = sqrt(np.sum(square(X_mf)) / N)
    C_y0 = sqrt(np.sum(square(Y_mf)) / N)
    #print square(C_x0), (C_x0 * C_y0), square(C_y0)
    #100k: 4.33951138129 1.62416117046 0.60787938453
    #10k: 0.365836997432 0.446187021232 0.544184593995

# Plot Correlations
def plot_corr():
    time, C_xx, C_xy, C_yx, C_yy = np.loadtxt('covariances/cov_N100k_tau3918.txt', usecols=(0,1,2,3,4), unpack=True)

    tau = time / 0.005 
    #C_xx = C_xx / 0.3658369974329        
    #C_xy = C_xy / 0.446187021232    #10k: 0.365836997432 0.446187021232 0.544184593995
    #C_yx = C_yx / 0.446187021232
    #C_yy = C_yy / 0.544184593995

    C_xx = C_xx / 4.33951138129     #100k: 4.33951138129 1.62416117046 0.60787938453   
    C_xy = C_xy / 1.62416117046     
    C_yx = C_yx / 1.62416117046
    C_yy = C_yy / 0.60787938453

    # Plot Covariance vs Lag Time, tau (time/0.005)
    plt.figure()
    plt.plot(tau, C_xx, color='b', linewidth=2.0, label=r'<X(0)X($\tau$)>')
    plt.plot(tau, C_xy, color='k', linewidth=2.0, label=r'<X(0)Y($\tau$)>')
    plt.plot(tau, C_yx, color='g', linewidth=2.0, label=r'<Y(0)X($\tau$)>')
    plt.plot(tau, C_yy, color='r', linewidth=2.0, label=r'<Y(0)Y($\tau$)>')

    plt.title('Correlation vs Lag Time')
    plt.xlabel('Lag Time, ' r'$\tau$')
    plt.ylabel('Correlation')
    plt.ylim((-0.2,1.5))
    plt.xlim((0.0,300))
    plt.legend(loc='upper right')
    
    #plt.savefig('tica_N100k_tau900_all_new.png', dpi=1080)
    # plt.savefig('temp.png', dpi=1080)
    # plt.show()

# norm()
plot_corr()

