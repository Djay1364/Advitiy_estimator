# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 23:33:30 2019

@author: tiwar
"""

import sensor 
import satellite as sat
from constants import *
import matplotlib.pyplot as plt
from numpy.random import multivariate_normal as mvg

sat = satellite.Satellite(v_state0,v_glob_est_state0,v_loc_est_state0,t0)

N = 10000
h = 0.1
gyroVarBias0=np.zeros(3)
sigma_u = 0.0001

sat.setGyroVarBias(gyroVarBias0)

var_bias = np.zeros((N,3))
time = np.zeros(N)
t=0

for i in range(N):
    v_bias_var_k = sat.getGyroVarBias()
    v_bias_var_k1 = v_bias_var_k + sigma_u*(h**0.5)*mvg(GYRO_F_BIAS,np.eye(3))
    var_bias[i,:] = v_bias_var_k1
    t = t+h
    time[i] = t
    sat.setGyroVarBias(v_bias_var_k1)

plt.plot(time,var_bias)