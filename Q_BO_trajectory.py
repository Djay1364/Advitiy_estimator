# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 01:58:52 2019

@author: tiwar
"""
import numpy as np
import matplotlib.pyplot as plt
from testfunction import propogate_quaternion
##Assuming W_bob = constant for the first case
N = 1000
pi = np.pi
v_w_bob = np.array([0.1,0.1,0.1])
h = 0.1
t = 0
q_BO = np.array([0,0,0,1])

quat = np.zeros([N,4])
time =[]
v_q_bob = np.zeros([N,4])
th = np.zeros([N,3])

for i in range(N):
    t=t+h
    #th = h*(i+1)*v_w_bob[1]
    #v_q_bob[i,:] = np.array([0,0,np.sin(th/2),np.cos(th/2)])
    #v_w_bob = np.array([np.sin(t),np.sin(t),np.sin(t)])
    q_BO = propogate_quaternion(v_w_bob,q_BO)
    q_BO = q_BO/np.sqrt(q_BO[0]**2+q_BO[1]**2+q_BO[2]**2+q_BO[3]**2)
    #print(q)
    print(np.sqrt(q_BO[0]**2+q_BO[1]**2+q_BO[2]**2+q_BO[3]**2))
    quat[i,:] = q_BO
    time.append(t)
    
#plt.plot(time,v_q_bob[:,3])
plt.plot(time,quat[:,3])
    

