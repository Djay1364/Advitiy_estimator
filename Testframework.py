# -*- coding: utf-8 -*-
"""
Created on Thu May  9 16:46:35 2019
@author: Sanskriti
"""

from constants_1U import *
import numpy as np
import satellite
import MEKFmaincode
import testfunction
import sensor
import matplotlib.pyplot as plt
import qnv
#import Q_BO_trajectory

control_step=1
model_step=1
sigma_u=0.01
t0=0

v_q0_BO = np.array([0,0,0,1])
v_w0_BOB = np.array([0,0,0])
v_est_q0_BO = np.array([0,1/2**0.5,0,1/2**0.5])
v_est_w0_BOB = np.array([0.1,0.1,0.1])
v_est_theta0 = np.array([0,0,0])
v_est_bias0 = np.array([0.001,0.001,0.001])
v_state0 = np.hstack((v_q0_BO,v_w0_BOB))
P_k=0.01*np.array([[0.1,0,0,0,0,0],
              [0,0.1,0,0,0,0],
              [0,0,0.1,0,0,0],
              [0,0,0,0.1,0,0],
              [0,0,0,0,0.1,0],
              [0,0,0,0,0,0.1]]) #initial state covariance matrix

v_glob_est_state0 = np.hstack((v_est_q0_BO,v_est_w0_BOB))
v_loc_est_state0 = np.hstack((v_est_theta0,v_est_bias0))

Advitiy = satellite.Satellite(v_state0,v_glob_est_state0,v_loc_est_state0,t0) 
Advitiy.setMag_i(np.array([1,1,1]))
Advitiy.setMag_b_m_c(np.array([1,1,1]))
Advitiy.setPos(np.array([1,1,1]))
Advitiy.setVel(np.array([1,2,1]))
Advitiy.setGyroEstBias(np.array([0.0005,0.0005,0.0005]))
Advitiy.setErrCovariance(P_k)
        
#Make satellite object
init=1
end=1000
k=0
t0=0
l=np.linspace(0,end,end)
m_sgp_output = np.genfromtxt('sgp_output.csv', delimiter=",")
m_magnetic_field_i = np.genfromtxt('mag_output_i.csv',delimiter=",") 
#Advitiy = satellite.Satellite(v_state0,t0) 
position=np.zeros((end,3))
velocity=np.zeros((end,3))

w_est=velocity.copy()
w_gyro=velocity.copy()
w_true=velocity.copy()

w_bias = velocity.copy()
w_est_bias = velocity.copy()

q=np.zeros((end,4))
q_true=np.zeros((end,4))
q_true_est = np.zeros((end,4))
RMSE=np.zeros((end,3))
p=np.zeros((end,6,6))
x=np.zeros((end,6))
xmod=np.zeros((end))
gyroVarBias0=np.zeros(3)
Advitiy.setGyroVarBias(gyroVarBias0)
q_true_bo = np.array([0,0,0,1])

for i in range(end):
    print(i+1)
    
    w_true_bob = np.array([2*np.sin(0.1*i), 2, 2])
    #w_true_bob = np.array([20*sin(0.01*i),20*sin(0.01*i),20*sin(0.01*i)])
    Advitiy.setState(np.hstack((q_true_bo,w_true_bob)))

    q_true_bo = testfunction.propogate_quaternion(w_true_bob,q_true_bo)
    q_true_bo = q_true_bo/np.linalg.norm(q_true_bo)
    
    w_true[i,:] = w_true_bob 
    q_true[i,:] = q_true_bo    
    
    Advitiy.setW_BO_b(w_true_bob)
    Advitiy.setQ_BO(q_true_bo)
    
    #N_u=np.random.multivariate_normal(np.array([0,0,0]),01e-4*np.eye(3))
    #gyroVarBias[i,:] = gyroVarBias[i,:] + sigma_u*((0.1*model_step)**0.5)*N_u#
    #Advitiy.setGyroVarBias(gyroVarBias[i,:])

    Advitiy.setMag_i(m_magnetic_field_i[i,1:4]*1e-9)
    mag,v_n_m = sensor.magnetometer(Advitiy)
    Advitiy.setMag_b_m_c(mag)
    Advitiy.setPos(m_sgp_output[i,1:4]*1e-6)
    Advitiy.setVel(m_sgp_output[i,4:7]*1e-6)
    #print(m_sgp_output[i,:])
    
    q[i,:], p[i,:], x[i, :] = MEKFmaincode.estimator(Advitiy)
    
    w_bias[i,:] = Advitiy.getGyroVarBias()
    w_est_bias[i,:] = Advitiy.getGyroEstBias()
    #print(w_est_bias[i,:])
    #RMSE[i,:]=np.sqrt(((b-f[3:6])**2)/end)
    q_kk=q[i,:]
    P_k=p[i,:]
    w_m_k=sensor.gyroscope(Advitiy)
    w_gyro[i,:]=testfunction.w_bob_calc(w_m_k,q_kk,-v_w_IO_o,np.array([0,0,0]))#w_est_bias[i,:])#w_m_k
    w_est[i,:]=w_gyro[i,:]-w_est_bias[i,:]
    xmod[i]=np.linalg.norm(x[i, :]) 

#plt.plot(l,p[:,0,0], 'r')
#plt.xlabel('number of iterations')  
#plt.ylabel('quaternion')  
#plt.plot(l,w_gyro[:,0], 'b')
#plt.plot(l,w_true[:,0], 'r')
#plt.plot(l,w_est[:,0], 'g')
#plt.plot(l,gyroVarBias[:,1], 'r')
#plt.plot(l,w_bias[:,1], 'y')
#plt.plot(l,w_est_bias[:,1], 'r')

plt.plot(l,q_true[:,1], 'b')
plt.plot(l,q[:,1], 'r')
#plt.plot(l,gyroVarBias[:,2])'''