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
import Q_BO_trajectory

control_step=1
model_step=1
t0=0
v_q0_BO = np.array([0,0,0,1])
v_w0_BOB = np.array([0,0,0])
v_est_q0_BO = np.array([0,1/2**0.5,0,1/2**0.5])
v_est_w0_BOB = np.array([0.1,0.1,0.1])
v_state0 = np.hstack((v_q0_BO,v_w0_BOB))
v_est_state0 = np.hstack((v_est_q0_BO,v_est_w0_BOB))
Advitiy = satellite.Satellite(v_state0,v_est_state0,t0) 
Advitiy.setMag_i(np.array([1,1,1]))
Advitiy.setMag_b_m_c(np.array([1,1,1]))
Advitiy.setPos(np.array([1,1,1]))
Advitiy.setVel(np.array([1,2,1]))
b=np.array([0.0001,0.0001,0.0001])
b_e=np.array([0,0,0])

#print("v_state0")
#print(v_state0)                         

#Make satellite object
init=1
end=200
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
q=np.zeros((end,4))
q_true=np.zeros((end,4))
q_true_bo = np.array([0,0,0,1])
q_true_est = np.zeros((end,4))
RMSE=np.zeros((end,3))
p=np.zeros((end,6,6))
x=np.zeros((end,6))
xmod=np.zeros((end))

for i in range(end):
    
    w_true_bob = np.array([0.2,0.2,0.2])
    w_true[i,:] = w_true_bob 
    q_true_bo = testfunction.propogate_quaternion(w_true_bob,q_true_bo)
    q_true_bo = q_true_bo/np.sqrt(q_true_bo[0]**2+q_true_bo[1]**2+q_true_bo[2]**2+q_true_bo[3]**2)
    q_true[i,:] = q_true_bo    
    Advitiy.setW_BO_b(w_true_bob)
    Advitiy.setQ_BO(q_true_bo)
    #Advitiy.setGyroVarBias(np.random.normal(0.001,0.1,3))
    Advitiy.setGyroVarBias(np.array([0.001,0.001,0.001]))
    w_m_k=sensor.gyroscope(Advitiy)
    #w_t_k=satellite.getW_BO_b(Advitiy)
    #w_t_bib=satellite.getW_BI_b(Advitiy)
    Advitiy.setMag_i(m_magnetic_field_i[i,1:4]*1e-9)
    mag = sensor.magnetometer(Advitiy)
    Advitiy.setMag_b_m_c(mag)
    Advitiy.setPos(m_sgp_output[i,1:4]*1e-6)
    Advitiy.setVel(m_sgp_output[i,4:7]*1e-6)
    position[i]=Advitiy.getPos()
    velocity[i]=Advitiy.getVel()
    #print(i)
    b_e=b-x[i-1,3:6]
    #print(b_e)
    #print(b_e)
    q[i,:], p[i :], f=MEKFmaincode.estimator(Advitiy)
    x[i, :]=f
    #print(f[0:3])
    RMSE[i,:]=np.sqrt(((b-f[3:6])**2)/end)
    q_kk=q[i,:]
    delta_x=[0,0,0,f[3:6]]
    delta_q = np.hstack([0.5*f[0:3],1])
    q_true_est[i,:] = qnv.quatMultiplyNorm(delta_q,q_kk)
    #print(i)
    #print(f[3:6])
    P_k=p[i,:]
    #w_true[i,:]=w_t_bib
    w_gyro[i,:]=testfunction.w_bob_calc(w_m_k,q_true_bo,-v_w_IO_o,b_e)#w_m_k
    w_est[i,:]=w_gyro[i,:]+x[i,3:6]-b
    #print("a",np.hstack((q_kk,w_est[i,:])))
    Advitiy.set_est_State(np.hstack((q_kk,w_est[i,:])))
        
    #print(w_est[i,:]- w_true[i,:])
    #print(w_gyro[i,:]+b)
    #print(w_est[:,0]-w_gyro[:,0]-b[0])
    xmod[i]=np.linalg.norm(f) 

#plt.plot(l,p[:,0,0], 'r')
plt.xlabel('number of iterations')  
plt.ylabel('quaternion')  
#plt.plot(l,w_gyro[:,0], 'b')
#plt.plot(l,w_est[:,0], 'r')

plt.plot(l,q_true[:,1], 'b')
plt.plot(l,q[:,1], 'r')
