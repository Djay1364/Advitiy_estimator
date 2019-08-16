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
control_step=1
model_step=1
t0=0
v_state0 = np.hstack((v_q0_BO,v_w0_BOB))
Advitiy = satellite.Satellite(v_state0,t0) 
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
end=5000
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
q_true=np.array([0,0,0,1])
q=np.zeros((end,4))
RMSE=np.zeros((end,3))
p=np.zeros((end,6,6))
x=np.zeros((end,6))
xmod=np.zeros((end))
for i in range(end):
####w_true_bob = np.array([0.1,0.1,0.1])
####q_true_bo = yet to be defined 
####Advitiy.setW_BO_b(w_true_bob)
####Advitiy.setQ_BO(q_true_bo)
    Advitiy.setGyroVarBias(np.array([0.001,0.001,0.001]))
    w_m_k=sensor.gyroscope(Advitiy)
    #w_t_k=satellite.getW_BO_b(Advitiy)
    #w_t_bib=satellite.getW_BI_b(Advitiy)
    Advitiy.setMag_i(m_magnetic_field_i[i,1:4]*1e-9)
    Advitiy.setMag_b_m_c(sensor.magnetometer(Advitiy))
    Advitiy.setPos(m_sgp_output[i,1:4]*1e-6)
    Advitiy.setVel(m_sgp_output[i,4:7]*1e-6)
    position[i]=Advitiy.getPos()
    velocity[i]=Advitiy.getVel()
    #print(i)
    b_e=b-x[i-1,3:6]
    #print(b_e)
    q[i,:], p[i :], f=MEKFmaincode.estimator(Advitiy)
    x[i, :]=f
    print('a',x[i,2])
    RMSE[i,:]=np.sqrt(((b-f[3:6])**2)/end)
    q_kk=q[i,:]
    delta_x=[0,0,0,f[3:6]]
    #print(i)
    #print(f[3:6])
    P_k=p[i,:]
    #w_true[i,:]=w_t_bib
    w_gyro[i,:]=w_m_k
    w_est[i,:]=w_m_k+x[i,3:6]-b
    #print(w_gyro[i,:]+b)
    #print(w_est[:,0]-w_gyro[:,0]-b[0])
    xmod[i]=np.linalg.norm(f) 

plt.plot(l,velocity)

#plt.plot(l,p[:,0,0], 'r')
#plt.xlabel('number of iterations')  
#plt.ylabel('z component of angular rates')  
#plt.plot(l,w_est[:,2], 'b')
#plt.plot(l,w_gyro[:,2], 'r')
    
#plt.plot(l,q[:,0],'b')
#plt.plot(l,q[:,1],'r')
#plt.plot(l,q[:,2],'g')
#plt.plot(l,q[:,3],'y')
#plt.plot(l,x[:,3])
#plt.plot(l,x[:,4])
#plt.plot(l,x[:,5])