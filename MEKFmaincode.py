# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 10:14:38 2019
@author: Sanskriti
"""

from constants_1U import *
import testfunction
import scipy
from scipy.linalg import expm
import numpy as np
import qnv
import satellite
import sensor
from constants_1U import *

v_state0 = np.hstack((v_q0_BO,v_w0_BOB))
print("v_state0")
print(v_state0)                         
t0 = 0
#Make satellite object
#Advitiy = satellite.Satellite(v_state0,t0) 
#Advitiy.setMag_i(np.array([1,1,1]))
#Advitiy.setMag_b_m_c(np.array([1,1,1]))
#Advitiy.setPos(np.array([1,1,1]))
#Advitiy.setVel(np.array([1,2,1]))
#Advitiy.setMag_o(np.array([1,2,1]))
  #t0 from line 42 of main_code

#b=np.array([1,1,1]) #bias
b_e=np.array([0,0,0]) #estimated bias
#I3=np.matrix('1,0,0;0,1,0;0,0,1') #defining 3x3 identity matrix 
#I4=np.matrix('1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1') #defining 4x4 identity matrix
#I6=np.matrix('1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1')
#P_k=np.matrix('1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1') #initial state covariance matrix
q_kk=np.array([0,0,0,1])
I3=np.array([[1,0,0],[0,1,0],[0,0,1]]) #defining 3x3 identity matrix 
I4=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) #defining 4x4 identity matrix
I6=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
P_k=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]) #initial state covariance matrix
delta_q_kk=np.array([0.1,0.1,0.1,0.985])

#w_m_k=np.matrix('1,1,1').T

R=I3
#quat=qnv.quat2rotm(np.squeeze(np.asarray(q_kk)))

#q=np.asmatrix(quat)
#delta_b=b-b_e #delta b, diffence between gyro bias and estimated bias
#w_bib=w_m_k-b_e #estimated omega
#w_bob=w_bib-(q*w_oio)
#delta_theta=q_kk[0:3,0]/q_kk[3,0] #vector part of error quaternion normalised such that scalar part is equated to 1
#delta_x=np.array([delta_theta[0],delta_theta[1],delta_theta[2],delta_b[0],delta_b[1],delta_b[2]])  #error state vector 
#A=I3-qnv.skew(w_bob)

t=1


#R=I3
#v_mag_o=Advitiy.getMag_o()
#v_mag_o=sat.getMag_o()
#v_mag_b_m=Advitiy.getMag_b_m_c()
b=np.array([0.0001,0.0001,0.0001])
b_e=np.array([0,0,0])


def estimator(sat):
    sat.setGyroVarBias(np.array([0,0,0]))
    v_mag_o=sat.getMag_o()
    v_mag_b_m=sat.getMag_b_m_c()
    w_m_k=sensor.gyroscope(sat)
    w_oio=-v_w_IO_o
    delta_b=testfunction.delta_b_calc(b,b_e)
    
    w_bob=testfunction.w_bob_calc(w_m_k,q_kk,w_oio,b_e)
    #print(w_bob)
    phi=testfunction.phi_calc(w_bob)
    #print(phi)
    delta_x=testfunction.delta_x_calc(q_kk,delta_b)
    #print(delta_x)
    q_k1k=testfunction.propogate_quaternion(w_bob,q_kk)
    #print(q_k1k)
    x_k1k=testfunction.propogate_state_vector(phi,delta_x)
    #print(x_k1k)
    P_k1k=testfunction.propogate_covariance(phi,P_k)
    #print(P_k1k)
    v_mag_b_e=testfunction.calc_v_mag_b_e(v_mag_b_m,v_mag_o,q_k1k)
    #print(v_mag_b_e)
    y=testfunction.calc_y(v_mag_b_m,v_mag_b_e)
    #print(y)
    M_m=testfunction.calc_M_m(v_mag_b_e,q_k1k)
    #print(M_m)
    K=testfunction.calc_K(P_k1k,M_m,R)
    #print("K=")
    #print(K)   
    x_k1k1=testfunction.update_state_vector(K,y,x_k1k,M_m)
    #print("x=" )
    #print(x_k1k1)
    P_k1k1=testfunction.update_covariance(I6,K,M_m,P_k1k,R)
    #print("p=")
    #print(P_k1k1)
    q_k1k1=testfunction.update_quaternion(x_k1k,q_k1k)
    #print("qk1k1")
    #x_k1k1 = testfunction.reset_error_quaternion(x_k1k1)
    #print(q_k1k1)
    return q_k1k1, P_k1k1, x_k1k1
    
#print("result")