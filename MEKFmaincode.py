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
q_kk=np.array([1/2**0.5,0,0,1/2**0.5])#it is supposed to come from QUEST
I3=np.array([[100,0,0],[0,100,0],[0,0,100]]) #defining 3x3 identity matrix 
I4=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) #defining 4x4 identity matrix
I6=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
P_k=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]) #initial state covariance matrix
delta_q_kk=np.array([0,0,0,1])

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
#b=np.array([0.0001,0.0001,0.0001])
#b_e=np.array([0,0,0])
 
def estimator(sat):

    v_glob_est_state = sat.get_glob_est_State() #### comes from the previous estimate
    q_kk = v_glob_est_state[0:4] 
    P_kk = sat.getErrCovariance() #### comes from the previous estimate
    b=sat.getGyroVarBias()#### gyro rate bias
    b_e=sat.getGyroEstBias()#### estimated gyro rate bias
    delta_b=testfunction.delta_b_calc(b,b_e)
    
    v_mag_o=sat.getMag_o() #### magnetic field in orbit frame
    v_mag_b_m, v_n_m=sensor.magnetometer(sat) #sat.getMag_b_m_c() 
    v_mag_b_m=v_mag_b_m/(np.linalg.norm(v_mag_b_m))
    w_m_k=sensor.gyroscope(sat) #### angular velocity of body w.r.t. inertial expressed in Body
    w_oio=-v_w_IO_o #### angular velocity of orbit w.r.t. inertial expressed in body
    
    x_kk=testfunction.delta_x_calc(delta_q_kk,delta_b)#### local state
    #q_true_bo = sat.getQ_BO() #### true quaternion
    w_bob = testfunction.w_bob_calc(w_m_k,q_kk,w_oio,b_e)
    
    phi=testfunction.phi_calc(w_bob) #### state transition matrix of local estimates
   
    q_k1k=testfunction.propogate_quaternion(w_bob,q_kk) #### quaternion propagation using quaternion kinematics
    
    x_k1k=np.dot(phi,x_kk)#testfunction.propogate_state_vector(phi,x_kk) #### local state propagation 
    
    P_k1k=testfunction.propogate_covariance(phi,P_kk) #### error covariance propagation

    v_mag_b_e=testfunction.calc_v_mag_b_e(v_mag_b_m,v_mag_o,q_k1k) ####
    
    y=testfunction.calc_y(v_mag_b_m,v_mag_b_e)
    
    M_m=testfunction.calc_M_m(v_mag_b_e,q_k1k) #### No need to pass q_k1k
    
    K=testfunction.calc_K(P_k1k,M_m,R)#np.matrix([[1,1,1],[1,1,1],[1,1,1],[1,1,1],[1,1,1],[1,1,1]]) 
    
    x_k1k1=testfunction.update_state_vector(K,y,x_k1k,M_m)
    
    P_k1k1=testfunction.update_covariance(I6,K,M_m,P_k1k,R)
    
    q_k1k1=testfunction.update_quaternion(x_k1k,q_k1k)
    
    w_est=w_m_k+x_k1k1[3:6]-b
    
    sat.set_loc_est_State(x_k1k1)
    sat.setErrCovariance(P_k1k1)
    sat.set_glob_est_State(np.hstack((q_k1k1,w_est)))
    sat.setGyroEstBias(b-x_k1k1[3:6])
    print(x_k1k1[3:6])
    
    return q_k1k1, P_k1k1, x_k1k1
    