# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 11:25:01 2019
@author: Sanskriti
"""

import scipy
from scipy.linalg import expm
import numpy as np
import qnv
import satellite as sat
from constants_1U import *
I3=np.array([[1,0,0],[0,1,0],[0,0,1]]) #defining 3x3 identity matrix 
I4=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) #defining 4x4 identity matrix
I6=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
P_k=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]) #initial state covariance matrix
q_kk=np.array([0,0,0,1])
b=np.array([0,0,0])
b_e=np.array([0,0,0])
#w_m_k=np.matrix('1,1,1').T
w_oio=np.array([1,1,1])
R=10*I3
#quat=qnv.quat2rotm(np.squeeze(np.asarray(q_kk)))

#q=np.asmatrix(quat)
#delta_b=b-b_e #delta b, diffence between gyro bias and estimated bias
#w_bib=w_m_k-b_e #estimated omega
#w_bob=w_bib-(q*w_oio)
#delta_theta=q_kk[0:3,0]/q_kk[3,0] #vector part of error quaternion normalised such that scalar part is equated to 1
#delta_x=np.array([delta_theta[0],delta_theta[1],delta_theta[2],delta_b[0],delta_b[1],delta_b[2]])  #error state vector 
#A=I3-qnv.skew(w_bob)
B=-I3
t=0.1
h=0.01
sigma_r_sq=1e-9
sigma_w_sq=1e-9
#C=np.matrix([[A[0,0],A[0,1],A[0,2],B[0,0],B[0,1],B[0,2]],[A[1,0],A[1,1],A[1,2],B[1,0],B[1,1],B[1,2]],[A[2,0],A[2,1],A[2,2],B[2,0],B[2,1],B[2,2]],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
Q_k=np.array([[sigma_r_sq,0,0,0,0,0],[0,sigma_r_sq,0,0,0,0],[0,0,sigma_r_sq,0,0,0],[0,0,0,sigma_w_sq,0,0],[0,0,0,0,sigma_w_sq,0],[0,0,0,0,0,sigma_w_sq]])
#phi=scipy.linalg.expm(A*t)  
def delta_b_calc(b,b_e):
    delta_b=b-b_e
    return delta_b

def w_bob_calc(w_m_k,q_true_bo,w_oio,b_e):
    R=qnv.quat2rotm(q_true_bo)
    w_bib=w_m_k-b_e #estimated omega
    #w_bib=sat.getW_BI_b()
    #v_bias_var = sat.getGyroVarBias()
    w_bob=w_bib-np.dot(R,w_oio)
    return w_bob

def phi_calc(w_bob):
    A=-qnv.skew(w_bob)
    C=np.matrix([[A[0,0],A[0,1],A[0,2],-1,0,0],[A[1,0],A[1,1],A[1,2],0,-1,0],[A[2,0],A[2,1],A[2,2],0,0,-1],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
    phi=scipy.linalg.expm(C*h)  
    return phi

def delta_x_calc(q_kk,delta_b):
    delta_theta=2*(q_kk[0:3]/q_kk[3]) #vector part of error quaternion normalised such that scalar part is equated to 1
    delta_x=np.array([delta_theta[0],delta_theta[1],delta_theta[2],delta_b[0],delta_b[1],delta_b[2]])  #error state vector 
    return delta_x  
                              
def propogate_quaternion(w_bob,q_kk):
    #state at t = t0	
	#rk-4 routine (updating satellite class state with obtained state at every step of rk4 routine)
	#first step of rk4 routine
    dq = h*qnv.quatDerBO(q_kk,w_bob)
    q_k1k = q_kk + dq
    return q_k1k

def propogate_state_vector(phi,x_kk): 
    
    R=np.dot(phi,x_kk)+h*np.dot()
    return R

def propogate_covariance(phi,P_k,w_est):
    w = np.linalg.norm(w_est)
    w_hat = qnv.skew(w_est)
    Q11 = (sigma_u**2)*h*I3 + (sigma_v**2)*(I3*(h**3/3) + ((w*h)**3/3 + 2*np.sin(w*h) - 2*w*h)*(w_hat**2)/w**5)
    Q12 = -(sigma_v**2)*(I3*(h**2/3) - (w*h-np.sin(w*h))*w_hat/w**3 + ((w*h)**2/2 + np.cos(w*h) -1)*(w_hat**2)/w**4)
    Q22 = (sigma_v**2)*h*I3
    
    Q_k = np.vstack((np.hstack((Q11,Q12)),np.hstack((Q12.T,Q22))))
    
    return np.dot(phi,np.dot(P_k,(phi.T)))+Q_k

#v_mag_b_m=sat.getMag_b_m_c()  
#v_mag_b_m=np.matrix('1,1,1').T
#v_mag_o=sat.getMag_o()
#v_mag_o=np.matrix('1,2,1').T
#q_k1k = np.matrix('1,0,0,0').T


def calc_v_mag_b_e(v_mag_b_m,v_mag_o,q_k1k):
    v_mag_o = v_mag_o/(np.linalg.norm(v_mag_o))
    quatk1km=qnv.quat2rotm(q_k1k)
    v_mag_b_e=np.dot(quatk1km,v_mag_o)
    v_mag_b_e=v_mag_b_e/(np.linalg.norm(v_mag_b_e))
    return v_mag_b_e

def calc_y(v_mag_b_m,v_mag_b_e):
    #print(v_mag_b_e)# v_mag_b_m)
    #print("v_mag_b_m",v_mag_b_m)
    #print("v_mag_b_e",v_mag_b_e)
    y=v_mag_b_m-v_mag_b_e
    return y

def calc_M_m(v_mag_b_e,q_k1k):
    D=qnv.skew(v_mag_b_e) 
    M_m=np.array([[D[0,0],D[0,1],D[0,2],0,0,0],[D[1,0],D[1,1],D[1,2],0,0,0],[D[2,0],D[2,1],D[2,2],0,0,0]]) #use np.hstack or np.vstack for elegance
    return M_m

def calc_K(P_k1k,M_m,R):
    K=np.dot(P_k1k,np.dot(M_m.T,np.linalg.inv(np.dot(M_m,np.dot(P_k1k,M_m.T))+ R))) #check
    return K
    
def update_quaternion(delta_x_k1k,q_kk):
    delta_q=np.hstack((delta_x_k1k[0:3]/2,np.array([1])))
    #print("delta_x_kk=")
    #print(delta_x_k1k)
    #print("q_kk=")
    #print(q_kk)
    return qnv.quatMultiplyNorm(delta_q,q_kk)

def update_state_vector(K,y,x_k1k,M_m):
    return np.dot(K,y-np.dot(M_m,x_k1k))

def update_covariance(I6,K,M_m,P_k1k,R):
    return np.dot(np.dot((I6-np.dot(K,M_m)),P_k1k),(I6-np.dot(K,M_m)).T)+np.dot(np.dot(K,R),K.T)

def reset_error_quaternion(x_k1k1):
    delta_th = np.array([0,0,0])
    x_reset = np.hstack((delta_th,x_k1k1[3:6]))
    return x_reset
    