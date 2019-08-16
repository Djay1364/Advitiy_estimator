from scipy.linalg import expm
import numpy as np
import qnv1
I3=np.matrix('1,0,0;0,1,0;0,0,1') #defining 3x3 identity matrix 
I4=np.matrix('1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1') #defining 4x4 identity matrix
P_k=np.matrix('[1,0,0],[0,1,0],[0,0,1]') #initial state covariance matrix
q_(k|k)=np.matrix('0,0,0,1')
delta_b=b-b_e #delta b, diffence between gyro bias and estimated bias
w_bib=w_m_k-b_e #estimated omega
w_bob=w_bib-w_oio*(qnv1.quat2rotm(q_(k|k)))
delta_theta=q_(k|k)[0:3]/q_(k|k)[3] #vector part of error quaternion normalised such that scalar part is equated to 1
delta_x=np.array([delta_theta[0],delta_theta[1],delta_theta[2],delta_b[0],delta_b[1],delta_b[2]])  #error state vector 
A=I3-qnv1.skew(w_bob)
B=-I3
T=1
C=np.matrix([[A[0,0],A[0,1],A[0,2],B[0,0],B[0,1],B[0,2]],[A[1,0],A[1,1],A[1,2],B[1,0],B[1,1],B[1,2]],[A[2,0],A[2,1],A[2,2],B[2,0],B[2,1],B[2,2]],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
Q_k=np.matrix([sigma_r_sq,0,0,0,0,0],[0,sigma_r_sq,0,0,0,0],[0,0,sigma_r_sq,0,0,0],[0,0,0,sigma_w_sq,0,0],[0,0,0,0,sigma_w_sq,0],[0,0,0,0,0,sigma_w_sq])
phi=np.linal.matrix_power(scipy.linalg.expm(A),T)                                      
def propogate_quaternion(w_bob,q_(k|k),T):
    
    
    return 
def propogate_state_vector(phi,delta_x):
    return phi*delta_x
def propogate_covariance(phi,Q_k,P_k):
    return phi*P_k*(phi.T)+Q_k
    
def update_quaternion(phi,delta_x,q_(k|k)):
    delta_q=(delta_x[0:3],0)
    return qnv1.quatMultiplyNorm(delta_q,q_(k|k)) 
    
    
    
                       