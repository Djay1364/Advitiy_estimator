import scipy
from scipy.linalg import expm
import numpy as np
import qnv
import satellite as sat
I3=np.matrix('1,0,0;0,1,0;0,0,1') #defining 3x3 identity matrix 
I4=np.matrix('1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1') #defining 4x4 identity matrix
I6=np.matrix('1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1')
P_k=np.matrix('1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1') #initial state covariance matrix
q_kk=np.matrix('0,0,0,1').T
b=np.matrix('1,1,1').T
b_e=np.matrix('0,0,0').T
w_m_k=np.matrix('1,1,1').T
w_oio=np.matrix('1,1,1').T
quat=qnv.quat2rotm(np.squeeze(np.asarray(q_kk)))

q=np.asmatrix(quat)
delta_b=b-b_e #delta b, diffence between gyro bias and estimated bias
w_bib=w_m_k-b_e #estimated omega
w_bob=w_bib-(q*w_oio)
delta_theta=q_kk[0:3,0]/q_kk[3,0] #vector part of error quaternion normalised such that scalar part is equated to 1
delta_x=np.array([delta_theta[0],delta_theta[1],delta_theta[2],delta_b[0],delta_b[1],delta_b[2]])  #error state vector 
A=I3-qnv.skew(w_bob)
B=-I3
t=1
sigma_r_sq=1
sigma_w_sq=1
C=np.matrix([[A[0,0],A[0,1],A[0,2],B[0,0],B[0,1],B[0,2]],[A[1,0],A[1,1],A[1,2],B[1,0],B[1,1],B[1,2]],[A[2,0],A[2,1],A[2,2],B[2,0],B[2,1],B[2,2]],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
Q_k=np.matrix([[sigma_r_sq,0,0,0,0,0],[0,sigma_r_sq,0,0,0,0],[0,0,sigma_r_sq,0,0,0],[0,0,0,sigma_w_sq,0,0],[0,0,0,0,sigma_w_sq,0],[0,0,0,0,0,sigma_w_sq]])
phi=scipy.linalg.expm(A*t)  
R=I3                                
def propogate_quaternion(w_bob,q_kk):
    #state at t = t0	
	
	#rk-4 routine (updating satellite class state with obtained state at every step of rk4 routine)
	#first step of rk4 routine
    dq = t*qnv.quatDerBO(q_kk,w_bob)
	
    q_k1k = q_kk + dq
    
    return q_k1k
def propogate_state_vector(phi,delta_x):
    return phi*delta_x
def propogate_covariance(phi,P_k):
    return phi*P_k*(phi.T)+Q_k

#v_mag_b_m=sat.getMag_b_m_c()  
v_mag_b_m=np.matrix('1,1,1').T
#v_mag_o=sat.getMag_o()
v_mag_o=np.matrix('1,2,1').T
q_k1k = np.matrix('1,0,0,0').T
quatk1k=qnv.quat2rotm(np.squeeze(np.asarray(q_kk)))
quatk1km=np.asmatrix(quatk1k)
v_mag_b_e=quatk1km*v_mag_o
y=v_mag_b_m-v_mag_b_e  #check again
D=qnv.skew(v_mag_b_e)
M_m=np.matrix([[D[0,0],D[0,1],D[0,2],0,0,0],[D[1,0],D[1,1],D[1,2],0,0,0],[D[2,0],D[2,1],D[2,2],0,0,0]])
K=P_k1k*M_m.T*(M_m*P_k1k*M_m.T+R).I 
   
def update_quaternion(phi,delta_x_k1k,q_kk):
    delta_q=(delta_x_k1k[0:3],0)
    return qnv.quatMultiplyNorm(delta_q,q_kk)

def update_state_vector(K,y,x_k1k,M_m):
    return x_k1k+K*(y-M_m*x_k1k)

def update_covariance(I6,K,M_m,P_k1k):
    return (I6-K*M_m)*P_k1k
 
    
