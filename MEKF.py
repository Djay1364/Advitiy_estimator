from scipy.linalg import expm
import numpy as np
import qnv1
I3=np.matrix('1,0,0;0,1,0;0,0,1') #defining 3x3 identity matrix 
I4=np.matrix('1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1') #defining 4x4 identity matrix
delt=0.1 #delta t, small time over which we integrate. 
w_m_k=np.array([1,2,3]) #omega measured by gyroscope for kth time
w_m_k1=np.array([1,2,3]) #omega measured by gyroscope for k+1th time
b=np.array([1,1,1])     #gyro bias
b_e=np.array([.9,.9,.9]) #estimated bias for first iteration
sigma_r_sq=1   #noise covariance for discrete system
sigma_w_sq=1   #noise covariance for discrete system
P_k=np.matrix('[1,0,0],[0,1,0],[0,0,1]') #initial state covariance matrix
q_k=np.array([0,0,0,1])  #initial quaternion taken (we will actually get this from quest)
def mekf_propogate(x,w_m_k,w_m_k1,q_k,b_e,P_k): #a function that takes in state vector,gyro readings,quaternion from previous iteration( from quest in case of first iteration), estimated bias and state covariance matrix from previous iteration.
                   delta_b=b-b_e #delta b, diffence between gyro bias and estimated bias
                   w_e=w_m_k-b_e #estimated omega
                   #delta_q=np.array([0,0,0,1]) #error quaternion
                   delta_theta=delta_q[0:3]/delta_q[3] #vector part of error quaternion normalised such that scalar part is equated to 1
                   delta_x=np.array([delta_theta[0],delta_theta[1],delta_theta[2],delta_b[0],delta_b[1],delta_b[2]])  #error state vector 
                   A=I3-qnv1.skew(w_e)*delt+0.5*pow(qnv1.skew(w_e),2)*pow((delt),2)    
                   B=-I3*delt+qnv1.skew(w_e)*pow(delt,2)*0.5-pow(qnv1.skew(w_e),2)*pow(delt,3)/6
                   phi=np.matrix([[A[0,0],A[0,1],A[0,2],B[0,0],B[0,1],B[0,2]],[A[1,0],A[1,1],A[1,2],B[1,0],B[1,1],B[1,2]],[A[2,0],A[2,1],A[2,2],B[2,0],B[2,1],B[2,2]],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
                   #computation of state transition matrix (phi)
                   Q11=sigma_r_sq*delt*I3+sigma_w_sq*(I3*pow(delt,3)/3+2*pow(delt,5)*pow(qnv1.skew(w_e),2)/120)
                   Q12=-sigma_w_sq*(I3*pow(delt,2)/2-pow(delt,3)*qnv1.skew(w_e)/6+pow(delt,4)*pow(qnv1.skew(w_e),2)/24)
                   Q22=sigma_w_sq*delt*I3
                   Q21=Q12.T
                   Q=np.matrix([[Q11[0,0],Q11[0,1],Q11[0,2],Q12[0,0],Q12[0,1],Q12[0,2]],[Q11[1,0],Q11[1,1],Q11[1,2],Q12[1,0],Q12[1,1],Q12[1,2]],[Q11[2,0],Q11[2,1],Q11[2,2],Q12[2,0],Q12[2,1],Q12[2,2]],[Q21[0,0],Q21[0,1],Q21[0,2],Q22[0,0],Q22[0,1],Q22[0,2]],[Q21[1,0],Q21[1,1],Q21[1,2],Q22[1,0],Q22[1,1],Q22[1,2]],[Q21[2,0],Q21[2,1],Q21[2,2],Q22[2,0],Q22[2,1],Q22[2,2]]])
                   #computation of discrete noise variance matrix
                   #following are the propogation equations
                   b_e_k1=b_e #propogation of bias estimate
                   w_e_k1=w_m_k1-b_e_k1 #estimated angular rate propogation
                   w_bar=(w_e+w_e_k1)/2 
                   q_k1=(I4+qnv1.omega(w_bar)*delt+0.5*qnv1.omega(w_bar)*qnv1.omega(w_bar)*pow(delt,2)+(qnv1.omega(w_e_k1)*qnv1.omega(w_e)-qnv1.omega(w_e)*qnv1.omega(w_e_k1))*pow(delt,2)/48)*q_k #propogated quaternion
                   P_k1=phi*P_k*(phi.T)+Q #propogated state covariance matrix
                   
                   return delta_x,b_e_k1,w_e_k1,q_k1,P_k1,delta_theta

                   

#measurement update
v_sun_b=np.array([1,1,1]) 
def measurement_update(delta_x,b_e_k1,w_e_k1,q_k1,P_k1,delta_theta):
    
    
   