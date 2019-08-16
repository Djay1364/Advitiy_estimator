# -*- coding: utf-8 -*-
"""
Created on Sat May 11 11:10:40 2019

@author: Sanskriti
"""
import numpy as np
import matplotlib.pyplot as plt

def predict(x, P, F, Q):
    xpred = np.dot(F,x)
    Ppred = np.dot(np.dot(F,P),(np.transpose(F))) + Q
    return xpred, Ppred

def innovation(xpred, Ppred, z, H, R):
    s = np.dot(H,xpred)
    nu = z - s
    S = R + np.dot(np.dot(H,Ppred),(np.transpose(H)))
    return nu, S

def innovation_update(xpred, Ppred, nu, S, H):
    K = np.dot(np.dot(Ppred,(np.transpose(H))),np.linalg.inv(S))
    xnew = xpred + np.dot(K,nu)
    Pnew = Ppred - np.dot(np.dot(K,S),(np.transpose(K)))
    return K, xnew, Pnew
    
delT = 1
t=5
l=np.linspace(0,t,t)
F = np.array([[1, delT],[0,1]])
H = np.array([1, 0])
x = np.array([0,10])
P = np.array([[10, 0],[0 ,10]])
Q = np.array([[1, 1],[1, 1]])
r = np.array([1])
z =np.array([2.5,1,4,2.5,5.5])
xpred=np.zeros((t,2))
xnew=np.zeros((t,2))
Pnew=np.zeros((t,2,2))
xmod=np.zeros((t))
mur=np.zeros((t))
xtmod=np.zeros((t))
m=np.zeros((t,2))
kalman=np.zeros((t,2))
xt=np.zeros((t,2))
for i in range(t):
    print(x)
    X, R = predict(x, P, F, Q)
    nu, S = innovation(X, R, z[i], H, R)
    K, x, P = innovation_update(X, R, nu, S, H)
    mur[i] = nu
    kalman[i,:]=K
    m[i,:]=np.dot(K,z[i])
    xpred[i,:]=X
    xnew[i,:]=x
    Pnew[i,:,:]=R
    xmod[i]=np.linalg.norm(x)
    xtmod[i]=np.linalg.norm(xpred)
print(i)
#print(xnew)
#print(Pnew)
#plt.plot(l,Pnew[:,0,0])
#plt.plot(l,xnew[:,0])
#plt.plot(l,z)
#plt.plot(l,xt[:,0])
#plt.plot(l,kalman)
#plt.show()    
    
    

    