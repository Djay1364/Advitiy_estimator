import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math 
import numpy as np

val=np.linspace(-3,3,30)
guess=[-2,1]
str=np.linspace(-1,1,30)
a=[0]
b=[0]
for i in range(0,9):
    a.append(a[i]+(1/9)*np.cos(10*np.pi/180))
    b.append(b[i]-(1/9)*np.sin(10*np.pi/180))
print(a[9],b[9])
plt.plot(a,b,linewidth=2.5,color='black')
for j in str:
    x_val=[]
    y_val=[]
    for i in val:
        def equations(p):
            X, Y = p
            K=.015
            return(K/2*(np.log((pow((X-a[0]),2)+pow(Y-b[0],2))*(pow((X-a[1]),2)+pow(Y-b[1],2))*(pow((X-a[2]),2)+pow(Y-b[2],2))*(pow((X-a[3]),2)+pow(Y-b[3],2))*(pow((X-a[4]),2)+pow(Y-b[4],2))*(pow((X-a[5]),2)+pow(Y-b[5],2))*(pow((X-a[6]),2)+pow(Y-b[6],2))*(pow((X-a[7]),2)+pow(Y-b[7],2))*(pow((X-a[8]),2)+pow(Y-b[8],2))*(pow((X-a[9]),2)+pow(Y-[9],2))))+Y+j,X+i)
    
        X,Y= fsolve(equations,guess)
        x_val.append(X)
        y_val.append(Y)
        guess[0]=X
        guess[1]=Y
        
    plt.plot (x_val, y_val)
    plt.show()
    plt.savefig('strmlines.png')