#from satellite import Satellite 
import math as mt
import numpy as np
#v_sun_o=sat.getSun_o()	#sun vector in orbit frame
#v_mag_o=sat.getMag_o()  #magnetic vector in orbit frame
#v_sun_b=sat.getSun_b_m()  #sun vector in body frame
#v_mag_b=sat.getMag_b_m()  #magnetic vector in body frame
def quest(v_sun_o,v_mag_o,v_sun_b_m,v_mag_b_m):

    r1=v_sun_o
    r2=v_mag_o
    b1=v_sun_b_m
    b2=v_mag_b_m
    
    a1=0.1  #weight of sun vector
    a2=0.9  #weight of magnetic vector
    r=np.cross(r1,r2)
    if r.all()==0:
        x=1;
    else:
        x=np.linalg.norm(r)
    r3=r/x              #define r3
    
    b=np.cross(b1,b2)
    if b.all()==0:
        y=1;
    else:
        y=np.linalg.norm(b)
    b3=b/y              #define b3
    
    p=np.cross(b3,r3)
    q=np.dot(b3,r3)
    r=a1*np.cross(b1,r1)
    s=a2*np.cross(b2,r2)
    rd=a1*np.dot(b1,r1)
    sd=a2*np.dot(b2,r2)
    
    alpha=(1+q)*(rd+sd)+np.dot(p,(r+s))
    beta=np.dot((b3+r3),(r+s))
    gamma=mt.sqrt((alpha**2)+(beta**2))
    q_bob = np.zeros(4)
    if alpha>0:      #for alpha greater than 0
        coef1=0.5/mt.sqrt(gamma*(gamma+alpha)*(1+q)) #normalising factor
        q_v=((gamma+alpha)*p)+(beta*(b3+r3))
        q_1=coef1*q_v[0]     #vector
        q_2=coef1*q_v[1]
        q_3=coef1*q_v[2]
        q_4=coef1*(gamma+alpha)*(1+q)      #scalar
        q_bob[0] = q_1
        q_bob[1] = q_2
        q_bob[2] = q_3
        q_bob[3] = q_4
        
        
       
    if alpha<=0:    #for alpha less than or equals 0
        coef2=0.5/mt.sqrt(gamma*(gamma-alpha)*(1+q))  #normalising factor
        q_v=((beta)*p)+((gamma-alpha)*(b3+r3))
        q_1=coef2*q_v[0]    #vector
        q_2=coef2*q_v[1]
        q_3=coef2*q_v[2]
        q_4=coef2*beta*(1+q)     #scalar
        q_bob[0] = q_1        
        q_bob[1] = q_2
        q_bob[2] = q_3
        q_bob[3] = q_4

    return q_bob
print(quest([1,0,0],[0,1,0],[0,1.09,0],[0,0,1]))    #to check test case
    
    
