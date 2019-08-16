import numpy as np
v_x_mag=input()  #taking input from igrf x-component
v_y_mag=input()                          #y-component
v_z_mag=input()                          #z-component
mag_bias=0.0003   #bias taken from data sheet of magnetometer used in pratham
def v_mag(v_x_mag,v_y_mag,v_z_mag,mag_bias):    
          x=v_x_mag-np.random.normal()-mag_bias  #adding random error and bias in each component.
          y=v_y_mag-np.random.normal()-mag_bias
          z=v_z_mag-np.random.normal()-mag_bias
          v_mag_m=np.array([x,y,z])
          return v_mag_m
