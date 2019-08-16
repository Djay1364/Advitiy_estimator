#constants used in sensor modelling
import numpy as np
#sunsensor
SS_GAIN=3.3 #gain factor as the readings from sunsensor are in voltage.
SS_QUANTIZER=pow(2,10)-1  #because we have 10 bit adc
MAX_ANGLE=55            #taken from pratham's simulink model,can be adjusted, it is the 1/2 of threshold value of angle which will be taken in FOV
v_S1=np.array([1,0,0])   #vectors perpendicular to sensors.
v_S2=np.array([-1,0,0])
v_S3=np.array([0,1,0])
v_S4=np.array([0,-1,0])
v_S5=np.array([0,0,1])
v_S6=np.array([0,0,-1])
#maximum voltage reading observed= 3.3
#minimum voltage reading observed= 0
#minimum voltage at which sensor is taken to be in light regiom is 1.894

#magnetometer
mag_bias=0.0003
