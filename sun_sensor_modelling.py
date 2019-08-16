import numpy as np
def ad_Convertor(sun_vector,v_S1,v_S2,v_S3,v_S4,v_S5,v_S6,SS_GAIN,SS_QUANTIZER):    #works as adc convertor, quantizes the readings of sensor
        u=1/(SS_QUANTIZER-1)                                                          #for more you can refer to documentation of this code.
                                                                       #takes input of sun vector, gain, quantizer, sensor vectors to give modelled quantized readings of sun sensor
        ss1=np.dot(sun_vector,v_S1)
        if ss1<0:
            ss1=0    
        ss1=(u)*(round(ss1/u))*SS_GAIN
        ss2=(np.dot(sun_vector,v_S2))
        if ss2<0:
            ss2=0    
        ss2=(u)*(round(ss2/u))*SS_GAIN
        ss3=(np.dot(sun_vector,v_S3))
        if ss3<0:
            ss3=0    
        ss3=(u)*(round(ss3/u))*SS_GAIN
        ss4=(np.dot(sun_vector,v_S4))
        if ss4<0:
            ss4=0    
        ss4=(u)*(round(ss4/u))*SS_GAIN
        ss5=(np.dot(sun_vector,v_S5))
        if ss5<0:
            ss5=0    
        ss5=(u)*(round(ss5/u))*SS_GAIN
        ss6=(np.dot(sun_vector,v_S6))
        if ss6<0:
            ss6=0    
        ss6=(u)*(round(ss6/u))*SS_GAIN
        return np.array([ss1,ss2,ss3,ss4,ss5,ss6])                            #stores quantized values of all readings.
f=np.array([0,0,0,0,0,0])
#a=np.cos(MAX_ANGLE*np.pi/180)*SS_GAIN
def light(ss,a):                #calucates the flag value-light, by giving boolean for each sun sensor according to thresold given as input.
    for i in range(6):
        if ss[i]>a:
            f[i]=1
        else:
            f[i]=0
    return f    
def calc_SV(ss,dark,light):             #calculates back the sun vector using the sensor readings.
      v_sun_m=np.array([0,0,0])            
      for i in range(6):
           if ss[i]==0:
               dark=dark+1
      if dark==6:
          light=0
      if light==1:     
         for n in range(3):
                  m=n*2;
                  if ss[m]>=ss[m+1]:
                        v_sun_m[n]=1.0*ss[m]
                  else:
                     v_sun_m[n]=-1.0*ss[m+1]
      mode=pow((pow(v_sun_m[0],2)+pow(v_sun_m[1],2)+pow(v_sun_m[2],2)),0.5)                       
      return v_sun_m/mode        #gives the unit sun vector to be used in quest.    
