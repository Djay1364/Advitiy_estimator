#v_sun_o=sat.getSun_o()	#sun vector in orbit frame
#v_mag_o=sat.getMag_o()  #magnetic vector in orbit frame
#v_sun_b=sat.getSun_b_m()  #sun vector in body frame
#v_mag_b=sat.getMag_b_m()  #magnetic vector in body frame
def estimator_default(v_sun_o,v_mag_o,v_sun_b_m,v_mag_b_m):

    r1=v_sun_o
    r2=v_mag_o
    b1=v_sun_b_m
    b2=v_mag_b_m
    q=[0,0,0,1]
    return q  
#print(estimator_default(0,0,0,0))      