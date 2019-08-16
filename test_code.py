import ddt
import sun_sensor_modelling
import unittest	#testing library
import numpy as np
from ddt import ddt,file_data,unpack,data
@ddt
class TestLatLon(unittest.TestCase):
	@file_data("test_ad_convertor.json")
	@unpack
	def test_ad_convertor(self,value):
	    given1,given2,given3,given4,given5,given6,given7,given8,given9,expected = np.asarray(value[0]),np.asarray(value[1]),np.asarray(value[2]),np.asarray(value[3]),np.asarray(value[4]),np.asarray(value[5]),np.asarray(value[6]),value[7],value[8],value[9]	#asarray to convert list to array
	    result = sun_sensor_modelling.ad_Convertor(given1,given2,given3,given4,given5,given6,given7,given8,given9)
	    self.assertTrue(np.allclose([result[0],result[1],result[2],result[3],result[4],result[5]],expected))
	
	@file_data("test_light.json")
	@unpack
	def test_light(self,value):
	    given1,given2,expected=np.asarray(value[0]),value[1],np.asarray(value[2])
	    result=sun_sensor_modelling.light(given1,given2)
	    self.assertTrue(np.allclose([result[0],result[1],result[2],result[3],result[4],result[5]],expected))
	    
	       	  	   	       	  	   	    
	@file_data("test_sv_calc.json")
	@unpack
	def test_sv_calc(self,value):
	    given1,given2,given3,expected = np.asarray(value[0]),value[1],value[2],np.asarray(value[3])	#asarray to convert list to array
	    result = sun_sensor_modelling.calc_SV(given1,given2,given3)
	    self.assertTrue(np.allclose([result[0],result[1],result[2]],expected))
	   	  	   	       	  	   	    

if __name__=='__main__':
	unittest.main(verbosity=2)	   	   	   	   
