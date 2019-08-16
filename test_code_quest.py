import ddt
import singlepoint_estimator
import unittest	#testing library
import numpy as np
from ddt import ddt,file_data,unpack,data
@ddt
class Testquest(unittest.TestCase):
        @file_data("test_quest.json")
        @unpack
        def test_quest(self,value):
            given1,given2,given3,given4,expected = np.asarray(value[0]),np.asarray(value[1]),np.asarray(value[2]),np.asarray(value[3]),np.asarray(value[4])
            result = singlepoint_estimator.quest(given1,given2,given3,given4)
            self.assertTrue(np.allclose([result[0],result[1],result[2],result[3]],expected))     #asarray to convert list to array

if __name__=='__main__':
	unittest.main(verbosity=2)	    	    

            
	
	