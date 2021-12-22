import unittest, os, sys, inspect
import numpy as np

curr_dir = os.path.dirname( os.path.abspath( inspect.getfile( inspect.currentframe()) ) );
parent_dir = os.path.dirname(curr_dir); 
sys.path.insert(0,parent_dir); 
import waves

class TestWaves(unittest.TestCase):
    
    def test_set_wvnum(self):
        mywave = waves.internal(); 
        mywave.medium = {'N2':1,'lat':0}; # internal wave with no rotation
        mywave.wave_vector = {'omega':np.nan,'kx':0,'ky':1e-4,'kz':1e-4}
        print(mywave.wave_vector)
        self.assertEqual(0,0)

if __name__ == '__main__':
    unittest.main()

