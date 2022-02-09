# TESTS FOR Coordinates Module


from pyacs.lib import coordinates as COO
from pyacs.lib import units as UNITS
import numpy as np

def test_geo2XYZ():
    
    X,Y,Z=COO.geo2xyz( -112.567,43.36,27.2 ,unit='dec_deg')
    # gamit convertc results  -1782430.03490 -4288974.18902  4356684.78382
    assert np.sqrt((X - -1782430.03490)**2) < 0.00002
    assert np.sqrt((Y - -4288974.18902)**2) < 0.00002
    assert np.sqrt((Z -  4356684.78382)**2) < 0.00002

def test_XYZ2geo():
    long,lat,he=COO.xyz2geo(918129.451, -4346071.255 , 4561977.839 )
    # gamit tform results 45.95580023  281.92863291      200.9022
    assert np.sqrt((np.remainder(np.degrees(long),360.) - 281.92863291)**2) < 0.000001
    assert np.sqrt((np.degrees(lat) - 45.95580023)**2) < 0.000001
    assert np.sqrt((he -  200.9022)**2) < 0.0001
    
def test_XYZ2geosphere_and_degminsec():
 
    # tform results for XYZ 918129.451, -4346071.255 , 4561977.839
    # N45 45 48.48702 W 78  4 16.92154 6367333.7313
  
    rlong,rlat,h=COO.xyz2geospheric(918129.451, -4346071.255 , 4561977.839 )
  
    deg,mn,sec=UNITS.radians2deg_mn_sec(rlong,angle_range='-180-180')
  
    assert deg == -78
    assert mn == 4
    assert np.sqrt((sec-16.92154)**2) < 0.00001
     
    deg,mn,sec=UNITS.radians2deg_mn_sec(rlat,angle_range='-180-180')
  
    assert deg == 45
    assert mn == 45
    assert np.sqrt((sec-48.48702)**2) < 0.00001
     
    assert np.sqrt((h-6367333.7313)**2)<0.0001

    