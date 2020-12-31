#test
import numpy as np
import netCDF4

test = netCDF4.Dataset("test_class_xy.nc",'r')
data = test['data'][:][:][:]
#data = test.variables['data'][:][:][:]
print(data)
