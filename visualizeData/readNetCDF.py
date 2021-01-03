#test
import numpy as np
import netCDF4
import matplotlib.pyplot as plt

test = netCDF4.Dataset("/Users/rmoll/Desktop/test_AdvDiff_xy.nc",'r')
dataMasked = test.variables['data']
data = dataMasked[:,:,:]
dataCol = data[0,:,0]

plt.figure(1)
twoDPlot = plt.imshow(data[0,:,:])
plt.figure(2)
oneDPlot = plt.plot(dataCol)
plt.show()
#plt.imshow(data[1,:,:])
#plt.plot(dataCol)
