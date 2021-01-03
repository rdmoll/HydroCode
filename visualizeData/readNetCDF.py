#test
import numpy as np
import netCDF4
import matplotlib.pyplot as plt

test = netCDF4.Dataset("/Users/rmoll/Desktop/test_AdvDiff_xy.nc",'r')
dataMasked = test.variables['data']
data = dataMasked[:,:,:]
dataCol0 = data[0,:,0]
dataCol1 = data[1,:,0]
dataCol2 = data[2,:,0]
dataCol3 = data[3,:,0]

#plt.figure(1)
#twoDPlot = plt.imshow(data[0,:,:])
#plt.figure(2)
oneDPlot0 = plt.plot(dataCol0)
oneDPlot1 = plt.plot(dataCol1)
oneDPlot2 = plt.plot(dataCol2)
oneDPlot3 = plt.plot(dataCol3)
plt.show()
#plt.imshow(data[1,:,:])
#plt.plot(dataCol)
