#test
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import os
from PIL import Image

# Read Data
test = netCDF4.Dataset("/Users/rmoll/Desktop/test_AdvDiff_xy.nc",'r')
dataMasked = test.variables['data']
data = dataMasked[:,:,:]
dataCol0 = data[0,:,0]
dataCol1 = data[9,:,0]
dataCol2 = data[19,:,0]
dataCol3 = data[29,:,0]
dataCol4 = data[39,:,0]

# Plot data
plt.figure(1)
oneDPlot0 = plt.plot(dataCol0)
oneDPlot1 = plt.plot(dataCol1)
oneDPlot2 = plt.plot(dataCol2)
oneDPlot3 = plt.plot(dataCol3)
oneDPlot4 = plt.plot(dataCol4)
plt.show()

#if not os.path.exists('images'):
#    os.makedirs('images')
#for i in range(0,40):
#    imageFileName ="./images/screen-" + str(i).zfill(4) + ".tif"
#    fig = plt.figure()
#    plt.imshow(data[i,:,:])
#    fig.savefig(imageFileName)
#    print(imageFileName)
#
#os.system("ffmpeg -r 30 -i ./images/screen-%04d.tif -crf 20 -pix_fmt yuv420p out.mp4")
