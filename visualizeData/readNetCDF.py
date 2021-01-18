#test
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import os
from PIL import Image

plot_1d = True
make_movie = False

# Read Data
test = netCDF4.Dataset("/Users/rmoll/Desktop/test_AdvDiff_xy.nc",'r')
dataMasked = test.variables['data']
data = dataMasked[:,:,:]
dataCol0 = data[0,:,0]
dataCol1 = data[9,:,0]
dataCol2 = data[19,:,0]
dataCol3 = data[29,:,0]
dataCol4 = data[39,:,0]

truth = netCDF4.Dataset("/Users/rmoll/Desktop/truth_AdvDiff_xy.nc",'r')
dataMaskedTruth = truth.variables['data']
dataTruth = dataMaskedTruth[:,:,:]
dataColTruth0 = dataTruth[0,:,0]
dataColTruth1 = dataTruth[9,:,0]
dataColTruth2 = dataTruth[19,:,0]
dataColTruth3 = dataTruth[29,:,0]
dataColTruth4 = dataTruth[39,:,0]

# Plot data
if plot_1d:
    plt.figure(1)
    #oneDPlot0 = plt.plot(dataCol0)
    #oneDPlot1 = plt.plot(dataCol1)
    #oneDPlot2 = plt.plot(dataCol2)
    #oneDPlot3 = plt.plot(dataCol3)
    #oneDPlot4 = plt.plot(dataCol4)
    #plt.show()

    #plt.figure(1)
    #oneDPlotTruth0 = plt.plot(dataColTruth0)
    #oneDPlotTruth1 = plt.plot(dataColTruth1)
    #oneDPlotTruth2 = plt.plot(dataColTruth2)
    #oneDPlotTruth3 = plt.plot(dataColTruth3)
    #oneDPlotTruth4 = plt.plot(dataColTruth4)
    err = plt.plot(dataColTruth4-dataCol1)
    #err = plt.plot(dataColTruth4-dataCol4)
    plt.show()

if make_movie:
    if not os.path.exists('images'):
        os.makedirs('images')
    for i in range(0,40):
        imageFileName ="./images/screen-" + str(i).zfill(4) + ".tif"
        fig = plt.figure()
        plt.imshow(data[i,:,:])
        fig.savefig(imageFileName)
        print(imageFileName)

    os.system("ffmpeg -r 30 -i ./images/screen-%04d.tif -crf 20 -pix_fmt yuv420p out.mp4")

test.close()
truth.close()
