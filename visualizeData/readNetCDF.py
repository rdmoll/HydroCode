#test
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import os
from PIL import Image

plot_1d = True
make_movie = False

advDiffFile_test = "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/test_AdvDiff.nc"
advDiffFile_truth = "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/truth_AdvDiff.nc"
advDiffNLFile = "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/test_AdvDiffNL.nc"

# Read Data
test = netCDF4.Dataset(advDiffFile_test,'r')
dataMasked = test.variables['data']
data = dataMasked[:,:,:]
dataColTest0 = data[0,:,0]
dataColTest1 = data[9,:,0]
dataColTest2 = data[19,:,0]
dataColTest3 = data[29,:,0]
dataColTest4 = data[39,:,0]

truth = netCDF4.Dataset(advDiffFile_truth,'r')
dataMaskedTruth = truth.variables['data']
dataTruth = dataMaskedTruth[:,:,:]
dataColTruth0 = dataTruth[0,:,0]
dataColTruth1 = dataTruth[9,:,0]
dataColTruth2 = dataTruth[19,:,0]
dataColTruth3 = dataTruth[29,:,0]
dataColTruth4 = dataTruth[39,:,0]

testNL = netCDF4.Dataset(advDiffNLFile,'r')
dataMaskedNL = testNL.variables['data']
dataNL = dataMaskedNL[:,:,:]
dataColNL0 = dataNL[0,:,0]
dataColNL1 = dataNL[9,:,0]
dataColNL2 = dataNL[19,:,0]
dataColNL3 = dataNL[29,:,0]
dataColNL4 = dataNL[39,:,0]

# Plot data
if plot_1d:
    plt.figure(1)
    oneDPlot0 = plt.plot(dataColTest0)
    #oneDPlot1 = plt.plot(dataColTest1)
    #oneDPlot2 = plt.plot(dataColTest2)
    #oneDPlot3 = plt.plot(dataColTest3)
    oneDPlot4 = plt.plot(dataColTest4)
    #oneDPlotTruth0 = plt.plot(dataColTruth0)
    #oneDPlotTruth1 = plt.plot(dataColTruth1)
    #oneDPlotTruth2 = plt.plot(dataColTruth2)
    #oneDPlotTruth3 = plt.plot(dataColTruth3)
    oneDPlotTruth4 = plt.plot(dataColTruth4)
    #err = plt.plot(dataColTruth4-dataCol1)

    plt.figure(2)
    oneDNLPlot0 = plt.plot(dataColNL0)
    oneDNLPlot1 = plt.plot(dataColNL1)
    oneDNLPlot2 = plt.plot(dataColNL2)
    oneDNLPlot3 = plt.plot(dataColNL3)
    oneDNLPlot4 = plt.plot(dataColNL4)

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
