#test
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import os

# Set which images are generated
plot_1d = False
make_movie = True

# Data file names
advDiffSimFile = "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/test_AdvDiff.nc"
advDiffTruthFile = "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/truth_AdvDiff.nc"
advDiffNLFile = "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/test_AdvDiffNL.nc"
testSolverFile = "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/testSolver.nc"

# Image indices
imageIdxs = [0, 9, 19, 29, 39]

# Read simple adv-diff data
#advDiffSim = netCDF4.Dataset(advDiffSimFile,'r')
#advDiffSimMasked = advDiffSim.variables['data_T']
#advDiffSimData = advDiffSimMasked[:,:,:]

# Read analytic solution to adv-diff solution
#advDiffTruth = netCDF4.Dataset(advDiffTruthFile,'r')
#advDiffTruthMasked = advDiffTruth.variables['data_T']
#advDiffTruthData = advDiffTruthMasked[:,:,:]

# Read simple nonlinear adv-diff data
#advDiffNL = netCDF4.Dataset(advDiffNLFile,'r')
#advDiffNLMasked = advDiffNL.variables['data_T']
#advDiffNLData = advDiffNLMasked[:,:,:]

testSolver = netCDF4.Dataset(testSolverFile,'r')
testSolverMaskedT = testSolver.variables['data_T']
testSolverMaskedU = testSolver.variables['data_u']
testSolverMaskedV = testSolver.variables['data_v']
testSolverDataT = testSolverMaskedT[:,:,:]
testSolverDataU = testSolverMaskedU[:,:,:]
testSolverDataV = testSolverMaskedV[:,:,:]

# Plot data
if plot_1d:
    #plt.figure(1)
    #advDiffPlot0 = plt.plot(advDiffSimData[imageIdxs[0],:,0])
    #advDiffPlot1 = plt.plot(advDiffSimData[imageIdxs[1],:,0])
    #advDiffPlot2 = plt.plot(advDiffSimData[imageIdxs[2],:,0])
    #advDiffPlot3 = plt.plot(advDiffSimData[imageIdxs[3],:,0])
    #advDiffPlot4 = plt.plot(advDiffSimData[imageIdxs[4],:,0])

    #advDiffTruthPlot0 = plt.plot(advDiffTruthData[imageIdxs[0],:,0])
    #advDiffTruthPlot1 = plt.plot(advDiffTruthData[imageIdxs[1],:,0])
    #advDiffTruthPlot2 = plt.plot(advDiffTruthData[imageIdxs[2],:,0])
    #advDiffTruthPlot3 = plt.plot(advDiffTruthData[imageIdxs[3],:,0])
    #advDiffTruthPlot4 = plt.plot(advDiffTruthData[imageIdxs[4],:,0])

    #err = plt.plot(advDiffTruthData[imageIdxs[4],:,0]-advDiffSimData[imageIdxs[4],:,0])

    #plt.figure(2)
    #advDiffNLPlot0 = plt.plot(advDiffNLData[imageIdxs[0],:,0])
    #advDiffNLPlot1 = plt.plot(advDiffNLData[imageIdxs[1],:,0])
    #advDiffNLPlot2 = plt.plot(advDiffNLData[imageIdxs[2],:,0])
    #advDiffNLPlot3 = plt.plot(advDiffNLData[imageIdxs[3],:,0])
    #advDiffNLPlot4 = plt.plot(advDiffNLData[imageIdxs[4],:,0])

    plt.figure(1)
    testSolverTPlot0 = plt.plot(testSolverDataU[imageIdxs[0],:,31])
    testSolverTPlot0 = plt.plot(testSolverDataU[imageIdxs[1],:,31])
    testSolverTPlot0 = plt.plot(testSolverDataU[imageIdxs[2],:,31])
    testSolverTPlot0 = plt.plot(testSolverDataU[imageIdxs[3],:,31])
    testSolverTPlot0 = plt.plot(testSolverDataU[imageIdxs[4],:,31])

    plt.show()

if make_movie:
    if not os.path.exists('images'):
        os.makedirs('images')
    for i in range(0,256):
        imageFileName ="./images/screen-" + str(i).zfill(4) + ".tif"
        fig = plt.figure()
        plt.imshow(testSolverDataT[i,:,:])
        fig.savefig(imageFileName)
        print(imageFileName)

    os.system("ffmpeg -r 30 -i ./images/screen-%04d.tif -crf 20 -pix_fmt yuv420p out.mp4")

#advDiffSim.close()
#advDiffTruth.close()
#advDiffNL.close()
testSolver.close()
