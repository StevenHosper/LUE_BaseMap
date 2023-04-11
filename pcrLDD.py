import configuration as config
import pcraster as pcr
from osgeo import gdal
import time
from pcraster import aguila
import matplotlib.pyplot as plt

# Timer to add some measure of functionality to the program
start_time = time.time()
print("Starting")
# LDD CREATION
# Attempt with pcraster to make ldd
pcr.setclone(config.arrayExtent, config.arrayExtent, config.resolution, 0, 0)
ds = gdal.Open(config.path + f'/data/{config.scenario}/dem.tiff')
print(ds)
raster = ds.GetRasterBand(1)
a = raster.ReadAsArray()

result = pcr.numpy2pcr(pcr.Scalar, a, 3.40282e+38)
print(result)
print("Starting the LDD create proces")
ldd = pcr.lddcreate(result, 9999999, 50000, 9999999, 9999999)
print(ldd)
pcr.report(ldd, config.path + f'/data/{config.scenario}/ldd_pcr_shaped.tiff')

print("--- %s seconds ---" % (time.time() - start_time))