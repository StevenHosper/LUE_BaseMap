import pcraster as pcr
from osgeo import gdal
import time
from pcraster import aguila
import matplotlib.pyplot as plt
import sys, os

sys.path.insert(0, "F:/Projecten intern (2023)/Stage Steven Hosper/Model/v1/model/")
from configuration_v2 import Configuration

configuration = Configuration("F:/Projecten intern (2023)/Stage Steven Hosper/Model/v1/config/config.ini")

# Timer to add some measure of functionality to the program
start_time = time.time()
print("Starting")
# LDD CREATION
# Attempt with pcraster to make ldd
workdir = configuration.generalSettings["inputDir"] + configuration.generalSettings["scenario"]

pcr.setclone(int(configuration.modelSettings["arrayExtent"]), int(configuration.modelSettings["arrayExtent"]), float(configuration.modelSettings["resolution"]), 0, 0)
ds = gdal.Open(workdir + configuration.dataSettings["dem"])
print(ds)
raster = ds.GetRasterBand(1)
a = raster.ReadAsArray()

result = pcr.numpy2pcr(pcr.Scalar, a, -3.40282e+38)
print(result)
print("Starting the LDD create proces")
ldd = pcr.lddcreate(result, 9999999, 500000, 9999999, 9999999)
modDEM = pcr.lddcreatedem(result, 9999999, 500000, 9999999, 9999999)
print(ldd)
pcr.report(ldd, workdir + '/ldd_f64_pcr_shaped.tiff')
pcr.report(modDEM, workdir + '/dem_f64_pcr_shaped.tiff')

print("--- %s seconds ---" % (time.time() - start_time))