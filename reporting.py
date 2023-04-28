# -*- coding: utf-8 -*-
"""
Created on Mon Mar 7 2023

@author: steven.hosper
"""
import lue.framework as lfr
import configuration as config
import datetime
import numpy as np
from osgeo import gdal
import pandas as pd

# Reporting for the HydrologicBaseModel
class report():
    def dynamic(date, timestep, variables, output_path):
        dateTime = date.strftime("%Y-%m-%d-%H%M")
        for variable, data in variables.items():
            lfr.to_gdal(data, output_path + "/{}_{}_{}.tiff".format(int(timestep), variable, dateTime))
        return 0
     
    def balanceReport(configuration):
        dir = configuration.generalSettings["outputDir"] + configuration.generalSettings["scenario"]
        
        p = pd.read_csv(configuration.generalSettings["inputDir"] + configuration.dataSettings["precipitationData"], names=["date", "p"])
        e = pd.read_csv(configuration.generalSettings["inputDir"] + configuration.dataSettings["evapotranspirationData"], names=["date", "e"])
        startIDX = 3
        endIDX   = 123
        meanPrecipitation = p.p[startIDX:endIDX].mean()
        meanEvapotranspiration = e.e[startIDX:endIDX].mean()
        print(meanPrecipitation)
        
        sD = list(map(int, configuration.modelSettings['startDate'].split(", ")))
        eD   = list(map(int, configuration.modelSettings['endDate'].split(", ")))
        startDate   = datetime.datetime(sD[0], sD[1], sD[2], sD[3], sD[4], sD[5])
        endDate     = datetime.datetime(eD[0], eD[1], eD[2], eD[3], eD[4], eD[5]) - datetime.timedelta(minutes=1)
        
        
        iniSurStor = report.tiff2npsum(dir + "/{}_height_{}.tiff".format(int(configuration.modelSettings["timestep"]), startDate.strftime("%Y-%m-%d-%H%M")))
        iniGroStor = dir + "/{}_Sgw_{}.tiff".format(int(configuration.modelSettings["timestep"]), startDate.strftime("%Y-%m-%d-%H%M"))
        iniIntStor = report.tiff2npsum(dir + "/{}_interceptionStorage_{}.tiff".format(int(configuration.modelSettings["timestep"]), startDate.strftime("%Y-%m-%d-%H%M")))
        endSurStor = report.tiff2npsum(dir + "/{}_height_{}.tiff".format(int(configuration.modelSettings["timestep"]), endDate.strftime("%Y-%m-%d-%H%M")))
        endGroStor = dir + "/{}_Sgw_{}.tiff".format(int(configuration.modelSettings["timestep"]), endDate.strftime("%Y-%m-%d-%H%M"))
        endIntStor = report.tiff2npsum(dir + "/{}_interceptionStorage_{}.tiff".format(int(configuration.modelSettings["timestep"]), endDate.strftime("%Y-%m-%d-%H%M")))
        
        resolution = (int(configuration.modelSettings["resolution"]))
        cellArea   = resolution ** 2
        
        delSurStor = (endSurStor - iniSurStor) * resolution
        delGroStor = report.tiff2npsumdifference(endGroStor, iniGroStor) * float(configuration.modelSettings["porosity"])
        delIntStor = (endIntStor - iniIntStor)
        netBalance = delSurStor + delIntStor + delGroStor
        precipitation       = (((endIDX - startIDX) / 12) * meanPrecipitation) / 1000 * cellArea * (int(configuration.modelSettings["arrayExtent"]) ** 2) * (float(configuration.modelSettings["validCellsPercentage"]))/100
        evapotranspiration  = (((endIDX - startIDX) / 12) * meanEvapotranspiration) / 1000 * cellArea * (int(configuration.modelSettings["arrayExtent"]) ** 2) * (float(configuration.modelSettings["validCellsPercentage"]))/100
        atmosphericBalance = precipitation - evapotranspiration
        
        print("Total Precipitation:         ", precipitation, "m3")
        print("Total Evapotranspiration:    ", evapotranspiration, "m3")
        print("Atmospheric balance:         ", atmosphericBalance, "m3 \n")
        
        print("iniTotalSurfaceHeight:       ", iniSurStor * resolution, "m3")
        print("endTotalSurfaceHeight:       ", endSurStor * resolution, "m3")
        print("delta Surface Storage:       ", delSurStor, "m3 \n")
        
        print("delta GroundWater Storage:   ", delGroStor, "m3 \n")
        
        print("iniTotalInterceptionStorage: ", iniIntStor, "m3")
        print("endTotalInterceptionStorage: ", endIntStor, "m3")
        print("delta Interception Storage:  ", delIntStor, " m3 \n")
        
        print("storage balance:             ", netBalance, "m3 \n")
        
        print("del atmos and storage:       ", atmosphericBalance - netBalance, "m3")
        print("time simulated:              ", (endIDX-startIDX)*5*60, "s")
        print("expected loss to outflow:    ", (netBalance - atmosphericBalance)/((endIDX-startIDX)*5*60), "m3/s \n")
        
        return 0
    
    def tiff2npsum(file):
        data = gdal.Open(file)
        img = data.GetRasterBand(1)
        raster = img.ReadAsArray()
        npArraySum = np.nansum(raster)
        return npArraySum
    
    def tiff2npsumdifference(file1, file2):
        data1 = gdal.Open(file1)
        img1 = data1.GetRasterBand(1)
        raster1 = img1.ReadAsArray()
        data2 = gdal.Open(file2)
        img2 = data2.GetRasterBand(1)
        raster2 = img2.ReadAsArray()
        npArraySum = np.nansum(raster1 - raster2)
        return npArraySum