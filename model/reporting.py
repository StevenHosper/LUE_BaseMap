# -*- coding: utf-8 -*-
"""
Created on Mon Mar 7 2023

@author: steven.hosper
"""
import lue.framework as lfr
import datetime
import numpy as np
from osgeo import gdal
import pandas as pd

# Reporting for the HydrologicBaseModel
class Report:
    def __init__(self, configuration):
        self.timestep = configuration.modelSettings['timestep']
        self.output_path = configuration.generalSettings['outputDir']
        
    def dynamic(self, date, variables):
        dateTime = date.strftime("%Y-%m-%d-%H%M")
        for variable, data in variables.items():
            lfr.to_gdal(data, self.output_path + "/{}_{}_{}.tiff".format(self.timestep,
                                                                         variable,
                                                                         dateTime
                                                                         ))
        return 0
     
    def balance_report(self, configuration):
        dir = configuration.generalSettings["outputDir"] + configuration.generalSettings["scenario"]
        
        start_date = self.string2datetime(configuration.modelSettings['startDate'], seperator= ", ")
        end_date   = self.string2datetime(configuration.modelSettings['endDate'], seperator= ", ")
        start_date_txt = start_date.strftime("%d/%m/%Y %H:%M")
        end_date_txt   = end_date.strftime("%d/%m/%Y %H:%M")
        
        if configuration.generalSettings['includePrecipitation'] == "True":
            p = pd.read_csv(configuration.generalSettings["inputDir"] + configuration.dataSettings["precipitationData"], names=["date", "p"])
            start_IDX= p.index[p["date"] == start_date_txt].tolist()[0]
            end_IDX  = p.index[p["date"] == end_date_txt].tolist()[0]
            mean_Precipitation = p.p[start_IDX:(end_IDX-1)].mean()
        else:
            mean_Precipitation = 0
        
        if configuration.generalSettings['includeEvapotranspiration'] == "True":
            e = pd.read_csv(configuration.generalSettings["inputDir"] + configuration.dataSettings["evapotranspirationData"], names=["date", "e"])
            start_IDX= e.index[e["date"] == start_date_txt].tolist()[0]
            end_IDX  = e.index[e["date"] == end_date_txt].tolist()[0]
            mean_Evapotranspiration = e.e[start_IDX:(end_IDX-1)].mean()
        else:
            mean_Evapotranspiration = 0
        
        end_date = end_date - datetime.timedelta(minutes=1)
        
        ini_sur_stor = self.tiff2npsum(dir + "/{}_height_{}.tiff".format(int(configuration.modelSettings["timestep"]), start_date.strftime("%Y-%m-%d-%H%M")))
        iniGroStor = dir + "/{}_Sgw_{}.tiff".format(int(configuration.modelSettings["timestep"]), start_date.strftime("%Y-%m-%d-%H%M"))
        ini_int_stor = self.tiff2npsum(dir + "/{}_interceptionStorage_{}.tiff".format(int(configuration.modelSettings["timestep"]), start_date.strftime("%Y-%m-%d-%H%M")))
        end_sur_stor = self.tiff2npsum(dir + "/{}_height_{}.tiff".format(int(configuration.modelSettings["timestep"]), end_date.strftime("%Y-%m-%d-%H%M")))
        endGroStor = dir + "/{}_Sgw_{}.tiff".format(int(configuration.modelSettings["timestep"]), end_date.strftime("%Y-%m-%d-%H%M"))
        end_int_stor = self.tiff2npsum(dir + "/{}_interceptionStorage_{}.tiff".format(int(configuration.modelSettings["timestep"]), end_date.strftime("%Y-%m-%d-%H%M")))
        
        resolution = (int(configuration.modelSettings["resolution"]))
        cell_area   = resolution ** 2
        
        del_sur_stor = (end_sur_stor - ini_sur_stor) * resolution
        del_gro_stor = self.tiff2npsumdifference(endGroStor, iniGroStor) * float(configuration.modelSettings["porosity"])
        del_int_stor = (end_int_stor - ini_int_stor)
        net_balance = del_sur_stor + del_int_stor + del_gro_stor
        precipitation       = (((end_IDX - start_IDX) / 12) * mean_Precipitation) / 1000 * cell_area * (int(configuration.modelSettings["arrayExtent"]) ** 2) * (float(configuration.modelSettings["validCellsPercentage"]))/100
        evapotranspiration  = (((end_IDX - start_IDX) / 12) * mean_Evapotranspiration) / 1000 * cell_area * (int(configuration.modelSettings["arrayExtent"]) ** 2) * (float(configuration.modelSettings["validCellsPercentage"]))/100
        atmospheric_balance = precipitation - evapotranspiration
        
        print("Total Precipitation:         ", precipitation, "m3")
        print("Total Evapotranspiration:    ", evapotranspiration, "m3")
        print("Atmospheric balance:         ", atmospheric_balance, "m3 \n")
        
        print("iniTotalSurfaceHeight:       ", ini_sur_stor * resolution, "m3")
        print("endTotalSurfaceHeight:       ", end_sur_stor * resolution, "m3")
        print("delta Surface Storage:       ", del_sur_stor, "m3 \n")
        
        print("delta GroundWater Storage:   ", del_gro_stor, "m3 \n")
        
        print("iniTotalInterceptionStorage: ", ini_int_stor, "m3")
        print("endTotalInterceptionStorage: ", end_int_stor, "m3")
        print("delta Interception Storage:  ", del_int_stor, " m3 \n")
        
        print("storage balance:             ", net_balance, "m3 \n")
        
        print("del atmos and storage:       ", atmospheric_balance - net_balance, "m3")
        print("time simulated:              ", (end_IDX-start_IDX)*5*60, "s")
        print("expected loss to outflow:    ", (net_balance - atmospheric_balance)/((end_IDX-start_IDX)*5*60), "m3/s \n")
        
        return 0
    
    def tiff2npsum(self, file):
        data = gdal.Open(file)
        img = data.GetRasterBand(1)
        raster = img.ReadAsArray()
        npArraySum = np.nansum(raster)
        return npArraySum
    
    def string2datetime(self, date_string: str, seperator: str):
        date_int_list = list(map(int, date_string.split(seperator)))
        datetime_date = datetime.datetime(date_int_list[0],
                                          date_int_list[1],
                                          date_int_list[2],
                                          date_int_list[3],
                                          date_int_list[4],
                                          date_int_list[5])
        return datetime_date
        
    
    def tiff2npsumdifference(self, file1, file2):
        data1 = gdal.Open(file1)
        img1 = data1.GetRasterBand(1)
        raster1 = img1.ReadAsArray()
        data2 = gdal.Open(file2)
        img2 = data2.GetRasterBand(1)
        raster2 = img2.ReadAsArray()
        np_array_sum = np.nansum(raster1 - raster2)
        return np_array_sum