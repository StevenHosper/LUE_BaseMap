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
from StandardArraysLUE import StandardArraysLUE

# Reporting for the HydrologicBaseModel
class Report:
    def __init__(self, configuration):
        self.standard_LUE   = StandardArraysLUE(configuration)
        self.timestep = configuration.modelSettings['timestep']
        self.output_dir = configuration.generalSettings['outputDir'] + configuration.generalSettings["scenario"]
        self.input_dir   = configuration.generalSettings['inputDir'] + configuration.generalSettings['scenario'] 
        
    def dynamic(self, date, variables):
        dateTime = date.strftime("%Y-%m-%d-%H%M")
        for variable, data in variables.items():
            lfr.to_gdal(data, self.output_dir + "/{}_{}_{}.tiff".format(self.timestep,
                                                                         variable,
                                                                         dateTime
                                                                         ))
        return 0
     
    def balance_report(self, configuration):
        start_date = self.string_to_datetime(configuration.modelSettings['startDate'], seperator= ", ")
        end_date   = self.string_to_datetime(configuration.modelSettings['endDate'], seperator= ", ")
        start_date_txt = start_date.strftime("%d/%m/%Y %H:%M")
        end_date_txt   = end_date.strftime("%d/%m/%Y %H:%M")
        
        
        if configuration.generalSettings['includePrecipitation'] == "True":
            p = pd.read_csv(configuration.generalSettings["inputDir"] + configuration.dataSettings["precipitationData"], names=["date", "p"])
            start_idx= p.index[p["date"] == start_date_txt].tolist()[0]
            end_idx  = p.index[p["date"] == end_date_txt].tolist()[0]
            mean_precipitation = p.p[start_idx:(end_idx-1)].mean()
        else:
            mean_precipitation = 0
            start_idx = 0
            end_idx   = (end_date - start_date).total_seconds() / 300
        
        if configuration.generalSettings['includeEvapotranspiration'] == "True":
            e = pd.read_csv(configuration.generalSettings["inputDir"] + configuration.dataSettings["evapotranspirationData"], names=["date", "e"])
            start_idx= e.index[e["date"] == start_date_txt].tolist()[0]
            end_idx  = e.index[e["date"] == end_date_txt].tolist()[0]
            mean_evapotranspiration = e.e[start_idx:(end_idx-1)].mean()
        else:
            mean_evapotranspiration = 0
        
        end_date = end_date - datetime.timedelta(minutes=1)
        
        # Load initial files, if they cannot be loaded, assume zero (same as model does)
        try:
            ini_sur_stor = self.tiff_to_np_sum(self.input_dir + configuration.dataSettings["iniWaterHeight"])
        except:
            self.standard_LUE.zero()
        
        try:
            ini_gro_stor = self.input_dir + configuration.dataSettings['iniGroundWaterStorage']
        except:
            self.standard_LUE.zero()
        
        try:
            ini_int_stor = self.tiff_to_np_sum(self.input_dir + configuration.dataSettings['iniInterceptionStorage'])
        except:
            self.standard_LUE.zero()
        
        end_sur_stor = self.tiff_to_np_sum(self.output_dir + "/{}_height_{}.tiff".format(int(configuration.modelSettings["timestep"]), end_date.strftime("%Y-%m-%d-%H%M")))
        end_gro_stor = self.output_dir + "/{}_gw_s_{}.tiff".format(int(configuration.modelSettings["timestep"]), end_date.strftime("%Y-%m-%d-%H%M"))
        end_int_stor = self.tiff_to_np_sum(self.output_dir + "/{}_int_s_{}.tiff".format(int(configuration.modelSettings["timestep"]), end_date.strftime("%Y-%m-%d-%H%M")))
        
        resolution = (int(configuration.modelSettings["resolution"]))
        cell_area   = resolution ** 2
        
        del_sur_stor = (end_sur_stor - ini_sur_stor) * resolution
        del_gro_stor = self.tiff_to_np_sum_difference(end_gro_stor, ini_gro_stor) * float(configuration.modelSettings["porosity"])
        del_int_stor = (end_int_stor - ini_int_stor)
        net_balance = del_sur_stor + del_int_stor + del_gro_stor
        precipitation       = (((end_idx - start_idx) / 12) * mean_precipitation) / 1000 * cell_area * (int(configuration.modelSettings["arrayExtent"]) ** 2) * (float(configuration.modelSettings["validCellsPercentage"]))/100
        evapotranspiration  = (((end_idx - start_idx) / 12) * mean_evapotranspiration) / 1000 * cell_area * (int(configuration.modelSettings["arrayExtent"]) ** 2) * (float(configuration.modelSettings["validCellsPercentage"]))/100
        atmospheric_balance = precipitation - evapotranspiration
        
        print("total precipitation:               ", precipitation, "m3")
        print("total evapotranspiration:          ", evapotranspiration, "m3 \n")
        
        print("iniTotalSurfaceHeight:             ", ini_sur_stor * resolution, "m3")
        print("endTotalSurfaceHeight:             ", end_sur_stor * resolution, "m3")
        print("delta Surface Storage:             ", del_sur_stor, "m3 \n")
        
        print("delta GroundWater Storage:         ", del_gro_stor, "m3 \n")
        
        print("iniTotalInterceptionStorage:       ", ini_int_stor, "m3")
        print("endTotalInterceptionStorage:       ", end_int_stor, "m3")
        print("delta Interception Storage:        ", del_int_stor, " m3 \n")
        
        print("atmospheric balance:               ", atmospheric_balance, "m3")
        print("storage balance:                   ", net_balance, "m3")
        print("del atmos and storage:             ", atmospheric_balance - net_balance, "m3 \n")
        
        print("time simulated:                    ", (end_idx-start_idx)*5*60, "s")
        print("waterbalance change in the system: ", (net_balance - atmospheric_balance)/((end_idx-start_idx)*5*60), "m3/s")
        
        ofdf            = pd.read_csv(self.output_dir + "/maximumDischarge.csv", sep=";", names=["timestep", "outflow"])
        average_outflow = ofdf["outflow"].mean() * -1
        print("measured loss to outflow:   ", average_outflow, "m3/s \n")
        
        return 0
    
    def tiff_to_np_sum(self, file):
        try:
            data = gdal.Open(file)
            img = data.GetRasterBand(1)
            raster = img.ReadAsArray()
            np_array_sum = np.nansum(raster)
        except:
            np_array_sum = 0
        return np_array_sum
    
    def string_to_datetime(self, date_string: str, seperator: str):
        date_int_list = list(map(int, date_string.split(seperator)))
        datetime_date = datetime.datetime(date_int_list[0],
                                          date_int_list[1],
                                          date_int_list[2],
                                          date_int_list[3],
                                          date_int_list[4],
                                          date_int_list[5])
        return datetime_date
        
    
    def tiff_to_np_sum_difference(self, file1, file2):
        try:
            data1 = gdal.Open(file1)
            img1 = data1.GetRasterBand(1)
            raster1 = img1.ReadAsArray()
            data2 = gdal.Open(file2)
            img2 = data2.GetRasterBand(1)
            raster2 = img2.ReadAsArray()
            np_array_sum = np.nansum(raster1 - raster2)
        except:
            np_array_sum = 0
        return np_array_sum