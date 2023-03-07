# -*- coding: utf-8 -*-
"""
Created on Mon Mar 7 2023

@author: steven.hosper
"""
import lue.framework as lfr

# Reporting for the HydrologicBaseModel
class report():
    def static(date, variables, output_path):
        for variable, data in variables.items():
            lfr.to_gdal(data, output_path + f'{variable}_{date}.tiff')
        return 0
    
    def dynamic(date, variables, fluxes, output_path):
        for variable, data in variables.items():
            lfr.to_gdal(data, output_path + f'{variable}_{date}.tiff')
        for flux, data in fluxes.items():
            lfr.to_gdal(data, output_path + f'{flux}_{date}.tiff')
        return 0