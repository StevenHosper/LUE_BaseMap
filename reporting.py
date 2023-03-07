# -*- coding: utf-8 -*-
"""
Created on Mon Mar 7 2023

@author: steven.hosper
"""
import lue.framework as lfr

# Reporting for the HydrologicBaseModel
class report():
    def static(date, variables, output_path):
        for variable in variables:
            lfr.to_gdal(variable, output_path + f'{variable}_{date}.tiff')
        return 0
    
    def dynamic(date, variables, fluxes, output_path):
        for variable in variables:
            lfr.to_gdal(variable, output_path + f'{variable}_{date}.tiff')
        for flux in fluxes:
            lfr.to_gdal(flux, output_path + f'{variable}_{date}.tiff')
        return 0