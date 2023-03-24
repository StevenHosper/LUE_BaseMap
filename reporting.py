# -*- coding: utf-8 -*-
"""
Created on Mon Mar 7 2023

@author: steven.hosper
"""
import lue.framework as lfr
import configuration as config

# Reporting for the HydrologicBaseModel
class report():
    def v2(date, time, variables, output_path):
        for variable, data in variables.items():
            lfr.to_gdal(data, output_path + f'{variable}_{date}_{time}_infil.tiff')
        return 0
    
    
    def static(date, variables, output_path):
        for variable, data in variables.items():
            lfr.to_gdal(data, output_path + f'{variable}_{date}.tiff')
        return 0
    
    def dynamic(date, second, variables, fluxes, output_path):
        if config.v2:
            rest = second % 15
            if rest == 0:
                for variable, data in variables.items():
                    lfr.to_gdal(data, output_path + f'{variable}_{date}_{second}.tiff')
                for flux, data in fluxes.items():
                    lfr.to_gdal(data, output_path + f'{flux}_{date}.tiff')
            else:
                pass
        else:
            for variable, data in variables.items():
                lfr.to_gdal(data, output_path + f'{variable}_{date}.tiff')
            for flux, data in fluxes.items():
                lfr.to_gdal(data, output_path + f'{flux}_{date}.tiff')
        return 0