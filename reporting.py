# -*- coding: utf-8 -*-
"""
Created on Mon Mar 7 2023

@author: steven.hosper
"""
import lue.framework as lfr
import configuration as config

# Reporting for the HydrologicBaseModel
class report():
    def dynamic(date, time, variables, output_path):
        for variable, data in variables.items():
            lfr.to_gdal(data, output_path + f'/{config.timestep}_{variable}_{date}_{time}.tiff')
        return 0