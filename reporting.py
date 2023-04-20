# -*- coding: utf-8 -*-
"""
Created on Mon Mar 7 2023

@author: steven.hosper
"""
import lue.framework as lfr
import configuration as config
import datetime

# Reporting for the HydrologicBaseModel
class report():
    def dynamic(date, timestep, variables, output_path):
        dateTime = date.strftime("%Y-%m-%d-%H%M")
        for variable, data in variables.items():
            lfr.to_gdal(data, output_path + "/{}_{}_{}.tiff".format(int(timestep), variable, dateTime))
        return 0