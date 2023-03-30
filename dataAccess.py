# -*- coding: utf-8 -*-
"""
Created on 7 Mar 2023

@author: steven.hosper
"""
# Functionality to access data for the HydrologicBaseModel
import lue.framework as lfr
import configuration
import numpy as np
import datetime
import uuid as uid
from osgeo import gdal
import configuration as config
from dataGen import generate as gen
import requests
from pcraster import aguila
import math
import pandas as pd

class get():
    def apiSession():
        session = requests.Session()
        session.headers = {'username': config.username,
                          'password': config.password,
                          'Content-Type': 'application/json',
                          }
        return session
    
    def apiTemporal(date: datetime.date, variable: str, session):
        """
        Summary:
            Get data with the use of an API from the Lizard server of Nelen & Schuurmans
            
        Input:
            date: the datetime.date of interest
            variable: the variable of interest
        
        Returns:
            data: the data from the tiff file in lue array format.
        """
        dates = f'{date}T11:00:00Z'
        
        # The different variables available with their corresponding Lizard rasters
        urls = {'precipitation'  : configuration.precipAPI, \
                'pot_evaporation': configuration.evapAPI,
                'dem'            : configuration.demAPI
                }
        
        url = f'https://demo.lizard.net/api/v4/rasters/{urls[variable]}/data/'
        # Make sure the variable input is correct
        if variable not in urls:
            raise Exception("The variable is not found. Available variables: 'precipitation', 'evaporation'.")
                
        # Format the API request
        request = {
                    "url": url,
                    "params":{
                       "cellsize": configuration.resolution,
                       "format": "geotiff",
                       "bbox": '3.330044347435246,50.723122219626666,7.278104205075899,53.80454395165623',          # bbox of the area of interest in coordinates
                       "width": configuration.arrayExtent,
                       "height": configuration.arrayExtent,
                       "srs": "EPSG:4326",
                       "target_srs": "EPSG:28992",
                       "start": dates,                                                                              # time
                    },}

        # Pull the request from the Lizard server
        pull = session.get(**request)
        pull.raise_for_status()
        
        # Assign a memory file to store the api content
        input_mem_file = f"/vsimem/input_{uid.uuid4()}.tif"
        output_mem_file = f"/vsimem/output_{uid.uuid4()}.tif"
        options = gdal.TranslateOptions(outputType= gdal.GDT_Float32)
        gdal.FileFromMemBuffer(input_mem_file, pull.content)
        gdal.Translate(output_mem_file, input_mem_file, options = options)
        
        data = lfr.from_gdal(output_mem_file, configuration.partitionShape)
        return data
    
    def localTemporal(path: str, date: datetime.date, variable: str) -> object:
        """
        Summary:
            Get data from the APIs to use in the model.
            Uses the format of "path/variable_date.tiff"
            
        Input:
            path: the directory where the files are located
            date: the datetime date
            variable: the variable name of interest
            
        Returns:
            data: the data from the tiff file to be used as an lue array
        """
        # Create the full path of the corresponding variable for the current
        variable_path = path + f'{variable}_{configuration.arrayExtent}_{date}.tiff'
        
        # Get the data from the .tif file stored in the path directory
        data = lfr.from_gdal(variable_path, configuration.partitionShape)
        return data
    
    def soil(soilType):
        # Assign K for de Wupsel
        Ks = gen.lue_one() * 0.05
        porosity = gen.lue_one() * 0.3
        scTable = pd.read_csv(config.path + f'/data/soil_conversion.csv')
        ID = scTable["ID"]
        KsValue = scTable["Ks"]
        for count, ID in enumerate(ID):
            Ks = lfr.where(soilType == ID, KsValue[count], Ks)
            
        return (Ks / 86400), porosity
    
    def landC(landUse):
        # Use the ID values given to the QGIS raster to determine which land-use types are assigned which values.
        land_c = gen.lue_one() * 0.8
        mannings = gen.lue_one() * 0.045
        soil_c = gen.lue_one() * 0.5
        for i in config.total:
            if i in configuration.concrete:
                land_c = lfr.where(landUse == i, 0.001, land_c)
                mannings = lfr.where(landUse == i, 0.015, mannings)
                soil_c = lfr.where(landUse == i, 0, soil_c)
                
            elif i in configuration.green:                             # Crops are given a multiplier of 1.2 as they also have pore structures \
                land_c = lfr.where(landUse == i, 1.1, land_c)            # but a different on can be assigned.
                mannings = lfr.where(landUse == i, 0.035, mannings)
                soil_c = lfr.where(landUse == i, 1, soil_c)
                
            elif i in configuration.water:                                              # Any type of water is assigned 1, as this should have the saturated hydraulic conductivity \
                land_c = lfr.where(landUse == i, 1.0, land_c)            # as precipitation value. --> However, these places probably have inflow from groundwater.
                mannings = lfr.where(landUse == i, 0.030, mannings)
                soil_c = lfr.where(landUse == i, 0, soil_c)
            
            elif i in configuration.compacted:
                land_c = lfr.where(landUse == i, 0.7, land_c)
                mannings = lfr.where(landUse == i, 0.025, mannings)
                soil_c = lfr.where(landUse == i, 0.5, soil_c)
                
            elif i in configuration.other_road:
                land_c = lfr.where(landUse == i, 0.3, land_c) 
                mannings = lfr.where(landUse == i, 0.015, mannings)
                soil_c = lfr.where(landUse == i, 0.2, soil_c)
                
        return land_c, mannings, soil_c
    
    def calculate_infiltration(Ks, landC):
        # Simplified version where the saturated hydraulic conductivity
        # is multiplied with a coefficient based on land-use
        return Ks * landC
    
    def precipitation(date, session):
        if config.v2:
            random = np.random.randint(0, config.dT)
            if random <= config.dT:
                precipitation = gen.lue_one()*0.05
                kernel = np.array(
                [
                    [1,1,1],
                    [1,1,1],
                    [1,1,1],
                ],
                dtype=np.uint8,
                )
                fraction_raincells = 0.1                                           # Determine the percentage of raincells in the array
                raincells = lfr.uniform(gen.lue_zero(), np.float32, 0, 1) <= fraction_raincells
                rainValue = lfr.focal_sum(raincells, kernel)
                for i in range(10):
                    precipitation = lfr.where(rainValue == i, precipitation * i, precipitation) 
                
            else:
                precipitation = gen.lue_zero()
        else:
            if configuration.includePrecipitation:
                if configuration.useAPI:
                    precipitation  = get.apiTemporal(date, 'precipitation', session)
                else:
                    precipitation  = get.localTemporal(f'{configuration.path}/data/generated/{configuration.arrayExtent}/', date, 'precipitation')
            else:
                precipitation = gen.lue_zero()
        return gen.lue_one() * 0 / (1000 * 3600) # Convert to meter / second rate
    
    def pot_evaporation(date, session):
        if config.v2:
            pot_evaporation = gen.lue_one() * 10 / 10 / 1000
        else:
            if configuration.includeEvaporation:
                if configuration.useAPI:
                    pot_evaporation   = get.apiTemporal(date, 'pot_evaporation', session)
                else:
                    pot_evaporation   = get.localTemporal(f'{configuration.path}/data/generated/{configuration.arrayExtent}/', date, 'potential_evaporation')
            else:
                pot_evaporation = gen.lue_zero()
        return pot_evaporation / 3600 # Convert to per second rate
    
    def infiltration(groundWaterHeight, Ks, land_c, porosity):
        if configuration.includeInfiltration:
            pot_infiltration = get.calculate_infiltration(Ks, land_c)
            pot_infiltration = lfr.where((groundWaterHeight - configuration.imperviousLayer)*porosity < pot_infiltration, \
                                          (configuration.imperviousLayer - groundWaterHeight)*porosity, pot_infiltration)
        else:
            pot_infiltration = gen.lue_zero()
        return pot_infiltration
    
    def EvaporationInfiltration(surfaceWaterHeight, pot_evaporation, pot_infiltration, e_ratio, i_ratio):
         # Check if the waterheight is more than the evaporation and infiltration combined
        ie = pot_evaporation + pot_infiltration
        
        # Calculate the actual evaporation and infiltration
        evaporation  = lfr.where(surfaceWaterHeight >= ie, pot_evaporation, surfaceWaterHeight*e_ratio)
        infiltration = lfr.where(surfaceWaterHeight >= ie, pot_infiltration, surfaceWaterHeight*i_ratio)
        
        return evaporation, infiltration
    
    def percolation(dem, groundWaterHeight, Ks):
        if configuration.includePercolation:
            part = (groundWaterHeight-0.5*dem)/dem
            part = lfr.where(part <= 0, gen.lue_zero(), part)
            percolation = lfr.where(dem < 0, gen.lue_zero(), \
                                    Ks * part)                                                # If the dem is negative, there is no percolation
            percolation = lfr.where(percolation < 0, gen.lue_zero(), percolation) / 100                       # If the percolation is negative, there is no percolation
        else:
            percolation = gen.lue_zero()
            part = gen.lue_zero()
        return percolation
    
    def ieRatio(pot_evaporation, pot_infiltration):
        if configuration.includeEvaporation and configuration.includeInfiltration:
            i_ratio = pot_infiltration / (pot_evaporation + pot_infiltration)
            e_ratio = pot_evaporation / (pot_evaporation + pot_infiltration)
        elif configuration.includeEvaporation:
            i_ratio = gen.lue_zero()
            e_ratio = gen.lue_one()
        elif configuration.includeInfiltration:
            i_ratio = gen.lue_one()
            e_ratio = gen.lue_zero()
        else:
            i_ratio = gen.lue_zero()
            e_ratio = gen.lue_zero()
        return i_ratio, e_ratio
    
    def interception(interceptionStorage, interceptionStorageMax, precipitation, landUse):
        # Im this simplified version the landUse type will be used for the LAI
        # This should be improved in a future version
        k = 0.5         # Constant taken from (Brolsma et al., 2010)
        
        f = gen.lue_zero()  # If not stated then there are no leafs
        
        for i in config.total:
            if i in config.green:
                f = lfr.where(landUse == i, math.exp(-k * 3), 0)
        
        interception = (gen.lue_one() - f) * precipitation        
        interception = lfr.where(interceptionStorageMax < interceptionStorage + interception, \
                                 interceptionStorageMax - interceptionStorage, interception)
        if config.includeInfiltration == False:
            interception = gen.lue_zero()            
        precipitation = precipitation - interception
        return interception, precipitation
    
    def interceptionStorageMax(landUse):
        LAI = gen.lue_zero()
        for i in config.total:
            if i in config.green:
                LAI = lfr.where(landUse == i, 3, LAI) # Currently everything within the 'green' category has an LAI of 3.
        interceptionStorageMax = (LAI / 3) / 1000
        return interceptionStorageMax
    
    def interceptionStorage(interceptionStorageOld, interception, evapotranspiration):
        interceptionStorage = interceptionStorageOld + interception - evapotranspiration
        return interceptionStorage
    
    def groundFlow():
        pass