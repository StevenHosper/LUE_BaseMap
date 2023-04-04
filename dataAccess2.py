# -*- coding: utf-8 -*-
"""
Created on 30 Mar 2023

@author: steven.hosper
"""
# Functionality to access data for the HydrologicBaseModel
import lue.framework as lfr
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
        urls = {'precipitation'  : config.precipAPI, \
                'pot_evaporation': config.evapAPI,
                'dem'            : config.demAPI
                }
        
        url = f'https://demo.lizard.net/api/v4/rasters/{urls[variable]}/data/'
        # Make sure the variable input is correct
        if variable not in urls:
            raise Exception("The variable is not found. Available variables: 'precipitation', 'evaporation'.")
                
        # Format the API request
        request = {
                    "url": url,
                    "params":{
                       "cellsize": config.resolution,
                       "format": "geotiff",
                       "bbox": '3.330044347435246,50.723122219626666,7.278104205075899,53.80454395165623',          # bbox of the area of interest in coordinates
                       "width": config.arrayExtent,
                       "height": config.arrayExtent,
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
        
        data = lfr.from_gdal(output_mem_file, config.partitionShape)
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
        variable_path = path + f'{variable}_{config.arrayExtent}_{date}.tiff'
        
        # Get the data from the .tif file stored in the path directory
        data = lfr.from_gdal(variable_path, config.partitionShape)
        return data
    
    def soil_csv(dataFile, soilType):
        # Assign K for de Wupsel
        Ks = gen.lue_one() * 0.05
        porosity = gen.lue_one() * 0.35
        scTable = pd.read_csv(dataFile)
        scTable = pd.read_csv(config.path + f'/data/soil_conversion.csv')
        ID = scTable["ID"]
        KsValue = scTable["Ks"]
        for count, ID in enumerate(ID):
            Ks = lfr.where(soilType == ID, KsValue[count], Ks)
            
        return (Ks / 86400), porosity
    
    def landCharacteristics_csv(dataFile, soilType):
        # Use the ID values given to the QGIS raster to determine which land-use types are assigned which values.
        # Standard values
        dummy                  = gen.lue_one()
        permeability           = dummy * 0.8
        mannings               = dummy * 0.045
        interceptionStorageMax = dummy * 0.001
        throughfallFraction    = dummy * 0.90
        LAI                     = dummy * 0.5
        cropFactor              = dummy * 1.0
        
        data = pd.read_csv(dataFile)
        ID = data["Code"]
        manningsFriction            = data["Friction"]
        permeabilityValue           = data["Permeability"]
        interceptionStorageMaxValue = data["Interception"]
        LAIValue                    = data["LAI"]
        throughfallFractionValue    = data["f"]
        cropFactorValue             = data["Crop_type"]
        
        for count, ID in enumerate(ID):
            mannings                = lfr.where(soilType == ID, manningsFriction[count], mannings)
            permeability            = lfr.where(soilType == ID, permeabilityValue[count], permeability)
            interceptionStorageMax  = lfr.where(soilType == ID, interceptionStorageMaxValue[count], interceptionStorageMax)
            LAI                     = lfr.where(soilType == ID, LAIValue[count], LAI)
            throughfallFraction     = lfr.where(soilType == ID, throughfallFractionValue[count], throughfallFraction)
            cropFactor              = lfr.where(soilType == ID, cropFactorValue[count], cropFactor)
                
        return mannings, permeability, interceptionStorageMax, LAI, throughfallFraction, cropFactor
    
    def landCharacteristics_old(landUse):
        # Use the ID values given to the QGIS raster to determine which land-use types are assigned which values.
        dummy                   = gen.lue_one()
        permeability            = dummy * 0.8
        mannings                = dummy * 0.045
        interceptionStorageMax  = dummy * 0.001
        throughfallFraction     = dummy * 0.90
        LAI                     = dummy * 0.5
        cropFactor              = dummy * 1.0
        
        for i in config.total:
            if i in config.concrete:
                permeability            = lfr.where(landUse == i, 0.001, permeability)
                mannings                = lfr.where(landUse == i, 0.015, mannings)
                interceptionStorageMax  = lfr.where(landUse == i, 0.0025, interceptionStorageMax)
                throughfallFraction     = lfr.where(landUse == i, 0.95, throughfallFraction)
                LAI                     = lfr.where(landUse == i, 0.1, LAI)
                cropFactor              = lfr.where(landUse == i, 8, cropFactor)
                
            elif i in config.green:                                                  # Crops are given a multiplier of 1.2 as they also have pore structures \
                permeability = lfr.where(landUse == i, 1.1, permeability)            # but a different on can be assigned.
                mannings = lfr.where(landUse == i, 0.035, mannings)
                interceptionStorageMax  = lfr.where(landUse == i, 0.003, interceptionStorageMax)
                throughfallFraction     = lfr.where(landUse == i, 0.47, throughfallFraction)
                LAI                     = lfr.where(landUse == i, 1.5, LAI)
                cropFactor              = lfr.where(landUse == i, 1, cropFactor)
                
            elif i in config.water:                                                  # Any type of water is assigned 1, as this should have the saturated hydraulic conductivity \
                permeability = lfr.where(landUse == i, 1.0, permeability)            # as precipitation value. --> However, these places probably have inflow from groundwater.
                mannings = lfr.where(landUse == i, 0.030, mannings)
                interceptionStorageMax  = lfr.where(landUse == i, 0, interceptionStorageMax)
                throughfallFraction     = lfr.where(landUse == i, 1.0, throughfallFraction)
                LAI                     = lfr.where(landUse == i, 0, LAI)
                cropFactor              = lfr.where(landUse == i, 1, cropFactor)
            
            elif i in config.compacted:
                permeability = lfr.where(landUse == i, 0.5, permeability)
                mannings = lfr.where(landUse == i, 0.025, mannings)
                interceptionStorageMax  = lfr.where(landUse == i, 0.003, interceptionStorageMax)
                throughfallFraction     = lfr.where(landUse == i, 0.95, throughfallFraction)
                LAI                     = lfr.where(landUse == i, 0.15, LAI)
                cropFactor              = lfr.where(landUse == i, 8, cropFactor)
               
            elif i in config.other_road:
                permeability = lfr.where(landUse == i, 0.25, permeability) 
                mannings = lfr.where(landUse == i, 0.015, mannings)
                interceptionStorageMax  = lfr.where(landUse == i, 0.003, interceptionStorageMax)
                throughfallFraction     = lfr.where(landUse == i, 0.95, throughfallFraction)
                LAI                     = lfr.where(landUse == i, 0.15, LAI)
                cropFactor              = lfr.where(landUse == i, 1, cropFactor)

        return mannings, permeability, interceptionStorageMax, LAI, throughfallFraction, cropFactor
    
    def precipitation(date, session):
        return gen.lue_one() * 1 / (1000 * 3600) * int(config.includePrecipitation) # Convert to meter / second rate
    
    def pot_evaporation(date, session):
        pot_evaporation = gen.lue_one() * 3 / 10 / (1000*3600) # Convert from mm/d to m/h 
        return pot_evaporation * int(config.includeEvaporation)
    
    def infiltration(dem, groundWaterHeight, Ks, permeability, porosity):
        pot_infiltration = Ks * permeability
        pot_infiltration = lfr.where((dem - groundWaterHeight)*porosity < pot_infiltration, \
                                        (dem - groundWaterHeight)*porosity, pot_infiltration)
        return pot_infiltration * int(config.includeInfiltration)
    
    def EvaporationInfiltration(precipitation, surfaceWaterHeight, pot_evaporation, pot_infiltration, e_ratio, i_ratio):
         # Check if the waterheight is more than the evaporation and infiltration combined
        ie = pot_evaporation + pot_infiltration
        waterAvailable = surfaceWaterHeight + precipitation
        enoughWater = waterAvailable >= ie
        
        # Calculate the actual evaporation and infiltration
        evaporation  = lfr.where(enoughWater, pot_evaporation, waterAvailable*e_ratio)
        infiltration = lfr.where(enoughWater, pot_infiltration, waterAvailable*i_ratio)
        
        return evaporation, infiltration
    
    def percolation(dem, groundWaterHeight, Ks):
        part = (groundWaterHeight-0.5*dem)/dem
        part = lfr.where(part <= 0, gen.lue_zero(), part)
        percolation = lfr.where(dem < 0, gen.lue_zero(), \
                                Ks * part)                                                # If the dem is negative, there is no percolation
        percolation = lfr.where(percolation < 0, gen.lue_zero(), percolation) / 100                       # If the percolation is negative, there is no percolation
        return percolation * int(config.includePercolation)
    
    def ieRatio(pot_evaporation, pot_infiltration):
        try:
            i_ratio = pot_infiltration / (pot_evaporation + pot_infiltration)
            e_ratio = pot_evaporation / (pot_evaporation + pot_infiltration)
        except ZeroDivisionError:
            i_ratio, e_ratio = gen.lue_zero()

        return i_ratio, e_ratio
    
    def interception(interceptionStorage, interceptionStorageMax, precipitation, landUse):
        # Im this simplified version the landUse type will be used for the LAI
        # This should be improved in a future version
        k = 0.5         # Constant taken from (Brolsma et al., 2010)
        
        f = gen.lue_zero()  # If not stated then there are no leafs
        
        if config.includeInterception == False:
            interception = gen.lue_zero()
        else:
            for i in config.total:
                if i in config.green:
                    f = lfr.where(landUse == i, math.exp(-k * 3), 0)
            
            interception = (gen.lue_one() - f) * precipitation        
            interception = lfr.where(interceptionStorageMax < interceptionStorage + interception, \
                                    interceptionStorageMax - interceptionStorage, interception)
                        
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