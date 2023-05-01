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
        options = gdal.TranslateOptions(outputType= gdal.GDT_Float64)
        gdal.FileFromMemBuffer(input_mem_file, pull.content)
        gdal.Translate(output_mem_file, input_mem_file, options = options)
        
        data = lfr.from_gdal(output_mem_file, config.partitionShape)
        return data
    
    def soil_csv(dataFile, soilType, configuration):
        # Assign K for de Wupsel
        Ks = gen.lue_one(configuration) * 0.05
        porosity = gen.lue_one(configuration) * 0.35
        wiltingPoint = gen.lue_one(configuration) * 0.15
        scTable = pd.read_csv(dataFile)
        ID = scTable["ID"]
        KsValue = scTable["Ks"]
        for count, ID in enumerate(ID):
            Ks = lfr.where(soilType == ID, KsValue[count], Ks)
        return (Ks / 86400), porosity, wiltingPoint
    
    def landCharacteristics_csv(dataFile, landUse, configuration):
        # Use the ID values given to the QGIS raster to determine which land-use types are assigned which values.
        # Standard values
        dummy                  = gen.lue_one(configuration)
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
            mannings                = lfr.where(landUse == ID, manningsFriction[count], mannings)
            permeability            = lfr.where(landUse == ID, permeabilityValue[count], permeability)
            interceptionStorageMax  = lfr.where(landUse == ID, interceptionStorageMaxValue[count], interceptionStorageMax)
            # LAI                     = lfr.where(soilType == ID, LAIValue[count], LAI)
            throughfallFraction     = lfr.where(landUse == ID, throughfallFractionValue[count], throughfallFraction)
            # cropFactor              = lfr.where(soilType == ID, cropFactorValue[count], cropFactor)
                
        return mannings, permeability, interceptionStorageMax, throughfallFraction
    
    def csvData(date, cellArea, dataFile, configuration):
        """
        Args:
            date (_type_): Gives the concurrent date of the model.
            cellArea (_type_): The area of the cell of the model.
            dataFile (_type_): The file out of which data is to be extracted, should be in mm/h

        Returns:
            _type_: A lue array with the data in a flux of m/s
        """
        roundedDate = get.roundDownDateTime(date)
        # Convert date to a string
        dateTime = roundedDate.strftime("%d/%m/%Y %H:%M")
        
        # Read the csv and set the index
        try:
            data = pd.read_csv(dataFile, sep=",", names=['datetime', 'dataValue'])
            data.set_index('datetime', inplace=True)
            dataValue = data.loc[f'{dateTime}']['dataValue']
        except:
            dataValue = 0
        
        return gen.lue_one(configuration) * (dataValue / (1000 * 3600)) * cellArea  # Convert to m3/s rate
    
    def interception(cellArea, interceptionStorage, interceptionStorageMax, precipitation, ref_evaporation, throughfallFraction, configuration):
        interception = (gen.lue_one(configuration) - throughfallFraction) * precipitation * int(configuration.generalSettings['includePrecipitation'])
        
        enoughWaterInt = (interception + interceptionStorage/int(configuration.modelSettings['iterationsBeforeReport'])) > ref_evaporation  

        interceptionStorage = lfr.where(enoughWaterInt, interceptionStorage + (interception - ref_evaporation) * int(configuration.modelSettings['iterationsBeforeReport']), 0)
        
        interceptionStorageSurplus = lfr.where(interceptionStorage > interceptionStorageMax, interceptionStorage - interceptionStorageMax, 0)
        interceptionStorage        = interceptionStorage - interceptionStorageSurplus
        precipitation = precipitation - (interception - interceptionStorageSurplus/int(configuration.modelSettings['iterationsBeforeReport']))
        
        evapotranspirationSurface = lfr.where(enoughWaterInt, 0, ref_evaporation - (interception + interceptionStorage/int(configuration.modelSettings['iterationsBeforeReport'])))
        return interceptionStorage, precipitation, evapotranspirationSurface
    
    def infiltration(Sgw, MaxSgw, cellArea, Ks, permeability, porosity, precipitation, evapotranspirationSurface):
        potInfiltration = Ks * permeability * int(config.includeInfiltration) # meters that can infiltrate the soil
        potInfiltration = lfr.where(((MaxSgw - Sgw)/cellArea)*porosity < potInfiltration, \
                                        ((MaxSgw - Sgw)/cellArea)*porosity, potInfiltration) * cellArea     # The amount that can infiltrate because of capacity times the area
        enoughWaterInf   = precipitation - evapotranspirationSurface > potInfiltration                      # If there is more water on the surface available than can infiltrate
        infiltrationSurface = lfr.where(enoughWaterInf, potInfiltration, precipitation - evapotranspirationSurface) # Either the potInfiltration will fully infiltrate
        potInfiltrationChannel = lfr.where(enoughWaterInf, 0, potInfiltration - infiltrationSurface) / cellArea
        return infiltrationSurface, potInfiltrationChannel                                                          # or the available water at the surface will.
    
    def evapotranspiration(precipitation, evapotranspirationSurface):
        enoughWaterE = precipitation > evapotranspirationSurface
        evapotranspirationSoil = lfr.where(enoughWaterE, 0, evapotranspirationSurface - precipitation)
        evapotranspirationSurface = lfr.where(enoughWaterE, evapotranspirationSurface, precipitation)
        
        return evapotranspirationSurface, evapotranspirationSoil
    
    def roundDownDateTime(dt):
        deltaMin = dt.minute % 5
        date = datetime.datetime(dt.year, dt.month, dt.day,
                                 dt.hour, dt.minute - deltaMin)
        return date