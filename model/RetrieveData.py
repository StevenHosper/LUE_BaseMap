# -*- coding: utf-8 -*-
"""
Created on 8 May 2023

@author: steven.hosper
"""

import lue.framework as lfr
import numpy as np
import datetime
import pandas as pd
from StandardArraysLUE import StandardArraysLUE

class RetrieveData():
    def __init__(self, configuration):
        """
        Initialize the class. 
        1) Initialize the standard arrays
        2) Set extent, input and output dir.
        """
        self.std_arr_lue        = StandardArraysLUE(configuration)
        
        self.array_extent       = int(configuration.modelSettings['arrayExtent'])
        self.partition_extent   = int(configuration.modelSettings['partitionExtent'])
        self.output_dir         = configuration.generalSettings['outputDir']
        self.input_dir          = configuration.generalSettings['inputDir']
        
    def soil_csv(self, dataFile, soilType):
        """Reads the soil properties dependent on the soil IDs of the lue array
        
        Args:
            dataFile (path):
            soilType (lpa*):
            
        Returns:
            Ks (lpa*):
            porosity (lpa*):
            wiltingPoint (lpa*)
            
        lpa*: lue partitioned array
        """
        # Assign standard values for de Wupsel
        Ks = self.std_arr_lue.one() * 0.05
        porosity = self.std_arr_lue.one() * 0.35
        wiltingPoint = self.std_arr_lue.one() * 0.15
        
        # Read pandas data table
        scTable = pd.read_csv(dataFile)
        
        # Split into ID and Ks value
        ID = scTable["ID"]
        KsValue = scTable["Ks"]
        
        # When the ID of the table meets an ID within the partitioned array, assign value
        for count, ID in enumerate(ID):
            Ks = lfr.where(soilType == ID, KsValue[count], Ks) / 86400      # To m/s 
        return Ks, porosity, wiltingPoint
    
    def land_characteristics_csv(self, dataFile, landUse):
        """Reads the land characteristics from a csv file
        
        Gives standard values to the entire array using the set array extent.
        Continues to read values from data file to the corresponding map IDs.
        
        Args:
            dataFile (path):    the path to the csv data file containing the land characteristics information
            landUse (lpa*):     array containing the id that matches every cell to its \
                                corresponding characteristic values.
                
        Returns:
            mannings (lpa*):                    Mannings Coefficient in a lue array
            permeability (lpa*):                Permeability in a lue array
            interception_storage_max (lpa*):    Interception storage max in a lue array
            throughfall_fraction (lpa*):        Throughfall fraction in a lue array
        
        lpa*: lue partitioned array
        """
        # Use the ID values given to the QGIS raster to determine which land-use types are assigned which values.
        # Standard values
        dummy                  = self.std_arr_lue.one()
        permeability           = dummy * 0.8
        mannings               = dummy * 0.045
        interception_storage_max = dummy * 0.001
        throughfall_fraction    = dummy * 0.90
        LAI                     = dummy * 0.5
        cropFactor              = dummy * 1.0
        
        # Open data table and load columns into variables
        data = pd.read_csv(dataFile)
        ID = data["Code"]
        mannings_friction               = data["Friction"]
        permeability_value              = data["Permeability"]
        interception_storage_max_value  = data["Interception"]
        LAIValue                        = data["LAI"]
        throughfall_fractionValue       = data["f"]
        cropFactorValue                 = data["Crop_type"]
        
        # When an ID matches the array ID, assign value
        for count, ID in enumerate(ID):
            mannings                    = lfr.where(landUse == ID, mannings_friction[count], mannings)
            permeability                = lfr.where(landUse == ID, permeability_value[count], permeability)
            interception_storage_max    = lfr.where(landUse == ID, interception_storage_max_value[count], interception_storage_max)
            # LAI                       = lfr.where(soilType == ID, LAIValue[count], LAI)
            throughfall_fraction        = lfr.where(landUse == ID, throughfall_fractionValue[count], throughfall_fraction)
            # cropFactor                = lfr.where(soilType == ID, cropFactorValue[count], cropFactor)
                
        return mannings, permeability, interception_storage_max, throughfall_fraction
    
    def roundDownDateTime(self, dt):
        """Rounds down the time to the previous 5 minutes
        
        Args:
            dt (datetime date):     datetime date that should be rounded to the prior 5 minutes
            
        Returns:
            date (datetime date):   datetime date rounded down
        
        """
        deltaMin = dt.minute % 5
        date = datetime.datetime(dt.year, dt.month, dt.day,
                                 dt.hour, dt.minute - deltaMin)
        return date
    
    def csvTimeseries2Flux(self, dataFile, refactor, date):
        """
        Args:
            dataFile (path):        The file out of which data is to be extracted
            refactor (float):       To refactor the value within the timeseries to a flux in m3/s to the model.
                                    Example: from mm/h to m3/s for cell area of 25 -> \
                                        ([cell area: 25]  /  [mm to m: 1000])  /  [h to s: 3600] 
            date (datetime date):   Gives the concurrent date of the model.
            
            

        Returns:
            dataValue (lue partitioned array): A lue array with the data in a flux of m/s
        """
        roundedDate = self.roundDownDateTime(date)
        # Convert date to a string
        dateTime = roundedDate.strftime("%d/%m/%Y %H:%M")
        
        # Read the csv and set the index
        try:
            data = pd.read_csv(dataFile, sep=",", names=['datetime', 'dataValue'])
            data.set_index('datetime', inplace=True)
            dataValue = data.loc[f'{dateTime}']['dataValue']
            dataValue = dataValue * refactor * self.std_arr_lue.one() # Convert to m/s rate from mm/h
        except:
            dataValue = self.std_arr_lue.zero()
        
        return dataValue  