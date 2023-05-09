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
        
    def soil_csv(self, data_file, soil_type):
        """Reads the soil properties dependent on the soil IDs of the lue array
        
        Args:
            data_file (path):
            soil_type (lpa*):
            
        Returns:
            Ks (lpa*):
            porosity (lpa*):
            wilting_point (lpa*)
            
        lpa*: lue partitioned array
        """
        # Assign standard values for de Wupsel
        Ks = self.std_arr_lue.one() * 0.05
        porosity = self.std_arr_lue.one() * 0.35
        wilting_point = self.std_arr_lue.one() * 0.15
        
        # Read pandas data table
        data_table = pd.read_csv(data_file)
        
        # Split into ID and Ks value
        ID = data_table["ID"]
        Ks_value = data_table["Ks"]
        
        # When the ID of the table meets an ID within the partitioned array, assign value
        for count, ID in enumerate(ID):
            Ks = lfr.where(soil_type == ID, Ks_value[count], Ks) / 86400      # To m/s 
        return Ks, porosity, wilting_point
    
    def land_characteristics_csv(self, data_file, land_use):
        """Reads the land characteristics from a csv file
        
        Gives standard values to the entire array using the set array extent.
        Continues to read values from data file to the corresponding map IDs.
        
        Args:
            data_file (path):    the path to the csv data file containing the land characteristics information
            land_use (lpa*):     array containing the id that matches every cell to its \
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
        crop_factor              = dummy * 1.0
        
        # Open data table and load columns into variables
        data = pd.read_csv(data_file)
        ID = data["Code"]
        mannings_friction               = data["Friction"]
        permeability_value              = data["Permeability"]
        interception_storage_max_value  = data["Interception"]
        LAI_value                        = data["LAI"]
        throughfall_fraction_value       = data["f"]
        crop_factor_value                 = data["Crop_type"]
        
        # When an ID matches the array ID, assign value
        for count, ID in enumerate(ID):
            mannings                    = lfr.where(land_use == ID, mannings_friction[count], mannings)
            permeability                = lfr.where(land_use == ID, permeability_value[count], permeability)
            interception_storage_max    = lfr.where(land_use == ID, interception_storage_max_value[count], interception_storage_max)
            # LAI                       = lfr.where(land_use == ID, LAI_value[count], LAI)
            throughfall_fraction        = lfr.where(land_use == ID, throughfall_fraction_value[count], throughfall_fraction)
            # crop_factor                = lfr.where(land_use == ID, crop_factor_value[count], crop_factor)
                
        return mannings, permeability, interception_storage_max, throughfall_fraction
    
    def rounddown_datetime(self, dt):
        """Rounds down the time to the previous 5 minutes
        
        Args:
            dt (datetime date):     datetime date that should be rounded to the prior 5 minutes
            
        Returns:
            date (datetime date):   datetime date rounded down
        
        """
        delta_min = dt.minute % 5
        date = datetime.datetime(dt.year, dt.month, dt.day,
                                 dt.hour, dt.minute - delta_min)
        return date
    
    def csv_timeseries_to_flux(self, data_file, refactor, date):
        """
        Args:
            data_file (path):        The file out of which data is to be extracted
            refactor (float):       To refactor the value within the timeseries to a flux in m3/s to the model.
                                    Example: from mm/h to m3/s for cell area of 25 -> \
                                        ([cell area: 25]  /  [mm to m: 1000])  /  [h to s: 3600] 
            date (datetime date):   Gives the concurrent date of the model.
            
            

        Returns:
            data_value (lue partitioned array): A lue array with the data in a flux of m/s
        """
        rounded_date = self.rounddown_datetime(date)
        # Convert date to a string
        date_time = rounded_date.strftime("%d/%m/%Y %H:%M")
        
        # Read the csv and set the index
        try:
            data = pd.read_csv(data_file, sep=",", names=['date_time', 'data_value'])
            data.set_index('date_time', inplace=True)
            data_value = data.loc[f'{date_time}']['data_value']
            data_value = data_value * refactor * self.std_arr_lue.one() # Convert to m/s rate from mm/h
        except:
            data_value = self.std_arr_lue.zero()
        
        return data_value  