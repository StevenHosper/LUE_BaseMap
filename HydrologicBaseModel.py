# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 2023

@author: steven.hosper
"""

import lue.framework as lfr
import lue.pcraster as lpr
import docopt
import numpy as np
import math as math
import os
import pandas as pd
import requests
import datetime
import sys
import time
import configuration as config
import getData as gD


# Timer to add some measure of functionality to the program
start_time = time.time()

usage = """\
Run the main model of the hydrologic base model.

Usage:
    {command}

Options:
    {command} : --hpx:thread = integer;
                The integer is the amount of cores used during the model run.
""".format(
    command=os.path.basename(sys.argv[0])
)

class mainModel(): 
    def __init__(self):
        """
        Summary:
            The initial calculation phase.
            Access to API session is granted.
            Command line arguments are translate to variables.
            Some useful arrays are created for calculations in the model.
        """

        self.s = requests.Session()
        self.s.headers = {'username': config.username,
                          'password': config.password,
                          'Content-Type': 'application/json',
                          }

        print(f'partition: {config.partitionShape}', f'array: {config.arrayShape}')
        
        # Create useful single value arrays
        # TO-DO: Incorporate this somewhere else
        # Zero, for empty cells or unincluded variables
        self.zero = lfr.create_array(config.arrayShape,
                                    config.partitionShape,
                                    dtype = np.dtype(np.float32),
                                    fill_value = 0,
                                    )
        # One, for additions
        self.ones = lfr.create_array(config.arrayShape,
                                    config.partitionShape,
                                    dtype=np.float32,
                                    fill_value=1,
                                    )
                
        # For the ldd direction being towards the cell itself (pit)
        self.sink = lfr.create_array(config.arrayShape,
                                    config.partitionShape,
                                    dtype = np.dtype(np.uint8),
                                    fill_value = 5,
                                    )
        
        print("__init__ done")
        
    
    def calculate_infiltration(self, Ks, land_c):
        """
        Summary:
            Calculate the rate of infiltration based on the soil type and the land use
        
        Input:
            Ks [m/day]: saturated hydraulic conductivity
            land_c [-]: land-use coefficient
        
        Returns:
            infil [m/day]: the amount of water that infiltrates the soil
        """
        # For now very basic formula, hydraulic conductivity times land-use coefficient/multiplier
        infil = lfr.multiply(Ks, land_c)
        return infil
    
    def kinematic_test(self, ldd, head, inflow):
        """
        Summary:
            Calculates the flux of material through each cell based on the ldd, inflow, velocity, channel_length and coefficients alpha & beta.
            
        Input:
            ldd: local drain direction map based on the Digital Elevation Map / current map height
            head: head height of the current situation
            inflow: currently only the precipitation
        
        Returns:
            kinematic [unit]: the kinematic wave flux for each cell
        """
        velocity = lfr.atan((head - lfr.downstream(ldd, head))/1)
            
        channel = lfr.create_array(config.arrayShape,
                                   config.partitionShape,
                                   dtype = np.float32,
                                   fill_value = 1,
                                   )
        
        kinematic = lfr.kinematic_wave(
                ldd,
                velocity,
                inflow,
                1.5,
                0.6,
                1.0,
                channel,
            )
        return kinematic
    
    
    def static(self, current_date: datetime.date, path: str, hydraulic_head: float, Ks, land_c, landUse):
        # Print the start date for logging purposes
        print(f'The startdate is: {current_date}')
        
        # Set the waterheight so that it matches the elevation head
        self.waterheight = lfr.where(self.dem < hydraulic_head, hydraulic_head - self.dem, 0)
        self.waterheight = lfr.where(landUse == 51, self.waterheight, 0 )
        
        # Create a groundwater table 
        # TO-DO : In the future we would preferable read the groundwater table from measuring stations \
        #         and apply these. However, that is not yet possible.
        self.groundwaterheight = self.dem - 0.1
        self.groundwaterheight = lfr.where(landUse == 51, self.dem, self.groundwaterheight)
        
        lfr.to_gdal(self.waterheight, config.path + f'/output/initial_water_height.tiff')

        # Access the data from the directory or the API dependend on the settings
        # Precipitation
        if config.includePrecipitation:
            if config.useAPI:
                precipitation     = gD.getData.get_api_data(current_date, 'precipitation', self.s)
            else:
                precipitation     = gD.getData.get_data(f'{config.path}/data/De Wupsel/', current_date, 'precipitation')
        else:
            precipitation = self.zero
        
        # Evaporation
        if config.includeEvaporation:
            if config.useAPI:
                pot_evaporation   = gD.getData.get_api_data(current_date, 'potential_evaporation', self.s)
            else:
                pot_evaporation   = gD.getData.get_data(f'{config.path}/data/De Wupsel/', current_date, 'potential_evaporation')
        else:
            pot_evaporation = self.zero
        
        
        # Calculate the infiltration rate
        if config.includeInfiltration:
            pot_infiltration = self.calculate_infiltration(Ks, land_c)
            pot_infiltration = lfr.where((self.dem - self.groundwaterheight) < pot_infiltration, \
                                         self.dem - self.groundwaterheight, pot_infiltration)
        else:
            pot_infiltration = self.zero
            
        # Calculate the percolation --> For now just a small percentage of Ks
        if config.includePercolation:
            percolation = lfr.where(self.dem < 0, self.zero, \
                                    Ks * 0.3 * ((self.groundwaterheight - 0.5 * self.dem) / self.dem))  # If the dem is negative, there is no percolation
            percolation = lfr.where(percolation < 0, self.zero, percolation)                            # If the percolation is negative, there is no percolation
        else:
            percolation = self.zero
            
        # Ratio of evaporation compared to infiltration
        if config.includeEvaporation and config.includeInfiltration:
            i_ratio = lfr.divide(pot_infiltration, lfr.add(pot_evaporation, pot_infiltration))
            e_ratio = lfr.divide(pot_evaporation, lfr.add(pot_evaporation, pot_infiltration))
        else:
            i_ratio = self.zero
            e_ratio = self.zero
        
        # Add precipitation to the watertable
        self.waterheight = self.waterheight + precipitation

        # Check if the waterheight is more than the evaporation and infiltration combined
        ie = pot_evaporation + pot_infiltration
        
        # Calculate the actual evaporation and infiltration
        evaporation  = lfr.where(self.waterheight >= ie, pot_evaporation, self.waterheight*e_ratio)
        infiltration = lfr.where(self.waterheight >= ie, pot_infiltration, self.waterheight*i_ratio)
        
        # Remove the evaporation and infiltration from the waterheight as it is lost to the atmosphere or groundwater.
        # *Note --> the rain now happens 'first', depending on when the rain falls during the day there might not be time to evaporate, but this is currently not taken into account.
        self.waterheight = self.waterheight - evaporation - infiltration
        self.groundwaterheight = self.groundwaterheight + infiltration - percolation
        
        # Make sure the waterheight is not below zero as this is not possible.
        # There can be no evaporation without water.
        self.waterheight = lfr.where(self.waterheight < self.zero, self.zero, self.waterheight)
        self.height = self.waterheight + self.dem
        
        # Create file with current situation of the water balance
        lfr.to_gdal(self.groundwaterheight, config.path + f'/output/groundwater_{current_date}.tiff')
        #lfr.to_gdal(self.height, config.path + f'/output/surfaceheight_{current_date}.tiff')
        lfr.to_gdal(self.waterheight, config.path + f'/output/waterheight_{current_date}.tiff')
        return 0
    
    

    def iterate(self, start_date: datetime.date, end_date: datetime.date, path: str, hydraulic_head: float, Ks, land_c, landUse):
        
        for i in range(int((end_date - start_date).days)):
            # Print the time to keep track while the program runs
            print(f'The date is: {start_date + datetime.timedelta(1+i)}')
            current_date = start_date + datetime.timedelta(1+i)
            
            # Recalculate if the water should still flow in the same direction
            ldd = lfr.d8_flow_direction(self.height)
            groundhead = self.groundwaterheight + self.waterheight
            groundldd = lfr.d8_flow_direction(groundhead)
            
            # Access the data from the directory or the API dependend on the settings
            # Precipitation
            if config.includePrecipitation:
                if config.useAPI:
                    precipitation     = gD.getData.get_api_data(current_date, 'precipitation', self.s)
                else:
                    precipitation     = gD.getData.get_data(f'{config.path}/data/De Wupsel/', current_date, 'precipitation')
            else:
                precipitation = self.zero
            
            # Evaporation
            if config.includeEvaporation:
                if config.useAPI:
                    pot_evaporation   = gD.getData.get_api_data(current_date, 'potential_evaporation', self.s)
                else:
                    pot_evaporation   = gD.getData.get_data(f'{config.path}/data/De Wupsel/', current_date, 'potential_evaporation')
            else:
                pot_evaporation = self.zero
            
            
            # Calculate the infiltration rate
            if config.includeInfiltration:
                pot_infiltration = self.calculate_infiltration(Ks, land_c)
                pot_infiltration = lfr.where((self.dem - self.groundwaterheight) < pot_infiltration, \
                                         self.dem - self.groundwaterheight, pot_infiltration)
            else:
                pot_infiltration = self.zero
            
            # Ratio of evaporation compared to infiltration
            if config.includeEvaporation and config.includeInfiltration:
                i_ratio = lfr.divide(pot_infiltration, lfr.add(pot_evaporation, pot_infiltration))
                e_ratio = lfr.divide(pot_evaporation, lfr.add(pot_evaporation, pot_infiltration))
            else:
                i_ratio = self.zero
                e_ratio = self.zero

            # Calculate the percolation --> For now just a small percentage of Ks
            if config.includePercolation:
                percolation = lfr.where(self.dem <= 0, self.zero, \
                                        Ks * 0.3 * ((self.groundwaterheight - 0.5 * self.dem) / self.dem))  # If the dem is negative, there is no percolation
                percolation = lfr.where(percolation < 0, self.zero, percolation)                            # If the percolation is negative, there is no percolation
            else:
                percolation = self.zero
            
            # Calculate the groundwater flow
            if config.includeGroundFlow:
                groundflow = Ks * (groundhead - lfr.downstream(groundldd, groundhead))                               # Q = k * i * A
            else:
                groundflow = self.zero
            
            # Helps with regulating the inflow of water
            upstream_cells = lfr.upstream(ldd, self.ones)
            
            # Check the difference in dem (as this determines the total height that should be filled to create an equal surface)
            height_difference = self.height - lfr.downstream(ldd, self.height)
            potential_flux = lfr.where(height_difference > self.waterheight, 0.5*self.waterheight, 0.5*height_difference)
            flux = lfr.where(ldd != 5, potential_flux, 0)
            
            # ROUTING
            self.waterheight = self.waterheight + lfr.upstream(ldd, flux) - flux
            self.groundwaterheight = self.groundwaterheight + lfr.upstream(groundldd, groundflow) - groundflow
            
            # The excess groundwater seeps upwards out of the soil
            groundwatersurplus = self.groundwaterheight - self.dem
            self.groundwaterheight = self.groundwaterheight - groundwatersurplus
            self.waterheight = self.waterheight + groundwatersurplus
            
            # Add precipitation to the watertable
            self.waterheight = self.waterheight + precipitation
            
            # Check if the waterheight is more than the evaporation and infiltration combined
            ie = pot_evaporation + pot_infiltration
            
            # Calculate the actual evaporation and infiltration
            evaporation  = lfr.where(self.waterheight >= ie, pot_evaporation, self.waterheight*e_ratio)
            infiltration = lfr.where(self.waterheight >= ie, pot_infiltration, self.waterheight*i_ratio)
            
            kinematic = self.kinematic_test(ldd, self.height, precipitation) 

            # Remove the evaporation and infiltration from the waterheight as it is lost to the atmosphere or groundwater.
            # *Note --> the rain now happens 'first', depending on when the rain falls during the day there might not be time to evaporate, but this is currently not taken into account.
            self.waterheight = self.waterheight - evaporation - infiltration
            self.groundwaterheight = self.groundwaterheight + infiltration - percolation
            
            # Waterheight can never be lower than zero.
            self.waterheight = lfr.where(self.waterheight < self.zero, self.zero, self.waterheight)
            
            # Adjust the concurrent height
            self.height = self.dem + self.waterheight
            
            # Save the variables to tiff files to a tiff file
            lfr.to_gdal(self.groundwaterheight, config.path + f'/output/groundwater_{current_date}.tiff')
            lfr.to_gdal(self.waterheight, config.path + f'/output/waterheight_{current_date}.tiff')
            #lfr.to_gdal(self.height, config.path + f'/output/surfaceheight_{current_date}.tiff')
        return 0
     
     
    @lfr.runtime_scope 
    def simulate(self):
        # Load initial variables
        self.dem  = lfr.from_gdal(config.path + f'/data/{config.scenario}/dem.tiff', config.partitionShape)
        landUse  = lfr.from_gdal(config.path + f'/data/{config.scenario}/landgebruik.tiff', config.partitionShape)
        soilType = lfr.from_gdal(config.path + f'/data/{config.scenario}/bodem.tiff', config.partitionShape)
        
        # Create initial ldd
        self.ldd = lfr.d8_flow_direction(self.dem)
        
        # Create the hydraulic conductivity variable
        Ks       = self.ones
        
        # Assign K for de Wupsel
        for i in range(20):
            if i == 2:
                Ks = lfr.where(soilType == i, 1.0*10**-1, Ks)           # Hydraulic conductivity in m/day
            if i == 9:
                Ks = lfr.where(soilType == i, 3.0*10**0, Ks)           # Hydraulic conductivity in m/day
            if i == 10:
                Ks = lfr.where(soilType == i, 4.0*10**0, Ks)           # Hydraulic conductivity in m/day
            if i == 11:
                Ks = lfr.where(soilType == i, 4.0*10**-1, Ks)           # Hydraulic conductivity in m/day
            if i == 12:
                Ks = lfr.where(soilType == i, 2.5*10**0, Ks)           # Hydraulic conductivity in m/day
            if i == 13:
                Ks = lfr.where(soilType == i, 2.0*10**-1, Ks)           # Hydraulic conductivity in m/day
            if i == 15:
                Ks = lfr.where(soilType == i, 5.0*10**-1, Ks)           # Hydraulic conductivity in m/day
            if i == 16:
                Ks = lfr.where(soilType == i, 1.0*10**-1, Ks)           # Hydraulic conductivity in m/day
            if i == 19:
                Ks = lfr.where(soilType == i, 3.0*10**-2, Ks)           # Hydraulic conductivity in m/day
            else:
                Ks = lfr.where(soilType == i, 5.0*10**-2, Ks)           # Hydraulic conductivity in m/day

        # Create the land-use coefficient variable
        land_c   = lfr.create_array(config.arrayShape,
                                    config.partitionShape,
                                    dtype=np.float32,
                                    fill_value=0.5,
                                    )
        
        # Create ID variables
        concrete       = [2, 4, 6, 8, 10, 14, 15, 16, 25, 28, 35, 166, 253]
        green          = [40, 41, 43, 44, 112, 157]
        water          = [51, 254]
        compacted      = [18]
        other_road     = [29]
        
        # Use the ID values given to the QGIS raster to determine which land-use types are assigned which values.
        for i in range(255):
            if i in concrete:
                land_c = lfr.where(landUse == i, 0.001, land_c)
                
            elif i in green or 44 > i > 157:                              # Crops are given a multiplier of 1.2 as they also have pore structures \
                land_c = lfr.where(landUse == i, 1.1, land_c)            # but a different on can be assigned.
                
            elif i in water:                                              # Any type of water is assigned 1, as this should have the saturated hydraulic conductivity \
                land_c = lfr.where(landUse == i, 1.0, land_c)            # as precipitation value. --> However, these places probably have inflow from groundwater.
            
            elif i in compacted:
                land_c = lfr.where(landUse == i, 0.7, land_c)
                
            elif i in other_road:
                land_c = lfr.where(landUse == i, 0.3, land_c) 
                
            else:
                land_c = lfr.where(landUse == i, 0.8, land_c)
        
        # Save the land-use coefficient so it can be checked
        lfr.to_gdal(land_c, config.path + f'/output/landUse_coefficients.tiff')
        lfr.to_gdal(Ks, config.path + f'/output/Ks.tiff')
        
        # Initialize the static
        self.static(config.startDate, config.path, config.groundWaterTable, Ks, land_c, landUse)
        print("static completed")
        
        self.iterate(config.startDate, config.endDate, config.path, config.groundWaterTable, Ks, land_c, landUse)
        print("iteration completed")
        return 0


# Initialize HPX runtime and run model, on the root locality -------------------
# General configuration options, which are valid on all
# platforms. Platform-specific options can be passed on the command line.
cfg = [
    # Make sure hpx_main is always executed
    "hpx.run_hpx_main!=1",
    # Allow for unknown command line options
    "hpx.commandline.allow_unknown!=1",
    # Disable HPX' short options
    "hpx.commandline.aliasing!=0",
    # Don't print diagnostics during forced terminate
    "hpx.diagnostics_on_terminate!=0",
    # Make AGAS clean up resources faster than by default
    "hpx.agas.max_pending_refcnt_requests!=50",
]

lfr.start_hpx_runtime(cfg)

# The root locality will distribute the work over all other
# localities. Never perform Python code on the other localities than the
# root locality unless you know what you are doing.
if lfr.on_root_locality():
    main = mainModel()
    main.simulate()
    
print("--- %s seconds ---" % (time.time() - start_time))