# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 2023

@author: steven.hosper
"""
# The HydrologicBaseModel
import lue.framework as lfr
import numpy as np
import math as math
import os
import requests
import datetime
import sys
import time
import configuration as config
import getData as data
import MakeGIF
import reporting


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
    
    
    def static(self, current_date: datetime.date, path: str, initialWaterTable: float, Ks, landC, landUse):
        # Print the start date for logging purposes
        print(f'The startdate is: {current_date}')
        
        # Set the waterheight so that it matches the elevation head
        self.waterheight = lfr.where(self.dem < initialWaterTable, initialWaterTable - self.dem, 0)
        self.waterheight = lfr.where(landUse == 51, self.waterheight, 0 )
        
        # Create a groundwater table 
        # TO-DO : In the future we would preferable read the groundwater table from measuring stations \
        #         and apply these. However, that is not yet possible.
        self.groundWaterHeight = self.dem - 0.1
        self.groundWaterHeight = lfr.where(landUse == 51, self.dem, self.groundWaterHeight)
        
        lfr.to_gdal(self.waterheight, config.path + f'/output/{config.scenario}/initial_water_height.tiff')

        # Access the data from the directory or the API dependend on the settings
        precipitation     = data.get.precipitation(current_date, self.s, self.zero)
        pot_evaporation   = data.get.pot_evaporation(current_date, self.s, self.zero)   
        pot_infiltration  = data.get.infiltration(self.dem, self.groundWaterHeight, Ks, landC, self.zero)
        percolation       = data.get.percolation(self.dem, self.groundWaterHeight, Ks, self.zero)
        i_ratio, e_ratio  = data.get.ieRatio(pot_evaporation, pot_infiltration, self.ones, self.zero)

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
        self.groundWaterHeight = self.groundWaterHeight + infiltration - percolation
        
        # Make sure the waterheight is not below zero as this is not possible.
        # There can be no evaporation without water.
        self.waterheight = lfr.where(self.waterheight < self.zero, self.zero, self.waterheight)
        self.height = self.waterheight + self.dem
        
        # Create file with current situation of the water balance
        variables = [self.waterheight, self.groundWaterHeight]
        reporting.report.static(current_date, variables, config.output_path)
        return 0
    
    

    def iterate(self, start_date: datetime.date, end_date: datetime.date, path: str, hydraulic_head: float, Ks, landC, landUse):
        
        for i in range(int((end_date - start_date).days)):
            # Print the time to keep track while the program runs
            print(f'The date is: {start_date + datetime.timedelta(1+i)}')
            current_date = start_date + datetime.timedelta(1+i)
            
            # Recalculate if the water should still flow in the same direction
            ldd = lfr.d8_flow_direction(self.height)
            groundhead = self.groundWaterHeight + self.waterheight
            groundldd = lfr.d8_flow_direction(groundhead)
            
            # Access the data from the directory or the API dependend on the settings
            precipitation     = data.get.precipitation(current_date, self.s, self.zero)
            pot_evaporation   = data.get.pot_evaporation(current_date, self.s, self.zero)
            pot_infiltration  = data.get.infiltration(self.dem, self.groundWaterHeight, Ks, landC, self.zero)
            percolation       = data.get.percolation(self.dem, self.groundWaterHeight, Ks, self.zero)
            
            
            # Ratio of evaporation compared to infiltration
            if config.includeEvaporation and config.includeInfiltration:
                i_ratio = lfr.divide(pot_infiltration, lfr.add(pot_evaporation, pot_infiltration))
                e_ratio = lfr.divide(pot_evaporation, lfr.add(pot_evaporation, pot_infiltration))
            else:
                i_ratio = self.zero
                e_ratio = self.zero
            
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
            self.groundWaterHeight = self.groundWaterHeight + lfr.upstream(groundldd, groundflow) - groundflow
            
            # The excess groundwater seeps upwards out of the soil
            groundWaterSurplus = self.groundWaterHeight - self.dem
            self.groundWaterHeight = self.groundWaterHeight - groundWaterSurplus
            self.waterheight = self.waterheight + groundWaterSurplus
            
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
            self.groundWaterHeight = self.groundWaterHeight + infiltration - percolation
            
            # Waterheight can never be lower than zero.
            self.waterheight = lfr.where(self.waterheight < self.zero, self.zero, self.waterheight)
            
            # Adjust the concurrent height
            self.height = self.dem + self.waterheight
            
            # Save the variables to tiff files to a tiff file
            lfr.to_gdal(self.groundWaterHeight, config.path + f'/output/{config.scenario}/groundwater_{current_date}.tiff')
            lfr.to_gdal(self.waterheight, config.path + f'/output/{config.scenario}/waterheight_{current_date}.tiff')
            #lfr.to_gdal(self.height, config.path + f'/output/{config.scenario}/surfaceheight_{current_date}.tiff')
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
        Ks       = data.get.Ks(soilType, self.ones)
        lfr.to_gdal(Ks, config.path + f'/output/{config.scenario}/Ks.tiff')
        
        # Create the land-use coefficient variable
        landC   = data.get.landC(landUse, self.ones)
        lfr.to_gdal(landC, config.path + f'/output/{config.scenario}/landUse_coefficients.tiff')
        
        # Initialize the static
        self.static(config.startDate, config.path, config.initialWaterTable, Ks, landC, landUse)
        print("static completed")
        
        self.iterate(config.startDate, config.endDate, config.path, config.initialWaterTable, Ks, landC, landUse)
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
    MakeGIF.main()
    
print("--- %s seconds ---" % (time.time() - start_time))