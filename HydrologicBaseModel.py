# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 2023

@author: steven.hosper
"""
# The HydrologicBaseModel
import lue.framework as lfr
import math as math
import os
import datetime
import sys
import time
import pcraster as pcr
from osgeo import gdal

# Own functions
# TO-DO: Reformat by the use of: from ... import ... as ...
import configuration as config
import dataAccess as dA
import MakeGIF
import reporting
import update
import dataGen as dG


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
        pass
    
    def static(self, current_date: datetime.date, initialWaterTable: float, Ks, landC, landUse):
        # Print the start date for logging purposes
        print(f'The startdate is: {current_date}')
        
        # Set the waterheight so that it matches the elevation head
        # TO-DO: Initialize  waterheight for rivers and water bodies
        iniSurfaceWaterHeight = lfr.where(self.dem < initialWaterTable, initialWaterTable - self.dem, 0)
        iniSurfaceWaterHeight = lfr.where(landUse == 51, iniSurfaceWaterHeight, 0 )
        
        # Create a groundwater table 
        # TO-DO : In the future we would preferable read the groundwater table from measuring stations \
        #         and apply these. However, that is not yet possible.
        iniGroundWaterHeight = self.dem - 0.1
        iniGroundWaterHeight = lfr.where(landUse == 51, self.dem, iniGroundWaterHeight)
        iniGroundWaterHeight = lfr.where(dG.generate.boundaryCell(), self.dem, iniGroundWaterHeight)
        
        lfr.to_gdal(iniSurfaceWaterHeight, config.path + f'/output/{config.scenario}/initial_water_height.tiff')

        # Access the data from the directory or the API dependend on the settings
        precipitation     = dA.get.precipitation(current_date, dA.get.apiSession())
        pot_evaporation   = dA.get.pot_evaporation(current_date, dA.get.apiSession())   
        pot_infiltration  = dA.get.infiltration(self.dem, iniGroundWaterHeight, Ks, landC)
        percolation       = dA.get.percolation(self.dem, iniGroundWaterHeight, Ks)
        i_ratio, e_ratio  = dA.get.ieRatio(pot_evaporation, pot_infiltration)
        
        # Add precipitation to the watertable
        iniSurfaceWaterHeight = iniSurfaceWaterHeight + precipitation

        # Check if the waterheight is more than the evaporation and infiltration combined
        ie = pot_evaporation + pot_infiltration
        
        # Calculate the actual evaporation and infiltration
        evaporation  = lfr.where(iniSurfaceWaterHeight >= ie, pot_evaporation, iniSurfaceWaterHeight*e_ratio)
        infiltration = lfr.where(iniSurfaceWaterHeight >= ie, pot_infiltration, iniSurfaceWaterHeight*i_ratio)
        
        
        # Remove the evaporation and infiltration from the waterheight as it is lost to the atmosphere or groundwater.
        # *Note --> the rain now happens 'first', depending on when the rain falls during the day there might not be time to evaporate, but this is currently not taken into account.
        iniSurfaceWaterHeight = iniSurfaceWaterHeight - evaporation - infiltration
        iniGroundWaterHeight = iniGroundWaterHeight + infiltration - percolation
        
        # Make sure the waterheight is not below zero as this is not possible.
        # There can be no evaporation without water.
        iniSurfaceWaterHeight = lfr.where(iniSurfaceWaterHeight < dG.generate.lue_zero(), dG.generate.lue_zero(), iniSurfaceWaterHeight)
        self.height = iniSurfaceWaterHeight + self.dem
        
        # Create file with current situation of the water balance
        variables = {"surfacewater": iniSurfaceWaterHeight, "groundwater": iniGroundWaterHeight}
        reporting.report.static(current_date, variables, config.output_path)
        return iniSurfaceWaterHeight, iniGroundWaterHeight
    
    

    def iterate(self, start_date: datetime.date, end_date: datetime.date, \
                iniSurfaceWaterHeight, iniGroundWaterHeight, Ks, landC):        # At this point we are splitting a day within 300 steps, aiming to go to seconds base.
        seconds = 3600
        surfaceWaterHeight = iniSurfaceWaterHeight
        groundWaterHeight = iniGroundWaterHeight
        for i in range(int((end_date - start_date).days * seconds)):
            current_date = start_date + datetime.timedelta(seconds=1+i)
            print(f'Second: {i}')
            # Recalculate if the water should still flow in the same direction
            ldd = lfr.d8_flow_direction(self.height)
            
            # Access the data from the directory or the API dependend on the settings
            precipitation     = dA.get.precipitation(current_date, dA.get.apiSession())
            pot_evaporation   = dA.get.pot_evaporation(current_date, dA.get.apiSession())
            pot_infiltration  = dA.get.infiltration(self.dem, groundWaterHeight, Ks, landC)
            percolation       = dA.get.percolation(self.dem, groundWaterHeight, Ks)
            i_ratio, e_ratio  = dA.get.ieRatio(pot_evaporation, pot_infiltration)

            
            # Add precipitation to the watertable
            surfaceWaterHeight = surfaceWaterHeight + precipitation
            
            # Check if the waterheight is more than the evaporation and infiltration combined
            ie = pot_evaporation + pot_infiltration
            
            # Calculate the actual evaporation and infiltration
            evaporation  = lfr.where(surfaceWaterHeight >= ie, pot_evaporation, surfaceWaterHeight*e_ratio)
            infiltration = lfr.where(surfaceWaterHeight >= ie, pot_infiltration, surfaceWaterHeight*i_ratio)
            
            # Check the difference in dem (as this determines the total height that should be filled to create an equal surface)
            height_difference = self.height - lfr.downstream(ldd, self.height)
            potential_flux = lfr.where(height_difference > surfaceWaterHeight, 0.5*surfaceWaterHeight, 0.5*height_difference)
            flux = lfr.where(ldd != 5, potential_flux, 0)
            
            # ROUTING
            surfaceWaterHeight = lfr.where(dG.generate.boundaryCell(), surfaceWaterHeight + lfr.upstream(ldd, flux) - flux, surfaceWaterHeight)
            
            groundWaterHeight = lfr.where(dG.generate.boundaryCell(), update.update.groundWaterHeight(
                self.dem, Ks, surfaceWaterHeight, groundWaterHeight, infiltration, \
                percolation, dG.generate.lue_zero()
                ), groundWaterHeight)

            # Remove the evaporation and infiltration from the waterheight as it is lost to the 
            # atmosphere or groundwater.
            surfaceWaterHeight = surfaceWaterHeight - evaporation - infiltration
            groundWaterHeight = groundWaterHeight + infiltration - percolation
            
            
            # Waterheight can never be lower than zero.
            surfaceWaterHeight = lfr.where(surfaceWaterHeight < dG.generate.lue_zero(), dG.generate.lue_zero(), surfaceWaterHeight)
            
            # Adjust the concurrent height
            self.height = self.dem + surfaceWaterHeight
            if config.v2: second = i;
            else: second = 1;
            
            # Create file with current situation of the water balance
            variables = {"surfacewater": surfaceWaterHeight, "groundwater": groundWaterHeight}
            fluxes    = {}
            reporting.report.dynamic(current_date, second, variables, fluxes, config.output_path)
            lfr.maximum(flux).get()
        return 0
     
     
    @lfr.runtime_scope 
    def simulate(self):
        # Load initial variables
        self.dem = lfr.from_gdal(config.path + f'/data/{config.scenario}/dem.tiff', config.partitionShape)
        landUse  = lfr.from_gdal(config.path + f'/data/{config.scenario}/landgebruik.tiff', config.partitionShape)
        soilType = lfr.from_gdal(config.path + f'/data/{config.scenario}/bodem.tiff', config.partitionShape)

        # Create the hydraulic conductivity variable
        Ks       = dA.get.Ks(soilType)
        lfr.to_gdal(Ks, config.path + f'/output/{config.scenario}/Ks.tiff')
        
        # Create the land-use coefficient variable
        landC   = dA.get.landC(landUse)
        lfr.to_gdal(landC, config.path + f'/output/{config.scenario}/landUse_coefficients.tiff')
        
        # Initialize the static
        iniSurfaceWaterHeight, iniGroundWaterHeight = self.static(config.startDate, config.initialWaterTable, Ks, landC, landUse)
        print("static completed")

        self.iterate(config.startDate, config.endDate, iniSurfaceWaterHeight, iniGroundWaterHeight, Ks, landC)
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
    # Run the main model
    main = mainModel()
    main.simulate()
    
    # Process the results into a gif
    MakeGIF.main()
    
print("--- %s seconds ---" % (time.time() - start_time))