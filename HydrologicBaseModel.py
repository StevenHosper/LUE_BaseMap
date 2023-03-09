# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 2023

@author: steven.hosper
"""
# The HydrologicBaseModel
import lue.framework as lfr
import lue.pcraster as lpr
import numpy as np
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
        self.waterheight = lfr.where(self.dem < initialWaterTable, initialWaterTable - self.dem, 0)
        self.waterheight = lfr.where(landUse == 51, self.waterheight, 0 )
        
        # Create a groundwater table 
        # TO-DO : In the future we would preferable read the groundwater table from measuring stations \
        #         and apply these. However, that is not yet possible.
        self.groundWaterHeight = self.dem - 0.1
        self.groundWaterHeight = lfr.where(landUse == 51, self.dem, self.groundWaterHeight)
        
        lfr.to_gdal(self.waterheight, config.path + f'/output/{config.scenario}/initial_water_height.tiff')

        # Access the data from the directory or the API dependend on the settings
        precipitation     = dA.get.precipitation(current_date, dA.get.apiSession(), dG.generate.lue_zero())
        pot_evaporation   = dA.get.pot_evaporation(current_date, dA.get.apiSession(), dG.generate.lue_zero())   
        pot_infiltration  = dA.get.infiltration(self.dem, self.groundWaterHeight, Ks, landC, dG.generate.lue_zero())
        percolation, part = dA.get.percolation(self.dem, self.groundWaterHeight, Ks, dG.generate.lue_zero())
        i_ratio, e_ratio  = dA.get.ieRatio(pot_evaporation, pot_infiltration, dG.generate.lue_one(), dG.generate.lue_zero())
        
        # Add precipitation to the watertable
        self.waterheight = self.waterheight + precipitation

        # Check if the waterheight is more than the evaporation and infiltration combined
        ie = pot_evaporation + pot_infiltration
        
        # Calculate the actual evaporation and infiltration
        evaporation  = lfr.where(self.waterheight >= ie, pot_evaporation, self.waterheight*e_ratio)
        infiltration = lfr.where(self.waterheight >= ie, pot_infiltration, self.waterheight*i_ratio)
        
        # Runoff test
        runoff = update.update.runoff(precipitation, evaporation, infiltration)
        
        self.ldd = lfr.where(runoff >= -100, self.ldd, 5)
        
        discharge = lfr.kinematic_wave(self.ldd, runoff, dG.generate.lue_one()*0.0001, 1.5, 0.6, 8000, dG.generate.lue_one())
        
        # Remove the evaporation and infiltration from the waterheight as it is lost to the atmosphere or groundwater.
        # *Note --> the rain now happens 'first', depending on when the rain falls during the day there might not be time to evaporate, but this is currently not taken into account.
        self.waterheight = self.waterheight - evaporation - infiltration
        self.groundWaterHeight = self.groundWaterHeight + infiltration - percolation
        
        # Make sure the waterheight is not below zero as this is not possible.
        # There can be no evaporation without water.
        self.waterheight = lfr.where(self.waterheight < dG.generate.lue_zero(), dG.generate.lue_zero(), self.waterheight)
        self.height = self.waterheight + self.dem
        
        # Create file with current situation of the water balance
        variables = {"waterheight": self.waterheight, "groundwater": self.groundWaterHeight, "runoff": runoff, "discharge": discharge}
        reporting.report.static(current_date, variables, config.output_path)
        return 0
    
    

    def iterate(self, start_date: datetime.date, end_date: datetime.date, Ks, landC):
        for i in range(int((end_date - start_date).days)):
            # Print the time to keep track while the program runs
            print(f'The date is: {start_date + datetime.timedelta(1+i)}')
            current_date = start_date + datetime.timedelta(1+i)
            
            # Recalculate if the water should still flow in the same direction
            ldd = lfr.d8_flow_direction(self.height)
            
            # Access the data from the directory or the API dependend on the settings
            precipitation     = dA.get.precipitation(current_date, dA.get.apiSession(), dG.generate.lue_zero())
            pot_evaporation   = dA.get.pot_evaporation(current_date, dA.get.apiSession(), dG.generate.lue_zero())
            pot_infiltration  = dA.get.infiltration(self.dem, self.groundWaterHeight, Ks, landC, dG.generate.lue_zero())
            percolation, part = dA.get.percolation(self.dem, self.groundWaterHeight, Ks, dG.generate.lue_zero())
            i_ratio, e_ratio  = dA.get.ieRatio(pot_evaporation, pot_infiltration, dG.generate.lue_one(), dG.generate.lue_zero())

            lfr.to_gdal(part, config.path + f'/output/{config.scenario}/percolation_{current_date}.tiff')
            
            # Add precipitation to the watertable
            self.waterheight = self.waterheight + precipitation
            
            # Check if the waterheight is more than the evaporation and infiltration combined
            ie = pot_evaporation + pot_infiltration
            
            # Calculate the actual evaporation and infiltration
            evaporation  = lfr.where(self.waterheight >= ie, pot_evaporation, self.waterheight*e_ratio)
            infiltration = lfr.where(self.waterheight >= ie, pot_infiltration, self.waterheight*i_ratio)
            
            # Check the difference in dem (as this determines the total height that should be filled to create an equal surface)
            height_difference = self.height - lfr.downstream(ldd, self.height)
            potential_flux = lfr.where(height_difference > self.waterheight, 0.5*self.waterheight, 0.5*height_difference)
            flux = lfr.where(ldd != 5, potential_flux, 0)
            
            # ROUTING
            self.waterheight = self.waterheight + lfr.upstream(ldd, flux) - flux
            
            self.groundWaterHeight = update.update.groundWaterHeight(
                self.dem, Ks, self.waterheight, self.groundWaterHeight, infiltration, \
                percolation, dG.generate.lue_zero()
                )

            # Remove the evaporation and infiltration from the waterheight as it is lost to the 
            # atmosphere or groundwater.
            self.waterheight = self.waterheight - evaporation - infiltration
            self.groundWaterHeight = self.groundWaterHeight + infiltration - percolation
            
            runoff = update.update.runoff(precipitation, evaporation, infiltration)
            discharge = lfr.kinematic_wave(self.ldd, runoff, dG.generate.lue_one()*0.0001, 1.5, 0.6, 8000, dG.generate.lue_one())
            
            # Waterheight can never be lower than zero.
            self.waterheight = lfr.where(self.waterheight < dG.generate.lue_zero(), dG.generate.lue_zero(), self.waterheight)
            
            # Adjust the concurrent height
            self.height = self.dem + self.waterheight
            
            # Create file with current situation of the water balance
            variables = {"waterheight": self.waterheight, "groundwater": self.groundWaterHeight, "runoff": runoff, "discharge": discharge}
            fluxes    = {}
            reporting.report.dynamic(current_date, variables, fluxes, config.output_path)
        return 0
     
     
    @lfr.runtime_scope 
    def simulate(self):
        # Load initial variables
        self.dem  = lfr.from_gdal(config.path + f'/data/{config.scenario}/dem.tiff', config.partitionShape)
        landUse  = lfr.from_gdal(config.path + f'/data/{config.scenario}/landgebruik.tiff', config.partitionShape)
        soilType = lfr.from_gdal(config.path + f'/data/{config.scenario}/bodem.tiff', config.partitionShape)

        # Attempt with pcraster to make ldd
        pcr.setclone(config.arrayExtent, config.arrayExtent, config.resolution, 0, 0)
        pcr.setglobaloption("lddin")
        ds = gdal.Open(config.path + f'/data/{config.scenario}/dem.tiff')
        raster = ds.GetRasterBand(1)
        a = raster.ReadAsArray()
        result = pcr.numpy2pcr(pcr.Scalar, a, 999)
        ldd = pcr.lddcreate(result, 9999999, 500, 9999999, 9999999)
        pcr.report(ldd, config.path + f'/output/{config.scenario}/ldd_pcr.tiff')
        
        self.ldd = lfr.from_gdal(config.path + f'/output/{config.scenario}/ldd_pcr.tiff', config.partitionShape)

        # Create the hydraulic conductivity variable
        Ks       = dA.get.Ks(soilType, dG.generate.lue_one())
        lfr.to_gdal(Ks, config.path + f'/output/{config.scenario}/Ks.tiff')
        
        # Create the land-use coefficient variable
        landC   = dA.get.landC(landUse, dG.generate.lue_one())
        lfr.to_gdal(landC, config.path + f'/output/{config.scenario}/landUse_coefficients.tiff')
        
        # Initialize the static
        self.static(config.startDate, config.initialWaterTable, Ks, landC, landUse)
        print("static completed")

        self.iterate(config.startDate, config.endDate, Ks, landC)
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