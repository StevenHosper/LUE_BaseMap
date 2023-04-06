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

# Own functions
# TO-DO: Reformat by the use of: from ... import ... as ...
import configuration as config
import dataAccess2 as dA
import MakeGIF
import reporting
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
        # Initialize data required from memory files
        # Get all constants
        self.dem        = lfr.from_gdal(config.path + f'/data/{config.scenario}/dem.tiff', config.partitionShape)               # DEM map of the study area
        self.dem        = lfr.where(self.dem < 0.1, 35, self.dem)
        self.landUse    = lfr.from_gdal(config.path + f'/data/{config.scenario}/landgebruik.tiff', config.partitionShape)       # Land-use, example: road
        soilType        = lfr.from_gdal(config.path + f'/data/{config.scenario}/bodem.tiff', config.partitionShape)             # example: sand or clay
        self.Ks, self.porosity    = dA.get.soil_csv(config.soilData, soilType)                                                  # soil characteristic
        self.mannings, self.permeability, self.interceptionStorageMax, self.throughfallFraction = dA.get.landCharacteristics_csv(config.landUseData, soilType)            # land-use characteristics
        
        self.imperviousLayer         = config.imperviousLayer
        
        # self.ldd = lfr.d8_flow_direction(self.dem)
        self.ldd        = lfr.from_gdal(config.path + f'/data/{config.scenario}/ldd_pcr_shaped.tiff', config.partitionShape)
        
        # Load initial values for waterheight
        self.iniSurfaceWaterHeight  = lfr.where(self.dem < config.initialWaterTable, config.initialWaterTable - self.dem, 0)
        self.iniGroundWaterHeight   = lfr.where(self.dem < self.imperviousLayer, self.dem - config.waterBelowDEM, self.imperviousLayer + ((self.dem - self.imperviousLayer)/2))
        # self.iniGroundWaterStorage  = self.iniGroundWaterHeight - self.imperviousLayer             # The height between the impervious bottom layer and the top of the groundwater table is the amount of water stored.
        self.iniDischarge           = lfr.from_gdal(config.path + f'/data/HupselOutput/1_discharge_2023-02-24_80.tiff', config.partitionShape) / 600 # Initial discharge through cell is zero (is speed of the water column in m/s)
        
        # self.notBoundaryCells       = dG.generate.boundaryCell() # Currently not working
        self.resolution             = config.resolution * dG.generate.lue_one()
        # self.cellArea               = self.resolution * self.resolution
        


    @lfr.runtime_scope
    def dynamicModel(self):
        dt = config.dt              # Amount of small timesteps for routing in seconds
        dT = config.dT              # Amount of large timesteps for loading and saving data
        
        surfaceWaterHeight  = self.iniSurfaceWaterHeight
        groundWaterHeight   = self.iniGroundWaterHeight
        discharge           = self.iniDischarge
        # Sgw                 = self.iniGroundWaterStorage
        # interceptionStorage    = dG.generate.lue_zero()
        
        date = config.startDate
        
        
        # Static, really small value because inflow = 0 is not accepted
        inflow = dG.generate.lue_one()*0.00000000001
        
        for i in range(dT):
            time = int((i * (dt*config.timestep)/60)) 
            
            # Load values
            precipitation               = dA.get.precipitation(date, dA.get.apiSession())
            ref_evaporation             = dA.get.pot_evaporation(date, dA.get.apiSession())
            pot_infiltration            = dA.get.infiltration(self.dem, groundWaterHeight, self.Ks, self.permeability, self.porosity)
            i_ratio, e_ratio            = dA.get.ieRatio(ref_evaporation, pot_infiltration)
            evaporation, infiltration   = dA.get.EvaporationInfiltration(precipitation, surfaceWaterHeight, ref_evaporation, pot_infiltration, e_ratio, i_ratio)
            # interception, precipitation = dA.get.interception(interceptionStorage, self.interceptionStorageMax, precipitation, self.throughfallFraction)
            # precipitation, interceptionEvaporation, surfaceEvaporation = dA.get.evapotranspiration(precipitation, ref_evaporation, interception, interceptionStorage)
            
            # Groundwater hydraulic gradients and corresponding discharge
            # gwLDD       = lfr.d8_flow_direction(groundWaterHeight)
            # gwGradient  = (groundWaterHeight - lfr.downstream(gwLDD, groundWaterHeight)) / self.resolution
            # Qgw         = Sgw * self.porosity * self.Ks * gwGradient * config.timestep                # Groundwater velocity in m/s
            
            for j in range(dt):
                discharge  = discharge + precipitation - evaporation - infiltration
                discharge  = lfr.where(discharge < 0.00000000001, 0.00000000001, discharge)
                
                # Kinematic Surface Water Routing 
                alpha               = 1.5
                beta                = 0.6
                channelLength       = self.resolution
                timestepduration    = 1.0 * config.timestep
                discharge           = lfr.kinematic_wave(self.ldd, discharge, inflow,\
                                               alpha, beta, timestepduration,\
                                               channelLength,)
            
            # Save / Report data
            print(f"Done: {i+1}/{dT}")
            
            variables = {"infiltration": infiltration, "discharge": discharge,}
            reporting.report.v2(date, time, variables, config.output_path)
            
            input("Press Enter to continue...")
            
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
    "hpx.diagnostics_on_terminate!=1",
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
    main.dynamicModel()

    # Process the results into a gif
    # MakeGIF.main()

print("--- %s seconds ---" % (time.time() - start_time))
