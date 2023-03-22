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
import dataAccess as dA
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
        self.dem        = lfr.from_gdal(config.path + f'/data/{config.scenario}/dem.tiff', config.partitionShape)           # DEM map of the study area
        landUse         = lfr.from_gdal(config.path + f'/data/{config.scenario}/landgebruik.tiff', config.partitionShape)        # Land-use, example: road
        soilType        = lfr.from_gdal(config.path + f'/data/{config.scenario}/bodem.tiff', config.partitionShape)              # example: sand or clay
        self.Ks         = dA.get.Ks(soilType)
        self.landC, self.mannings = dA.get.landC(landUse)
        self.ldd = lfr.d8_flow_direction(self.dem)
        # Load initial values for waterheight
        self.iniSurfaceWaterHeight = lfr.where(self.dem < config.initialWaterTable, config.initialWaterTable - self.dem, 0)
        self.iniGroundWaterHeight  = self.dem - config.waterBelowDEM
        for i in config.water:
            self.iniGroundWaterHeight = lfr.where(landUse == i, self.dem, self.iniGroundWaterHeight)

        # Reporting
        variables       = {"iniSurfaceWaterHeight": self.iniSurfaceWaterHeight}
        reporting.report.v2(config.startDate, 0, variables, config.output_path)


    @lfr.runtime_scope
    def dynamicModel(self):
        dt = config.dt              # Amount of small timesteps for routing in seconds
        dT = config.dT              # Amount of large timesteps for loading and saving data
        surfaceWaterHeight = self.iniSurfaceWaterHeight
        groundWaterHeight  = self.iniGroundWaterHeight
        discharge          = dG.generate.lue_zero() # Initial discharge through cell is zero (is speed of the water column in m/s)

        for i in range(int((config.endDate - config.startDate).days * dT)):
            # Timing and date
            currentDate = config.startDate + datetime.timedelta(minutes = i * dt/60)                # The actual date
            time = int(i * dt / 60)                                                                 # Time in minutes

            # Load data
            precipitation     = dA.get.precipitation(currentDate, dA.get.apiSession()) / 3600       # Divide by 3600 to convert from hour rate to second rate
            pot_evaporation   = dA.get.pot_evaporation(currentDate, dA.get.apiSession()) / 3600     # Divide by 3600 to convert from hour rate to second rate
            pot_infiltration  = dA.get.infiltration(self.dem, (self.dem - 0.1), self.Ks, self.landC)
            percolation       = dA.get.percolation(self.dem, groundWaterHeight, self.Ks)
            i_ratio, e_ratio  = dA.get.ieRatio(pot_evaporation, pot_infiltration)
            evaporation, infiltration = dA.get.EvaporationInfiltration(surfaceWaterHeight, pot_evaporation,\
                                                                       pot_infiltration, e_ratio, i_ratio)

            variables = {"precipitation": precipitation, "pot_evaporation": pot_evaporation, "pot_infiltration": pot_infiltration,\
                         "percolation": percolation}
            reporting.report.v2(currentDate, time, variables, config.output_path)

            # Loop
            for j in range(dt):
                # surfaceWaterHeight update based on vertical fluxes
                surfaceWaterHeight = discharge / (config.resolution ** 2)
                inflow             = precipitation - evaporation - infiltration
                inflow             = lfr.where(inflow > 0.00001, inflow, 0.00001) # Cannot be lower than zero

                # surfaceWaterRouting
                alpha              = 1.5
                beta               = 0.6
                channelLength      = dG.generate.lue_one()*config.resolution
                timestepduration   = 1.0
                discharge          = lfr.kinematic_wave(self.ldd, discharge, inflow,\
                                                        alpha, beta, timestepduration,\
                                                        channelLength,)

                lfr.maximum(precipitation).get()


            # Save / Report data
            print(f"Done: {i+1}/{dT}")
            variables = {"surfaceWaterHeight": surfaceWaterHeight,\
                         "discharge": discharge}
            reporting.report.v2(currentDate, time, variables, config.output_path)
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
    main.dynamicModel()

    # Process the results into a gif
    # MakeGIF.main()

print("--- %s seconds ---" % (time.time() - start_time))