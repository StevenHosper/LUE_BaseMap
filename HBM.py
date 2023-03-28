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
        self.landUse    = lfr.from_gdal(config.path + f'/data/{config.scenario}/landgebruik.tiff', config.partitionShape)        # Land-use, example: road
        soilType        = lfr.from_gdal(config.path + f'/data/{config.scenario}/bodem.tiff', config.partitionShape)              # example: sand or clay
        self.Ks, self.porosity    = dA.get.soil(soilType)
        self.landC, self.mannings = dA.get.landC(self.landUse)
        # self.ldd = lfr.d8_flow_direction(self.dem)
        self.ldd        = lfr.from_gdal(config.path + f'/data/{config.scenario}/ldd_pcr_shaped.tiff', config.partitionShape)
        # Load initial values for waterheight
        self.iniSurfaceWaterHeight  = lfr.where(self.dem < config.initialWaterTable, config.initialWaterTable - self.dem, 0)
        self.iniGroundWaterStorage  = dG.generate.lue_one()*(config.imperviousLayer - config.waterBelowDEM)             # The height between the impervious bottom layer and the top of the groundwater table is the amount of water stored.
        self.iniGroundWaterHeight   = self.dem - config.waterBelowDEM

        # Reporting
        variables       = {"iniSurfaceWaterHeight": self.iniSurfaceWaterHeight}
        reporting.report.v2(config.startDate, 0, variables, config.output_path)


    @lfr.runtime_scope
    def dynamicModel(self):
        dt = config.dt              # Amount of small timesteps for routing in seconds
        dT = config.dT              # Amount of large timesteps for loading and saving data
        surfaceWaterHeight = self.iniSurfaceWaterHeight
        groundWaterHeight  = self.iniGroundWaterHeight
        #discharge          = dG.generate.lue_zero()
        discharge          = lfr.from_gdal(config.path + f'/data/{config.scenario}/discharge_1mm1t2h40.tiff', config.partitionShape) # Initial discharge through cell is zero (is speed of the water column in m/s)
        interceptionStorageMax = dA.get.interceptionStorageMax(landUse=self.landUse)
        interceptionStorage    = dG.generate.lue_zero()
        Sgw                = self.iniGroundWaterStorage                                             # In groundwaterheight in meters (not accounting for porosity)
        
        for i in range(int((config.endDate - config.startDate).days * dT)):
            # Timing and date
            currentDate = config.startDate + datetime.timedelta(minutes = i * dt/60)                # The actual date
            time = int(i * dt/60)                                                                   # Time in minutes

            # Load data
            precipitation     = dA.get.precipitation(currentDate, dA.get.apiSession())
            pot_evaporation   = dA.get.pot_evaporation(currentDate, dA.get.apiSession())
            pot_infiltration  = dA.get.infiltration(self.dem, (self.dem - 0.1), self.Ks, self.landC)
            percolation       = dA.get.percolation(self.dem, groundWaterHeight, self.Ks)
            i_ratio, e_ratio  = dA.get.ieRatio(pot_evaporation, pot_infiltration)
            evaporation, infiltration = dA.get.EvaporationInfiltration(surfaceWaterHeight, pot_evaporation,\
                                                                       pot_infiltration, e_ratio, i_ratio)

            variables = {"precipitation": precipitation, "pot_evaporation": pot_evaporation, "pot_infiltration": pot_infiltration,\
                         "percolation": percolation}
            reporting.report.v2(currentDate, time, variables, config.output_path)

            gwLDD = lfr.d8_flow_direction(groundWaterHeight)
            # gwGradient = lfr.slope(groundWaterHeight, config.resolution)
            gwGradient = (groundWaterHeight - lfr.downstream(gwLDD, groundWaterHeight)) / config.resolution
            Qgw        = lfr.where(dG.generate.boundaryCell(), Sgw * self.porosity * self.Ks * gwGradient, 0)                                                # Groundwater velocity in m/s
            variables  = {"gwGradient": gwGradient, "Sgw": Sgw, "Qgw": Qgw}
            reporting.report.v2(currentDate, time, variables, config.output_path) 
                       
            # Loop
            for j in range(dt):
                # surfaceWaterHeight update based on vertical fluxes
                surfaceWaterHeight = discharge / (config.resolution ** 2)
                surfaceWaterHeight = surfaceWaterHeight - evaporation - infiltration
                surfaceWaterHeight = lfr.where(surfaceWaterHeight > 0, surfaceWaterHeight, 0)       # Once surface water height reaches zero, it cannot be lower
                discharge          = surfaceWaterHeight * config.resolution ** 2                    # Convert surface water height back to discharge
                # interception, precipitation = dA.get.interception(interceptionStorage, interceptionStorageMax,\
                #                                                   precipitation, self.landUse)
                
                # Groundwater storage
                Sgw        = lfr.where(dG.generate.boundaryCell(), Sgw + (infiltration/self.porosity) - Qgw/self.porosity + lfr.upstream(gwLDD, Qgw)/self.porosity, Sgw) # Groundwater gets the infiltration
                seepage    = lfr.where(Sgw > (config.imperviousLayer), (Sgw - config.imperviousLayer)*self.porosity, 0)             # The excess groundwater height * porosity is the outflow
                Sgw        = lfr.where(Sgw > (config.imperviousLayer), config.imperviousLayer, Sgw)                                 # After the seepage the groundwater is back to the maximum groundwater height
                
                inflow             = precipitation + seepage
                inflow             = lfr.where(inflow > 0.00000001, inflow, 0.00000001) # Cannot be lower than zero

                # surfaceWaterRouting
                alpha              = 1.5
                beta               = 0.6
                channelLength      = dG.generate.lue_one()*config.resolution
                timestepduration   = 1.0
                discharge          = lfr.kinematic_wave(self.ldd, discharge, inflow,\
                                                        alpha, beta, timestepduration,\
                                                        channelLength,)

                lfr.maximum(precipitation).get()

            # Groundwaterheight map
            groundWaterHeight = self.dem - config.imperviousLayer + Sgw
            
            # Save / Report data
            print(f"Done: {i+1}/{dT}")
            variables = {"surfaceWaterHeight": surfaceWaterHeight, "inflow": inflow,\
                         "groundWaterHeight": groundWaterHeight, "discharge": discharge,}
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