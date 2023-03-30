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
        self.dem        = lfr.where(self.dem < 0.1, 35, self.dem)
        self.landUse    = lfr.from_gdal(config.path + f'/data/{config.scenario}/landgebruik.tiff', config.partitionShape)        # Land-use, example: road
        soilType        = lfr.from_gdal(config.path + f'/data/{config.scenario}/bodem.tiff', config.partitionShape)              # example: sand or clay
        self.Ks, self.porosity    = dA.get.soil(soilType)
        self.landC, self.mannings, self.soilC = dA.get.landC(self.landUse)
        
        # self.ldd = lfr.d8_flow_direction(self.dem)
        self.ldd        = lfr.from_gdal(config.path + f'/data/{config.scenario}/ldd_pcr_shaped.tiff', config.partitionShape)
        
        # Load initial values for waterheight
        self.iniSurfaceWaterHeight  = lfr.where(self.dem < config.initialWaterTable, config.initialWaterTable - self.dem, 0)
        self.iniGroundWaterHeight   = lfr.where(self.dem < config.imperviousLayer, self.dem - config.waterBelowDEM, config.imperviousLayer + ((self.dem - config.imperviousLayer)/2))
        self.iniGroundWaterStorage  = self.iniGroundWaterHeight - config.imperviousLayer             # The height between the impervious bottom layer and the top of the groundwater table is the amount of water stored.
        #self.iniGroundWaterHeight = lfr.from_gdal(config.path + f'/data/HupselOutput/1_groundWaterHeight_2023-02-24_80.tiff', config.partitionShape)
        #self.iniDischarge         = dG.generate.lue_zero()
        self.iniDischarge         = lfr.from_gdal(config.path + f'/data/HupselOutput/1_discharge_2023-02-24_80.tiff', config.partitionShape) / 600 # Initial discharge through cell is zero (is speed of the water column in m/s)
        
        # Reporting
        # variables       = {"iniSurfaceWaterHeight": self.iniSurfaceWaterHeight, "iniGroundWaterHeight": self.iniGroundWaterHeight}
        # reporting.report.v2(config.startDate, 0, variables, config.output_path)
        


    @lfr.runtime_scope
    def dynamicModel(self):
        dt = config.dt              # Amount of small timesteps for routing in seconds
        dT = config.dT              # Amount of large timesteps for loading and saving data
        
        surfaceWaterHeight = self.iniSurfaceWaterHeight
        groundWaterHeight  = self.iniGroundWaterHeight
        discharge          = self.iniDischarge
        Sgw                = self.iniGroundWaterStorage                                         # In groundwaterheight in meters (not accounting for porosity)
        
        interceptionStorageMax = dA.get.interceptionStorageMax(landUse=self.landUse)
        interceptionStorage    = dG.generate.lue_zero()
        
        for i in range(int((config.endDate - config.startDate).days * dT)):
            # Timing and date
            currentDate = config.startDate + datetime.timedelta(seconds = i * dt * config.timestep)                # The actual date
            time = int((i * (dt*config.timestep)/60))                                                                # Time in minutes

            # Load data
            precipitation     = dA.get.precipitation(currentDate, dA.get.apiSession())      * config.timestep
            ref_evaporation   = dA.get.pot_evaporation(currentDate, dA.get.apiSession())    * config.timestep
            pot_infiltration  = dA.get.infiltration(self.dem, Sgw, self.Ks, self.landC, self.porosity)* config.timestep
            percolation       = dA.get.percolation(self.dem, groundWaterHeight, self.Ks)    * config.timestep
            i_ratio, e_ratio  = dA.get.ieRatio(ref_evaporation, pot_infiltration)
            evaporation, infiltration = dA.get.EvaporationInfiltration(surfaceWaterHeight, ref_evaporation,\
                                                                       pot_infiltration, e_ratio, i_ratio)

            # variables = {"precipitation": precipitation, "evaporation": pot_evaporation, "infiltration": pot_infiltration,\
            #             }
            # reporting.report.v2(currentDate, time, variables, config.output_path)

            gwLDD = lfr.d8_flow_direction(groundWaterHeight)
            # gwGradient = lfr.slope(groundWaterHeight, config.resolution)
            gwGradient = lfr.where(dG.generate.boundaryCell(), (groundWaterHeight - lfr.downstream(gwLDD, groundWaterHeight)) / config.resolution, 0)
            Qgw        = (Sgw - config.imperviousLayer) * self.porosity * self.Ks * gwGradient * config.timestep                # Groundwater velocity in m/s
            # variables  = {"gwGradient": gwGradient, "Sgw": Sgw, "Qgw": Qgw}
            # reporting.report.v2(currentDate, time, variables, config.output_path) 
                       
            # Loop
            for j in range(dt):
                # surfaceWaterHeight update based on vertical fluxes
                surfaceWaterHeight = discharge / (config.resolution ** 2)
                
                # interception, precipitation = dA.get.interception(interceptionStorage, interceptionStorageMax,\
                #                                                   precipitation, self.landUse)
                
                # Groundwater storage
                Sgw        = lfr.where(dG.generate.boundaryCell(),\
                             Sgw + (infiltration/self.porosity) - Qgw/self.porosity + lfr.upstream(gwLDD, Qgw)/self.porosity \
                                 - ((Sgw - config.imperviousLayer)*0.8) * ref_evaporation * self.soilC,\
                             Sgw) # Groundwater gets the infiltration
                seepage    = lfr.where((Sgw + config.imperviousLayer) > self.dem, ((Sgw + config.imperviousLayer) - self.dem)*self.porosity, 0)             # The excess groundwater height * porosity is the outflow
                Sgw        = lfr.where((Sgw + config.imperviousLayer) > self.dem, self.dem - config.imperviousLayer, Sgw)                                 # After the seepage the groundwater is back to the maximum groundwater height
                
                surfaceWaterHeight = surfaceWaterHeight - evaporation - infiltration + precipitation + seepage
                surfaceWaterHeight = lfr.where(surfaceWaterHeight > 0, surfaceWaterHeight, 0)       # Once surface water height reaches zero, it cannot be lower
                discharge          = surfaceWaterHeight * config.resolution ** 2                    # Convert surface water height back to discharge
                
                inflow             = dG.generate.lue_one()*0.0000000001

                # surfaceWaterRouting
                alpha              = 1.5
                beta               = 0.6
                channelLength      = dG.generate.lue_one()*config.resolution
                timestepduration   = 1.0 * config.timestep
                discharge          = lfr.kinematic_wave(self.ldd, discharge, inflow,\
                                                        alpha, beta, timestepduration,\
                                                        channelLength,)
                lfr.maximum(precipitation).get()

            # Groundwaterheight map
            groundWaterHeight = self.dem - config.imperviousLayer + Sgw
            
            # Save / Report data
            print(f"Done: {i+1}/{dT}")
            variables = {"infiltration": infiltration, "groundWaterHeight": groundWaterHeight,\
                         "discharge": discharge, "Qgw": Qgw}
            reporting.report.v2(currentDate, time, variables, config.output_path)
            
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