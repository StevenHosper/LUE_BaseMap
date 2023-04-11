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
import csv

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
        self.dem        = lfr.from_gdal(config.path + f'/data/{config.scenario}/dem3.tif', config.partitionShape)               # DEM map of the study area
        self.dem        = lfr.where(self.dem < 0.1, 35, self.dem)
        self.landUse    = lfr.from_gdal(config.path + f'/data/{config.scenario}/landgebruik.tiff', config.partitionShape)       # Land-use, example: road
        soilType        = lfr.from_gdal(config.path + f'/data/{config.scenario}/bodem.tiff', config.partitionShape)             # example: sand or clay
        self.Ks, self.porosity    = dA.get.soil_csv(config.soilData, soilType)                                                  # soil characteristic
        self.mannings, self.permeability, self.interceptionStorageMax, self.throughfallFraction = dA.get.landCharacteristics_csv(config.landUseData, soilType)            # land-use characteristics
        
        self.groundwaterBase         = config.groundwaterBase
        
        # self.ldd = lfr.d8_flow_direction(self.dem)
        self.ldd        = lfr.from_gdal(config.path + f'/data/{config.scenario}/ldd_pcr_shaped.tiff', config.partitionShape)
        
        # Load initial values for waterheight
        self.iniSurfaceWaterHeight  = dG.generate.lue_zero()
        
        # iniGroundWaterHeight generated
        self.iniGroundWaterHeight   = lfr.where(self.dem > self.groundwaterBase + config.waterBelowDEM, self.dem - config.waterBelowDEM, self.groundwaterBase)
        self.iniGroundWaterHeight   = lfr.where(self.iniGroundWaterHeight > self.dem, self.dem, self.iniGroundWaterHeight)
        
        # iniGroundWaterHeight from memory
        # self.iniGroundWaterHeight   = lfr.from_gdal(config.path + f'/data/{config.scenario}/1_groundWaterHeight_2023-02-24_10hours.tiff', config.partitionShape)
        
        self.iniGroundWaterStorage  = self.iniGroundWaterHeight - (self.dem - config.imperviousLayerBelowDEM)
        #self.iniDischarge           = lfr.from_gdal(config.path + f'/data/{config.scenario}/1_discharge_2023-02-24_10hours.tiff', config.partitionShape) / 5 # Initial discharge through cell is zero (is speed of the water column in m/s)
        self.iniDischarge           = dG.generate.lue_zero()
        # self.notBoundaryCells       = dG.generate.boundaryCell() # Currently not working
        self.resolution             = config.resolution * dG.generate.lue_one()
        self.cellArea               = self.resolution * self.resolution
        


    @lfr.runtime_scope
    def dynamicModel(self):
        dt = config.dt              # Amount of small timesteps for routing in seconds
        dT = config.dT              # Amount of large timesteps for loading and saving data
        
        surfaceWaterHeight  = self.iniSurfaceWaterHeight
        groundWaterHeight   = self.iniGroundWaterHeight
        discharge           = self.iniDischarge
        Sgw                 = self.iniGroundWaterStorage
        MaxSgw              = config.imperviousLayerBelowDEM
        interceptionStorage = dG.generate.lue_zero()
        
        date = config.startDate
        
        Qmax = []
        
        # Static, really small value because inflow = 0 is not accepted
        inflow = dG.generate.lue_one()*0.000000000001
        with open(config.output_path + "maximumDischarge.csv", "w", newline="") as f:
            writer = csv.writer(f, delimiter=';')
            for i in range(dT):
                time = int((i * (dt*config.timestep)/60)) 
                
                # Load values
                precipitation               = dA.get.precipitation(date, self.cellArea, dA.get.apiSession())
                ref_evaporation             = dA.get.pot_evaporation(date, self.cellArea, dA.get.apiSession())
                interceptionStorage, interception, precipitation, evapotranspirationSurface = \
                                              dA.get.interception(self.cellArea, interceptionStorage, self.interceptionStorageMax, precipitation, \
                                                                  ref_evaporation, self.throughfallFraction)
                pot_infiltration            = dA.get.infiltration(self.dem, self.cellArea, groundWaterHeight, self.Ks, self.permeability, self.porosity)
                i_ratio, e_ratio            = dA.get.ieRatio(ref_evaporation, pot_infiltration)
                evaporation, infiltration   = dA.get.EvaporationInfiltration(precipitation, surfaceWaterHeight, ref_evaporation, pot_infiltration, e_ratio, i_ratio)
                
                # Groundwater hydraulic gradients and corresponding discharge
                gwLDD       = lfr.d8_flow_direction(groundWaterHeight)
                gwGradient  = (groundWaterHeight - lfr.downstream(gwLDD, groundWaterHeight)) / self.resolution
                Qgw         = Sgw * self.Ks * gwGradient * config.timestep * self.cellArea               # Groundwater velocity in m/s
                
                for j in range(dt):
                    Sgw         = Sgw + (infiltration * self.porosity) + lfr.upstream(gwLDD, Qgw) - Qgw
                    seepage     = lfr.where(Sgw > MaxSgw, (Sgw - MaxSgw)*self.porosity, 0)
                    discharge   = discharge + precipitation - evaporation - infiltration + seepage
                    discharge   = lfr.where(discharge < 0.000000000001, 0.000000000001, discharge)
                    
                    # Kinematic Surface Water Routing 
                    alpha               = 1.5
                    beta                = 0.6
                    channelLength       = self.resolution
                    timestepduration    = 1.0 * config.timestep
                    discharge           = lfr.kinematic_wave(self.ldd, discharge, inflow,\
                                                alpha, beta, timestepduration,\
                                                channelLength,)
                    
                    Sgw         = Sgw - (seepage/self.porosity)
                    
                    Qmaxvalue = lfr.maximum(discharge).get()
                    print(Qmaxvalue)
                    writer.writerow([i*60 + j, Qmaxvalue])
                
                Sgw = lfr.where(Sgw <= 0, 0, Sgw)
                groundWaterHeight = self.dem - config.imperviousLayerBelowDEM + Sgw
                
                # Save / Report data
                print(f"Done: {i+1}/{dT}")
                
                variables = {"discharge": discharge, "groundWaterHeight": groundWaterHeight, "Sgw": Sgw, "Qgw": Qgw}
                reporting.report.v2(date, time, variables, config.output_path)
        
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
