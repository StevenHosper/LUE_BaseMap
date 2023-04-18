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
import numpy as np

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
        self.dem        = lfr.from_gdal(config.path + f'/data/{config.scenario}/v2/dem.tiff', config.partitionShape)               # DEM map of the study area
        self.dem        = lfr.where(self.dem < 0.1, 35, self.dem)
        landUse         = lfr.from_gdal(config.path + f'/data/{config.scenario}/landgebruik.tiff', config.partitionShape)       # Land-use, example: road
        soilType        = lfr.from_gdal(config.path + f'/data/{config.scenario}/bodem.tiff', config.partitionShape)             # example: sand or clay
        self.Ks, self.porosity, self.wiltingPoint = dA.get.soil_csv(config.soilData, soilType)                                  # soil characteristic
        self.mannings, self.permeability, self.interceptionStorageMax, self.throughfallFraction = dA.get.landCharacteristics_csv(config.landUseData, landUse)            # land-use characteristics
        
        self.groundwaterBase         = config.groundwaterBase
        
        # self.ldd = lfr.d8_flow_direction(self.dem)
        self.ldd        = lfr.from_gdal(config.path + f'/data/{config.scenario}/ldd_pcr_shaped.tiff', config.partitionShape)
        
        # iniGroundWaterHeight generated
        #self.iniGroundWaterHeight   = lfr.where(self.dem > self.groundwaterBase + config.waterBelowDEM, self.dem - config.waterBelowDEM, self.groundwaterBase)
        #self.iniGroundWaterHeight   = lfr.where(self.iniGroundWaterHeight > self.dem, self.dem, self.iniGroundWaterHeight)
        
        # iniGroundWaterHeight from memory
        self.iniGroundWaterHeight   = lfr.from_gdal(config.path + f'/data/{config.scenario}/v2/1_groundWaterHeight_2023-02-24_599.tiff', config.partitionShape)
        self.iniGroundWaterHeight   = lfr.where(self.iniGroundWaterHeight > self.dem, self.dem - config.waterBelowDEM, self.iniGroundWaterHeight)
        
        # Initial discharge
        self.iniDischarge           = lfr.from_gdal(config.path + f'/data/{config.scenario}/v2/1_discharge_2023-02-24_599.tiff', config.partitionShape) # Initial discharge through cell is zero (is speed of the water column in m/s)
        #self.iniDischarge           = dG.generate.lue_zero()
        
        # Initial InterceptionStorage and groundWaterStorage
        self.iniInterceptionStorage = dG.generate.lue_zero()
        self.iniGroundWaterStorage  = (self.iniGroundWaterHeight - (self.dem - config.impermeableLayerBelowDEM)) * self.cellArea
        
        # self.notBoundaryCells       = dG.generate.boundaryCell() # Currently not working
        
        # Set constants
        self.resolution             = config.resolution * dG.generate.lue_one()
        self.cellArea               = self.resolution * self.resolution
        self.slope          	    = lfr.slope(self.dem, config.resolution)
        
        


    @lfr.runtime_scope
    def dynamicModel(self):
        dt = config.dt              # Amount of small timesteps for routing in seconds
        dT = config.dT              # Amount of large timesteps for loading and saving data
        
        # Loading initial conditions
        groundWaterHeight   = self.iniGroundWaterHeight
        discharge           = self.iniDischarge
        Sgw                 = self.iniGroundWaterStorage
        interceptionStorage = self.iniInterceptionStorage
        
        # Setting min and max soil values
        MaxSgw              = config.impermeableLayerBelowDEM * self.cellArea       # Full storage of porosity
        MinSgw              = MaxSgw * (self.wiltingPoint / self.porosity)          # Minimum storage because of wilting point
        
        # Values for discharge to height calculation
        sqrtSlope           = lfr.sqrt(self.slope)
        sqrtSlope           = lfr.where(sqrtSlope < 0.001, 0.001, sqrtSlope)
        sqrtSlope           = lfr.where(sqrtSlope > 0.05, 0.05, sqrtSlope)
        width               = 1
        coefficient         = self.mannings / (sqrtSlope * width)
        
        # Channel length and area
        channelLength       = self.resolution
        channelArea         = width * channelLength
        
        # Kinematic Surface Water Routing Constants
        alpha               = 1.5
        beta                = 0.6
        timestepduration    = 1.0 * config.timestep

        # Static, really small value because inflow = 0 is not accepted
        inflow = dG.generate.lue_one()*0.000000000001
        
        # Current date
        date = config.startDate
        
        # Open file to write maximum discharge values to for post simulation validation.
        with open(config.output_path + "maximumDischarge.csv", "w", newline="") as f:
            writer = csv.writer(f, delimiter=';')
            
            # Start model for dT large periods
            for i in range(dT):
                # Time in minutes is the small iteration multiplied with the timestep (both in seconds) divided by 60 seconds.
                time = int((i * (dt*config.timestep)/60)) 
                
                # Load flux and storage values
                precipitation               = dA.get.precipitation(time, self.cellArea, dA.get.apiSession()) # m/s
                ref_evaporation             = dA.get.pot_evaporation(date, self.cellArea, dA.get.apiSession()) # m/s
                interceptionStorage, precipitation, evapotranspirationSurface = \
                                              dA.get.interception(self.cellArea, interceptionStorage, self.interceptionStorageMax, precipitation, \
                                                                  ref_evaporation, self.throughfallFraction)
                evapotranspirationSurface, evapotranspirationSoil = \
                                              dA.get.evapotranspiration(precipitation, evapotranspirationSurface, discharge)
                infiltration                = dA.get.pot_infiltration(Sgw, MaxSgw, self.cellArea, self.Ks, self.permeability, self.porosity, \
                                                                      discharge, precipitation, evapotranspirationSurface)


                # Groundwater LDD, gradient and flow flux
                gwLDD       = lfr.d8_flow_direction(groundWaterHeight)
                gwGradient  = (groundWaterHeight - lfr.downstream(gwLDD, groundWaterHeight)) / self.resolution
                Qgw         = self.Ks * gwGradient * config.timestep * self.resolution                        # Groundwater velocity in m/s
                
                # If the groundwater flow because of the impermeable layer is larger than the amount of water available, than it should be set so only the stored water will move.
                Qgw         = lfr.where(Qgw * config.dt > Sgw - MinSgw, Sgw - MinSgw, Qgw)
                Qgw         = lfr.where(Qgw < 0.000000000001, 0.000000000001, Qgw)
                
                # Add all vertical processes for the surfacewater and all processes groundwater
                gwFlux      = ((infiltration - evapotranspirationSoil)/self.porosity) + lfr.upstream(gwLDD, Qgw) - Qgw          # Is now in cubic meters
                swFlux      =  precipitation - evapotranspirationSurface - infiltration                                         # Is now in cubic meters
                
                for j in range(dt):
                    # The groundwater is adjusted by the fluxes
                    Sgw         = Sgw + gwFlux                                                                     
                    
                    height = lfr.pow(coefficient*discharge, 0.6)
                    
                    # If the groundwater table surpases the digital elevation map, groundwater is turned into runoff.
                    seepage     = lfr.where(Sgw > MaxSgw, (Sgw - MaxSgw)*self.porosity, 0)
                    
                    # Discharge is affected by the surfacewater fluxes, and seepage is added
                    height   = height + (swFlux + seepage)/channelArea
                    
                    discharge = lfr.pow(height, 1.6666666667) / coefficient
                    
                    # Because the kinematic wave has difficulties working with zero's, we have opted for a very small value. This will impact model results.
                    discharge   = lfr.where(discharge < 0.000000000001, 0.000000000001, discharge)

                    # Water routing based on the kinematic wave function, currently alpha is a float. Hopefully mannings raster can be used in the future.
                    discharge           = lfr.kinematic_wave(self.ldd, discharge, inflow,\
                                                alpha, beta, timestepduration,\
                                                channelLength,)
                    
                    # Any water that is moved from groundwater to discharge has to be removed from the groundwaterStorage
                    Sgw         = Sgw - (seepage/self.porosity)
                    
                    # Get the maximum value of the discharge raster (to limit the amount of tasks created by HPX)
                    Qmaxvalue = lfr.maximum(discharge).get()
                    print(Qmaxvalue)
                    
                    # Write value to csv for later validation
                    writer.writerow([i*60 + j, Qmaxvalue])
                
                # Make sure the storage of water can never be smaller than zero. There cannot be negative storage. Because we do not adjust fluxes every second,
                # but every 60 seconds it is possible a cell is depleted and becomes slightly negative.
                Sgw = lfr.where(Sgw <= 0, 0, Sgw)
                
                # Adjust the GW Table for the LDD creation of the next timestep.
                groundWaterHeight = self.dem - config.impermeableLayerBelowDEM + Sgw/self.cellArea
                
                # Save / Report data
                print(f"Done: {i+1}/{dT}")
                variables = {"discharge": discharge, "seepage": seepage, "Qgw": Qgw, "groundWaterHeight": groundWaterHeight, "Sgw": Sgw, "evapoSoil": evapotranspirationSoil, "infiltration": infiltration,
                             "gwFlux": gwFlux}
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
