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
from configuration_v2 import Configuration
import dataAccess2 as dA
import MakeGIF
from reporting import report
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
    def __init__(self, configuration):
        # Set directories
        self.inputDir   = configuration.generalSettings['inputDir'] + configuration.generalSettings['scenario'] 
        self.outputDir  = configuration.generalSettings['outputDir'] + configuration.generalSettings['scenario']
        
        partitionShape  = 2 * (configuration.modelSettings['partitionExtent'],)
        arrayShape      = 2 * (configuration.modelSettings['arrayExtent'],)
        
        # Initialize data required from memory files
        # Get all constants        
        self.dem        = lfr.from_gdal(self.inputDir + configuration.dataSettings['dem'], partitionShape)               # DEM map of the study area
        self.dem        = lfr.where(self.dem < 0.1, 35, self.dem)
        landUse         = lfr.from_gdal(self.inputDir + configuration.dataSettings['landUseMap'], partitionShape)       # Land-use, example: road
        soilType        = lfr.from_gdal(self.inputDir + configuration.dataSettings['soilMap'], partitionShape)             # example: sand or clay
        self.Ks, self.porosity, self.wiltingPoint = dA.get.soil_csv(configuration.generalSettings['inputDir'] + configuration.dataSettings['soilData'], soilType)                                  # soil characteristic
        self.mannings, self.permeability, self.interceptionStorageMax, self.throughfallFraction = \
            dA.get.landCharacteristics_csv(configuration.generalSettings['inputDir'] + configuration.dataSettings['landUseData'],
                                           landUse)            # land-use characteristics
        
        self.groundWaterBase         = float(configuration.modelSettings['groundWaterBase']) * dG.generate.lue_one()
        
        # self.ldd = lfr.d8_flow_direction(self.dem)
        self.ldd        = lfr.from_gdal(self.inputDir + configuration.dataSettings['ldd'], partitionShape)
        
        # Set constants
        self.resolution                 = float(configuration.modelSettings['resolution'])
        self.cellArea                   = self.resolution * self.resolution
        self.slope          	        = lfr.slope(self.dem, self.resolution)
        self.impermeableLayerBelowDEM   = float(configuration.modelSettings['impermeableLayerBelowDEM'])
        self.impermeableLayerHeight     = self.dem - self.impermeableLayerBelowDEM
        self.waterBelowDEM              = float(configuration.modelSettings['waterBelowDEM'])
        # self.notBoundaryCells       = dG.generate.boundaryCell() # Currently not working

        # Load initial groundWaterHeight, if no raster is supplied, use the waterBelowDEM in combination with DEM to create a initialGroundWaterHeight layer.
        try:
            self.iniGroundWaterHeight   = lfr.from_gdal(self.inputDir + configuration.dataSettings['iniGroundWaterHeight'], partitionShape)
            self.iniGroundWaterHeight   = lfr.where(self.iniGroundWaterHeight > self.dem, self.dem - self.waterBelowDEM, self.iniGroundWaterHeight)
        except:
            self.iniGroundWaterHeight   = lfr.where(self.dem > self.groundWaterBase + self.waterBelowDEM, self.dem - self.waterBelowDEM,
                                                    self.groundWaterBase)
            self.iniGroundWaterHeight   = lfr.where(self.iniGroundWaterHeight > self.dem, self.dem, self.iniGroundWaterHeight)
        
        # Load initial discharge, if no raster is supplied, set to zero.
        try:
            self.iniDischarge           = lfr.from_gdal(self.inputDir + configuration.dataSettings['iniDischarge'], partitionShape)
        except:
            self.iniDischarge           = dG.generate.lue_zero()
        
        # Initial InterceptionStorage and groundWaterStorage
        try:
            self.iniInterceptionStorage = lfr.from_gdal(self.inputDir + configuration.dataSettings['iniInterceptionStorage'], partitionShape)
        except:
            self.iniInterceptionStorage = dG.generate.lue_zero()

        self.iniGroundWaterStorage  = (self.iniGroundWaterHeight - (self.impermeableLayerHeight)) * self.cellArea
        

    @lfr.runtime_scope
    def dynamicModel(self, configuration):
        dt = int(configuration.modelSettings['iterationsBeforeReport'])
        sD = list(map(int, configuration.modelSettings['startDate'].split(", ")))
        eD = list(map(int, configuration.modelSettings['endDate'].split(", ")))
        startDate   = datetime.datetime(sD[0], sD[1], sD[2], sD[3], sD[4], sD[5])
        endDate     = datetime.datetime(eD[0], eD[1], eD[2], eD[3], eD[4], eD[5])
        dT = int((endDate - startDate).seconds / dt)
        
        # Loading initial conditions
        groundWaterHeight   = self.iniGroundWaterHeight
        discharge           = self.iniDischarge
        height              = height = lfr.pow(coefficient*discharge, 0.6)
        Sgw                 = self.iniGroundWaterStorage
        interceptionStorage = self.iniInterceptionStorage
        
        # Setting min and max soil values
        MaxSgw              =self.impermeableLayerBelowDEM * self.cellArea       # Full storage of porosity
        MinSgw              = MaxSgw * (self.wiltingPoint / self.porosity)          # Minimum storage because of wilting point
        
        # Values for discharge to height calculation
        sqrtSlope           = lfr.sqrt(self.slope)
        sqrtSlope           = lfr.where(sqrtSlope < 0.001, 0.001, sqrtSlope)
        sqrtSlope           = lfr.where(sqrtSlope > 0.05, 0.05, sqrtSlope)
        width               = 1
        coefficient         = self.mannings / (sqrtSlope * width)
        
        # Channel length and area
        channelLength       = self.resolution * dG.generate.lue_one()
        channelArea         = width * channelLength
        
        # Kinematic Surface Water Routing Constants
        alpha       = 1.5
        beta        = 0.6
        timestep    = 1.0 * float(configuration.modelSettings['timestep'])

        # Static, really small value because inflow = 0 is not accepted
        inflow = dG.generate.lue_one()*0.000000000001
        
        # Open file to write maximum discharge values to for post simulation validation.
        with open(self.outputDir + "maximumDischarge.csv", "w", newline="") as f:
            writer = csv.writer(f, delimiter=';')
            
            # Start model for dT large periods
            for i in range(dT):
                # Time in minutes is the small iteration multiplied with the timestep (both in seconds) divided by 60 seconds.
                date = startDate + datetime.timedelta(seconds = i * (dt*timestep)) 
                
                # Load flux and storage values
                precipitation               = dA.get.csvData(date, self.cellArea,
                                                                   configuration.generalSettings['inputDir'] + configuration.dataSettings['precipitationData']) # m/s
                ref_evaporation             = dA.get.csvData(date, self.cellArea, 
                                                                     configuration.generalSettings['inputDir'] + configuration.dataSettings['evapotranspirationData']) # m/s
                interceptionStorage, precipitation, evapotranspirationSurface = \
                                              dA.get.interception(self.cellArea, interceptionStorage, self.interceptionStorageMax, precipitation, \
                                                                  ref_evaporation, self.throughfallFraction)
                evapotranspirationSurface, evapotranspirationSoil = \
                                              dA.get.evapotranspiration(precipitation, evapotranspirationSurface, height)
                infiltration                = dA.get.pot_infiltration(Sgw, MaxSgw, self.cellArea, self.Ks, self.permeability, self.porosity, \
                                                                      height, precipitation, evapotranspirationSurface)


                # Groundwater LDD, gradient and flow flux
                gwLDD       = lfr.d8_flow_direction(groundWaterHeight)
                dHgw        = groundWaterHeight - lfr.downstream(gwLDD, groundWaterHeight)
                gwGradient  = (dHgw) / self.resolution
                dHgw        = lfr.where(dHgw < Sgw, dHgw, Sgw)
                Qgw         = self.Ks * gwGradient * timestep * self.resolution * dHgw * self.resolution                       # Groundwater velocity in m/s
                
                # If the groundwater flow because of the impermeable layer is larger than the amount of water available, than it should be set so only the stored water will move.
                Qgw         = lfr.where(Qgw * dt > Sgw - MinSgw, (Sgw - MinSgw)/dt, Qgw)
                Qgw         = lfr.where(Sgw < MinSgw, 0.000000000001, Qgw)
                
                # Add all vertical processes for the surfacewater and all processes groundwater
                gwFlux      = ((infiltration - evapotranspirationSoil)/self.porosity) + lfr.upstream(gwLDD, Qgw) - Qgw          # Is now in cubic meters
                swFlux      =  precipitation - evapotranspirationSurface - infiltration                                         # Is now in cubic meters
                
                for j in range(dt):
                    # The groundwater is adjusted by the fluxes
                    Sgw         = Sgw + gwFlux                                                                     
                    
                    # If the groundwater table surpases the digital elevation map, groundwater is turned into runoff.
                    seepage     = lfr.where(Sgw > MaxSgw, (Sgw - MaxSgw)*self.porosity, 0)
                    
                    # Discharge is affected by the surfacewater fluxes, and seepage is added
                    height   = height + (swFlux + seepage)/channelArea
                    
                    discharge = lfr.pow(height, 1.6666666667) / coefficient
                    
                    # Because the kinematic wave has difficulties working with zero's, we have opted for a very small value. This will impact model results.
                    discharge   = lfr.where(discharge < 0.000000000001, 0.000000000001, discharge)

                    # Water routing based on the kinematic wave function, currently alpha is a float. Hopefully mannings raster can be used in the future.
                    discharge           = lfr.kinematic_wave(self.ldd, discharge, inflow,\
                                                alpha, beta, timestep,\
                                                channelLength,)
                    
                    height = lfr.pow(coefficient*discharge, 0.6)
                    
                    # Any water that is moved from groundwater to discharge has to be removed from the groundwaterStorage
                    Sgw         = Sgw - (seepage/self.porosity)
                    
                    # Get the maximum value of the discharge raster (to limit the amount of tasks created by HPX)
                    Qmaxvalue = lfr.maximum(discharge).get()
                    print(Qmaxvalue)
                    
                    # Write value to csv for later validation
                    writer.writerow([i*60 + j, Qmaxvalue])
                
                # Adjust the GW Table for the LDD creation of the next timestep.
                groundWaterHeight = self.impermeableLayerHeight + Sgw/self.cellArea
                
                # Save / Report data
                print(f"Done: {i+1}/{dT}")
                variables = {"discharge": discharge, "seepage": seepage, "Qgw": Qgw, "groundWaterHeight": groundWaterHeight, "Sgw": Sgw, "evapoSoil": evapotranspirationSoil, "infiltration": infiltration,
                             "gwFlux": gwFlux, "swFlux": swFlux}

                report.dynamic(date, timestep, variables, self.outputDir)
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
    configuration = Configuration("F:/Projecten intern (2023)/Stage Steven Hosper/Model/v1/config.ini")
    main = mainModel(configuration)
    main.dynamicModel(configuration)

    # Process the results into a gif
    if configuration.generalSettings['makeGIF'] == 'True':
        print(f"Creating a GIF for: {configuration.gifSettings['variables']}.")
        MakeGIF.run(configuration)

print("--- %s seconds ---" % (time.time() - start_time))
