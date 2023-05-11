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

# Submodules
from configuration_v2 import Configuration
from reporting import Report
from StandardArraysLUE import StandardArraysLUE
from RetrieveData import RetrieveData
from CalculateFlux import CalculateFlux
from utilityFunctionsHBM import utilityFunctions

# Other tools
import tools.MakeGIF

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

class mainModel:
    def __init__(self, configuration):
        print("Initializing the program...")
        # Initialize submodules
        self.standard_LUE   = StandardArraysLUE(configuration)
        self.retrieve_data  = RetrieveData(configuration)
        self.calculate_flux = CalculateFlux(configuration)
        
        # Set directories
        self.input_dir   = configuration.generalSettings['inputDir'] + configuration.generalSettings['scenario'] 
        self.output_dir  = configuration.generalSettings['outputDir'] + configuration.generalSettings['scenario']
        
        partition_shape  = 2 * (configuration.modelSettings['partitionExtent'],)
        
        # Initialize data required from memory files
        # Get all constants        
        self.dem        = lfr.from_gdal(self.input_dir + configuration.dataSettings['dem'], partition_shape)              # DEM map of the study area
        self.dem        = lfr.where(self.dem < 0.1, 35, self.dem)
        land_use         = lfr.from_gdal(self.input_dir + configuration.dataSettings['landUseMap'], partition_shape)       # Land-use, example: road
        soil_type        = lfr.from_gdal(self.input_dir + configuration.dataSettings['soilMap'], partition_shape)          # example: sand or clay
        
        # Retrieve the soil properties
        self.Ks, self.porosity, self.wilting_point = self.retrieve_data.soil_csv(
            configuration.generalSettings['inputDir'] + configuration.dataSettings['soilData'], soil_type)               # soil characteristic
        
        # Retrieve the land-use properties
        self.mannings, self.permeability, self.max_int_s, self.throughfall_frac = \
            self.retrieve_data.land_characteristics_csv(
                configuration.generalSettings['inputDir'] + configuration.dataSettings['landUseData'], land_use)         # land-use characteristics
        
        self.gw_base = float(configuration.modelSettings['groundWaterBase']) * self.standard_LUE.one()
        # self.ldd = lfr.d8_flow_direction(self.dem)
        self.ldd        = lfr.from_gdal(self.input_dir + configuration.dataSettings['ldd'], partition_shape)
        # Set constants
        self.resolution                 = float(configuration.modelSettings['resolution'])
        self.cell_area                  = self.resolution * self.resolution
        self.slope          	        = lfr.slope(self.dem, self.resolution)
        self.imperm_below_dem           = float(configuration.modelSettings['impermeableLayerBelowDEM'])
        self.imperm_lay_height          = self.dem - self.imperm_below_dem
        self.water_below_dem            = float(configuration.modelSettings['waterBelowDEM'])
        self.max_gw_s                   = self.imperm_below_dem * self.cell_area            # Full storage of porosity
        self.min_gw_s                   = self.max_gw_s * (self.wilting_point / self.porosity)        # Minimum storage because of wilting point
        # self.notBoundaryCells       = generate.boundaryCell() # Currently not working
        
        # Load initial groundWaterStorage, if no raster is supplied, use the waterBelowDEM in combination with DEM to create a initialGroundWaterStorage layer.
        try:
            self.ini_gw_s   = lfr.from_gdal(self.input_dir + configuration.dataSettings['iniGroundWaterStorage'], partition_shape)
        except:
            print("Did not find a initial groundwater height file, looked at: {}".format(configuration.dataSettings['iniGroundWaterStorage']))
            self.ini_gw_s   = lfr.where(self.dem > (self.gw_base + self.water_below_dem),
                                        ((self.dem - self.water_below_dem)-self.imperm_lay_height) * self.cell_area,
                                        (self.gw_base - self.imperm_below_dem) * self.cell_area)
            self.ini_gw_s   = lfr.where(self.ini_gw_s > self.max_gw_s, self.max_gw_s, self.ini_gw_s)
        
        # Load initial discharge, if no raster is supplied, set to zero.
        try:
            self.ini_water_h           = lfr.from_gdal(self.input_dir + configuration.dataSettings['iniWaterHeight'], partition_shape)
        except:
            print("Did not find a initial discharge file, looked at: {}".format(configuration.dataSettings['iniWaterHeight']))
            self.ini_water_h           = self.standard_LUE.zero()
        
        # Initial InterceptionStorage and groundWaterStorage
        try:
            self.ini_int_s = lfr.from_gdal(self.input_dir + configuration.dataSettings['iniInterceptionStorage'], partition_shape)
        except:
            print("Did not find a initial interception storage file, looked at: {}".format(configuration.dataSettings['iniInterceptionStorage']))
            self.ini_int_s = self.standard_LUE.zero()

        print("\n")
        
    def update_and_route(self):
        pass

    @lfr.runtime_scope
    def dynamic_model(self, configuration, report):
        dt = int(configuration.modelSettings['iterationsBeforeReport'])
        start_date   = utilityFunctions.string_to_datetime(configuration.modelSettings['startDate'], ", ")
        end_date     = utilityFunctions.string_to_datetime(configuration.modelSettings['endDate'], ", ")
        dT = int((end_date - start_date).seconds / dt)
        
        # Loading initial conditions
        height      = self.ini_water_h
        gw_s        = self.ini_gw_s
        int_s       = self.ini_int_s
        gw_height   = self.imperm_lay_height + gw_s/self.cell_area
        
        # Values for discharge to height calculation
        slope_sqrd  = utilityFunctions.calculate_sqrd_slope(self.slope, 0.05, 0.00001)
        width       = 1
        coefficient = self.mannings / (slope_sqrd * width)
        
        # Channel length and area
        channel_length      = self.resolution * self.standard_LUE.one()
        channel_area        = width * channel_length
        channel_rat         = channel_area / self.cell_area
        infil_to_gw_s       = channel_area / self.porosity
        
        # Refactorings value from mm/hour to m/h times the cell area.
        refactor            = (self.cell_area / 1000) / 3600         
        
        # Kinematic Surface Water Routing Constants
        alpha       = 1.5
        beta        = 0.6
        timestep    = 1.0 * float(configuration.modelSettings['timestep'])
        c           = 5/3

        # Static, really small value because inflow = 0 is not accepted
        inflow = self.standard_LUE.one()*1E-20
        
        # Open file to write maximum discharge values to for post simulation validation.
        with open(self.output_dir + "/maximumDischarge.csv", "w", newline="") as f:
            writer = csv.writer(f, delimiter=';')
            
            # Start model for dT large periods
            for i in range(dT):
                # Time in minutes is the small iteration multiplied with the timestep (both in seconds) divided by 60 seconds.
                date = start_date + datetime.timedelta(seconds = i * (dt*timestep)) 
                
                # Load flux and storage values
                precipitation = self.retrieve_data.csv_timeseries_to_flux(configuration.generalSettings['inputDir'] +
                                                                          configuration.dataSettings['precipitationData'],
                                                                          refactor, date) # m/s
                
                ref_evaporation = self.retrieve_data.csv_timeseries_to_flux(configuration.generalSettings['inputDir'] +
                                                                            configuration.dataSettings['evapotranspirationData'],
                                                                            refactor, date) # m/s
                
                int_s, precipitation, evapotranspiration_surface = self.calculate_flux.interception(int_s,
                                                                                                    self.max_int_s,
                                                                                                    precipitation,
                                                                                                    ref_evaporation,
                                                                                                    self.throughfall_frac)
                
                evapotranspiration_surface, evapotranspiration_soil = self.calculate_flux.evapotranspiration(precipitation,
                                                                                                             evapotranspiration_surface)
                
                direct_infiltration, pot_channel_infiltation = self.calculate_flux.infiltration(gw_s,
                                                                                                self.max_gw_s,
                                                                                                self.Ks,
                                                                                                self.permeability,
                                                                                                self.porosity,
                                                                                                precipitation,
                                                                                                evapotranspiration_surface)
                
                # The infiltration happens only in the region that is used by the channel and therefore this factor should be accounted for
                pot_channel_infiltation = pot_channel_infiltation * channel_rat  # is in m/s

                # Groundwater LDD, gradient and flow flux
                gw_ldd          = lfr.d8_flow_direction(gw_height)
                del_h_gw        = gw_height - lfr.downstream(gw_ldd, gw_height)
                gw_grad         = (del_h_gw) / self.resolution
                gw_flow         = self.Ks * gw_grad * timestep * (gw_height - self.imperm_lay_height) * self.resolution                       # Groundwater velocity in m2/s
                
                # If the groundwater flow because of the impermeable layer is larger than the amount of water available, than it should be set so only the stored water will move.
                gw_flow         = lfr.where(gw_flow * dt > gw_s - self.min_gw_s, (gw_s - self.min_gw_s)/dt, gw_flow)
                gw_flow         = lfr.where(gw_s < self.min_gw_s, 1E-20, gw_flow)
                
                # Add all vertical processes for the surfacewater and all processes groundwater
                gw_flux      = ((direct_infiltration - evapotranspiration_soil)/self.porosity) + lfr.upstream(gw_ldd, gw_flow) - gw_flow          # Is now in cubic meters
                sw_flux      =  precipitation - evapotranspiration_surface - direct_infiltration                                         # Is now in cubic meters
                
                for j in range(dt):
                    # The groundwater is adjusted by the fluxes
                    # channel_infiltation = lfr.where(height > pot_channel_infiltation, pot_channel_infiltation, height)
                    gw_s         = gw_s + gw_flux                                #+ channel_infiltation*infil_to_gw_s                                                                 
                    
                    # If the groundwater table surpases the digital elevation map, groundwater is turned into runoff.
                    seepage     = lfr.where(gw_s > self.max_gw_s, (gw_s - self.max_gw_s)*self.porosity, 0)
                    
                    # Discharge is affected by the surfacewater fluxes, and seepage is added
                    height   = height + ((sw_flux + seepage)/channel_area)            #- channel_infiltation
                    
                    discharge = lfr.pow(height, c) / coefficient
                    
                    # Because the kinematic wave has difficulties working with zero's, we have opted for a very small value. This will impact model results.
                    discharge   = lfr.where(discharge < 1E-20, 1E-20, discharge)

                    # Water routing based on the kinematic wave function, currently alpha is a float. Hopefully mannings raster can be used in the future.
                    discharge           = lfr.kinematic_wave(self.ldd, discharge, inflow,\
                                                alpha, beta, timestep,\
                                                channel_length,)
                    
                    height = lfr.pow(coefficient*discharge, 0.6)
                    
                    # Any water that is moved from groundwater to discharge has to be removed from the groundwaterStorage
                    gw_s         = gw_s - (seepage / self.porosity)
                    
                    # Get the maximum value of the discharge raster (to limit the amount of tasks created by HPX)
                    outflow = lfr.minimum(lfr.zonal_sum(discharge, self.ldd == 5)).get()
                    print("outflow: ", outflow)
                    
                    # Write value to csv for later validation
                    writer.writerow([i*60 + j, outflow])
                
                # Adjust the GW Table for the LDD creation of the next timestep.
                gw_height = self.imperm_lay_height + gw_s/self.cell_area
                
                # Save / Report data
                print(f"Done: {i+1}/{dT}")
                variables = {"discharge": discharge, "int_s": int_s, "height": height, "gw_s": gw_s,
                             }
                report.dynamic(date, variables)   
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
    configuration = Configuration("F:/Projecten intern (2023)/Stage Steven Hosper/Model/v1/config/config.ini")
    report        = Report(configuration)
    main = mainModel(configuration)
    main.dynamic_model(configuration, report)
    report.balance_report(configuration)  
    
    # Process the results into a gif
    if configuration.generalSettings['makeGIF'] == 'True':
        print(f"Creating a GIF for: {configuration.gifSettings['variables']}.")
        tools.MakeGIF.run(configuration)

print("--- %s seconds ---" % (time.time() - start_time))
