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
    def __init__(self, configuration, report):
        print("Initializing the program...")
        # Initialize submodules
        self.standard_LUE   = StandardArraysLUE(configuration)
        self.retrieve_data  = RetrieveData(configuration)
        self.calculate_flux = CalculateFlux(configuration)
        self.report         = report
        
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
        
        
        #### SETTING CONSTANTS ####
        ## GENERAL ##
        self.resolution                 = float(configuration.modelSettings['resolution'])
        self.cell_area                  = self.resolution * self.resolution
        self.stdmin                     = float(configuration.modelSettings['standard_minimal_value'])
        self.slope          	        = lfr.slope(self.dem, self.resolution)
        self.slope_sqrd                 = utilityFunctions.calculate_sqrd_slope(self.slope, 0.05, 0.00001)
        self.channel_width              = int(configuration.modelSettings['channel_width'])
        self.coefficient                = self.mannings / (0.008 * self.channel_width)
        self.imperm_below_dem           = float(configuration.modelSettings['impermeableLayerBelowDEM'])
        self.imperm_lay_height          = self.dem - self.imperm_below_dem
        self.water_below_dem            = float(configuration.modelSettings['waterBelowDEM'])
        self.max_gw_s                   = self.imperm_below_dem * self.cell_area                        # Full storage of porosity
        self.min_gw_s                   = self.max_gw_s * (self.wilting_point / self.porosity)          # Minimum storage because of wilting point
        # Refactorings value from mm/hour to m/h times the cell area.
        self.refactor                   = (self.cell_area / 1000) / 3600 
        
        # Channel length and area
        self.channel_length      = self.resolution * self.standard_LUE.one()
        self.channel_area        = self.channel_width * self.channel_length
        self.channel_rat         = self.channel_area / self.cell_area
        self.infil_to_gw_s       = self.channel_area / self.porosity
        # self.notBoundaryCells       = generate.boundaryCell() # Currently not working
        
        # Kinematic Surface Water Routing Constants
        self.alpha       = 1.5
        self.beta        = 0.6
        self.timestep    = 1.0 * float(configuration.modelSettings['timestep'])
        self.c           = 5/3

        # Static, really small value because inflow = 0 is not accepted
        self.inflow = self.standard_LUE.one()*self.stdmin
        
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
            self.ini_sur_h           = lfr.from_gdal(self.input_dir + configuration.dataSettings['iniWaterHeight'], partition_shape)
        except:
            print("Did not find a initial discharge file, looked at: {}".format(configuration.dataSettings['iniWaterHeight']))
            self.ini_sur_h           = lfr.where(self.dem > 0, self.standard_LUE.zero(), self.standard_LUE.zero())
        
        # Initial InterceptionStorage and groundWaterStorage
        try:
            self.ini_int_s = lfr.from_gdal(self.input_dir + configuration.dataSettings['iniInterceptionStorage'], partition_shape)
        except:
            print("Did not find a initial interception storage file, looked at: {}".format(configuration.dataSettings['iniInterceptionStorage']))
            self.ini_int_s = lfr.where(self.dem > 0, self.standard_LUE.zero(), self.standard_LUE.zero())

        # REPORTING INITIAL CONDITIONS
        variables = {"ini_gw_s": self.ini_gw_s, "ini_int_s": self.ini_int_s, "ini_sur_h": self.ini_sur_h,
                             }
        self.report.initial(variables)
        
        print("\n")
        
    def update_and_route_both(self, gw_s, gw_flux, height, sw_flux):
        # The groundwater is adjusted by the fluxes
        #channel_infiltation = lfr.where(height > pot_channel_infiltation, pot_channel_infiltation, height)
        gw_s         = gw_s + gw_flux                                #+ channel_infiltation*infil_to_gw_s                                                                 
        
        # If the groundwater table surpases the digital elevation map, groundwater is turned into runoff.
        seepage     = lfr.where(gw_s > self.max_gw_s, (gw_s - self.max_gw_s)*self.porosity, 0)
        
        # Discharge is affected by the surfacewater fluxes, and seepage is added
        height   = height + ((sw_flux + seepage)/self.channel_area)            #- channel_infiltation
        
        discharge = lfr.pow(height, self.c) / self.coefficient
        
        # Because the kinematic wave has difficulties working with zero's, we have opted for a very small value. This will impact model results.
        discharge   = lfr.where(discharge < self.stdmin, self.stdmin, discharge)

        # Water routing based on the kinematic wave function, currently alpha is a float. Hopefully mannings raster can be used in the future.
        discharge           = lfr.kinematic_wave(self.ldd, discharge, self.inflow,\
                                    self.alpha, self.beta, self.timestep,\
                                    self.channel_length)
        
        height = lfr.pow(self.coefficient*discharge, 0.6)
        
        # Any water that is moved from groundwater to discharge has to be removed from the groundwaterStorage
        gw_s         = gw_s - (seepage / self.porosity)
        
        return height, discharge, gw_s

    def update_and_route_surface(self, height, sw_flux):
        # Discharge is affected by the surfacewater fluxes, and seepage is added
        height   = height + ((sw_flux)/self.channel_area)            #- channel_infiltation
        
        discharge = lfr.pow(height, self.c) / self.coefficient
        
        # Because the kinematic wave has difficulties working with zero's, we have opted for a very small value. This will impact model results.
        discharge   = lfr.where(discharge < self.stdmin, self.stdmin, discharge)

        # Water routing based on the kinematic wave function, currently alpha is a float. Hopefully mannings raster can be used in the future.
        discharge           = lfr.kinematic_wave(self.ldd, discharge, self.inflow,\
                                    self.alpha, self.beta, self.timestep,\
                                    self.channel_length)
        
        height = lfr.pow(self.coefficient*discharge, 0.6)
        
        # Any water that is moved from groundwater to discharge has to be removed from the groundwaterStorage
        return height, discharge
    
    def update_and_route_subsurface(self, gw_s, gw_flux):
        # The groundwater is adjusted by the fluxes
        gw_s         = gw_s + gw_flux                                                             
        
        # If the groundwater table surpases the digital elevation map, groundwater is turned into runoff.
        seepage     = lfr.where(gw_s > self.max_gw_s, (gw_s - self.max_gw_s)*self.porosity, 0)
        
        gw_s        = gw_s - (seepage / self.porosity)
        
        # Any water that is moved from groundwater to discharge has to be removed from the groundwaterStorage
        
        return gw_s, seepage
    
    @lfr.runtime_scope
    def dynamic_model(self, configuration):
        dt = int(configuration.modelSettings['iterationsBeforeReport'])
        start_date   = utilityFunctions.string_to_datetime(configuration.modelSettings['startDate'], ", ")
        end_date     = utilityFunctions.string_to_datetime(configuration.modelSettings['endDate'], ", ")
        dT = int((end_date - start_date).total_seconds() / dt)
        
        # Loading initial conditions
        height      = self.ini_sur_h
        gw_s        = self.ini_gw_s
        int_s       = self.ini_int_s
        gw_height   = self.imperm_lay_height + gw_s/self.cell_area        
        discharge   = lfr.pow(height, self.c) / self.coefficient
        
        # Open file to write maximum discharge values to for post simulation validation.
        with open(self.output_dir + "/maximumDischarge.csv", "w", newline="") as f:
            writer = csv.writer(f, delimiter=';')
            
            # Start model for dT large periods
            for i in range(dT):
                # Time in minutes is the small iteration multiplied with the timestep (both in seconds) divided by 60 seconds.
                date = start_date + datetime.timedelta(seconds = i * (dt*self.timestep)) 
                
                # Load flux and storage values
                precipitation = self.retrieve_data.csv_timeseries_to_flux(configuration.generalSettings['inputDir'] +
                                                                          configuration.dataSettings['precipitationData'],
                                                                          self.refactor, date) # m/s
                
                ref_evaporation = self.retrieve_data.csv_timeseries_to_flux(configuration.generalSettings['inputDir'] +
                                                                            configuration.dataSettings['evapotranspirationData'],
                                                                            self.refactor, date) # m/s
                
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
                pot_channel_infiltation = pot_channel_infiltation * self.channel_rat  # is in m/s

                # Groundwater LDD, gradient and flow flux
                gw_ldd          = lfr.d8_flow_direction(gw_height)
                del_h_gw        = gw_height - lfr.downstream(gw_ldd, gw_height)
                gw_grad         = (del_h_gw) / self.resolution
                gw_flow         = self.Ks * gw_grad * self.timestep * (gw_height - self.imperm_lay_height) * self.resolution                       # Groundwater velocity in m2/s
                
                # If the groundwater flow because of the impermeable layer is larger than the amount of water available, than it should be set so only the stored water will move.
                gw_flow         = lfr.where(gw_flow * dt > gw_s - self.min_gw_s, (gw_s - self.min_gw_s)/dt, gw_flow)
                gw_flow         = lfr.where(gw_s < self.min_gw_s, self.stdmin, gw_flow)
                gw_flow         = lfr.where(gw_ldd == 5, self.standard_LUE.zero(), 0.000001)
                
                out_flow = gw_flow
                in_flow  = lfr.upstream(gw_ldd, gw_flow)
                # Add all vertical processes for the surfacewater and all processes groundwater
                gw_flux      =  in_flow - out_flow
                sw_flux      =  precipitation - evapotranspiration_surface - direct_infiltration                                         # Is now in cubic meters
                
                
                ############## ROUTING AND UPDATE PER TIMESTEP ##############
                
                if configuration.modelSettings['routing'] == 'both':            
                    for j in range(dt):
                        # Use the gw_flux and sw_flux to adjust the surface- and groundwater height accordingly
                        # Then route the water lateraly and repeat this for the amount of iterations required.
                        height, discharge, gw_s = self.update_and_route_both(gw_s, gw_flux, height, sw_flux)
                        
                        # Get the maximum value of the discharge raster (to limit the amount of tasks created by HPX)
                        outflow = lfr.minimum(lfr.zonal_sum(discharge, self.ldd == 5)).get()
                        print("outflow: ", outflow)
                        
                        # Write value to csv for later validation
                        writer.writerow([i*60 + j, outflow])
                        
                elif configuration.modelSettings['routing'] == 'surface':
                    for j in range(dt):
                        # Use the gw_flux and sw_flux to adjust the surface- and groundwater height accordingly
                        # Then route the water lateraly and repeat this for the amount of iterations required.
                        height, discharge = self.update_and_route_surface(height, sw_flux)
                        
                        # Get the maximum value of the discharge raster (to limit the amount of tasks created by HPX)
                        outflow = lfr.minimum(lfr.zonal_sum(discharge, self.ldd == 5)).get()
                        total_volume = lfr.minimum(lfr.zonal_sum(height, self.dem>0)).get()
                        print("outflow: ", outflow, "total: ", total_volume)
                        
                        # Write value to csv for later validation
                        writer.writerow([i*60 + j, outflow])
                    
                elif configuration.modelSettings['routing'] == 'subsurface':
                    for j in range(dt):
                        # Use the gw_flux and sw_flux to adjust the surface- and groundwater height accordingly
                        # Then route the water lateraly and repeat this for the amount of iterations required.
                        gw_s, seepage = self.update_and_route_subsurface(gw_s, gw_flux)
                        
                        outflow = lfr.maximum(lfr.zonal_sum(seepage, self.dem>0)).get()
                        total_volume = lfr.minimum(lfr.zonal_sum(gw_s, self.dem>0)).get()
                        print("outflow: ", outflow, "total_volume: ", total_volume)
                        
                        # Write value to csv for later validation
                        writer.writerow([i*60 + j, outflow])
                
                else:
                    raise("Unknown 'routing' mode, change configuration setting to one of the following: both, surface, subsurface")
                    break
                
                
                
                # Adjust the GW Table for the LDD creation of the next timestep.
                gw_height = self.imperm_lay_height + gw_s/self.cell_area
                
                
                
                
                # Save / Report data
                print(f"Done: {i+1}/{dT}")
                variables = {"discharge": discharge, "int_s": int_s, "height": height, "gw_s": gw_s, "gw_flow": gw_flow, "gw_flux": gw_flux
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
    main = mainModel(configuration, report)
    main.dynamic_model(configuration)
    report.balance_report(configuration)  
    
    # Process the results into a gif
    if configuration.generalSettings['makeGIF'] == 'True':
        print(f"Creating a GIF for: {configuration.gifSettings['variables']}.")
        tools.MakeGIF.run(configuration)

print("--- %s seconds ---" % (time.time() - start_time))
