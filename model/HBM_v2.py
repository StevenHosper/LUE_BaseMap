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
#import tools.MakeGIF

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


class MyModel(lfr.Model):
    def __init__(self, configuration, report):
        lfr.Model.__init__(self)
        self.configuration  = configuration
        self.standard_LUE   = StandardArraysLUE(configuration)
        self.retrieve_data  = RetrieveData(configuration)
        self.calculate_flux = CalculateFlux(configuration)
        self.report         = report
        
    def initialize(self):
        # Set directories
        self.input_dir   = self.configuration.generalSettings['inputDir'] + self.configuration.generalSettings['scenario'] 
        
        # Timesteps (all should have same unit, seconds in this case)
        self.timestep        = float(self.configuration.modelSettings['timestep'])
        self.update_timestep = float(self.configuration.modelSettings['iterationsBeforeData'])
        self.report_timestep = float(self.configuration.modelSettings['iterationsBeforeReport'])
        
        partition_shape  = 2 * (self.configuration.modelSettings['partitionExtent'],)
        
        # Initialize data required from memory files
        # Get all constants        
        self.dem         = lfr.from_gdal(self.input_dir + self.configuration.dataSettings['dem'], partition_shape)              # DEM map of the study area
        land_use         = lfr.from_gdal(self.input_dir + self.configuration.dataSettings['landUseMap'], partition_shape)       # Land-use, example: road
        soil_type        = lfr.from_gdal(self.input_dir + self.configuration.dataSettings['soilMap'], partition_shape)          # example: sand or clay
        self.resolution  = float(self.configuration.modelSettings['resolution'])
        self.cell_area   = self.resolution * self.resolution
        
        # Retrieve the soil properties
        self.Ks, self.porosity, self.wilting_point, self.Ki = self.retrieve_data.soil_csv(
            self.configuration.generalSettings['inputDir'] + self.configuration.dataSettings['soilData'], soil_type)               # soil characteristic
        
        # Retrieve the land-use properties
        self.mannings, self.permeability, self.max_int_s, self.throughfall_frac = \
            self.retrieve_data.land_characteristics_csv(
                self.configuration.generalSettings['inputDir'] + self.configuration.dataSettings['landUseData'], land_use,
                self.cell_area)         # land-use characteristics
        
        self.gw_base = float(self.configuration.modelSettings['groundWaterBase']) * self.standard_LUE.one()
        # self.ldd = lfr.d8_flow_direction(self.dem)
        self.ldd        = lfr.from_gdal(self.input_dir + self.configuration.dataSettings['ldd'], partition_shape)
        
        
        #### SETTING CONSTANTS ####
        ## GENERAL ##
        self.routing                    = self.configuration.modelSettings['routing']
        self.start_date                 = utilityFunctions.string_to_datetime(self.configuration.modelSettings['startDate'], ", ")
        self.stdmin                     = float(self.configuration.modelSettings['standard_minimal_value'])
        slope          	                = self.standard_LUE.one() * 0.008
        self.slope_sqrd                 = utilityFunctions.calculate_sqrd_slope(slope, 0.05, 0.00001)
        self.channel_width              = int(self.configuration.modelSettings['channel_width'])
        self.coefficient                = self.mannings / (0.008 * self.channel_width)
        aquifer_height                  = float(self.configuration.modelSettings['impermeableLayerBelowDEM'])
        self.imperm_lay_height          = self.dem - aquifer_height
        self.water_below_dem            = float(self.configuration.modelSettings['waterBelowDEM'])
        self.max_gw_s                   = aquifer_height * self.cell_area                             # Full storage of porosity
        self.min_gw_s                   = self.max_gw_s * (self.wilting_point / self.porosity)          # Minimum storage because of wilting point
        
        # Refactorings value from mm/hour to m/s times the cell area.
        self.refactor                   = (self.cell_area / 1000) / 3600 
        self.d                          = (5**(2/3))/(4**(2/3))
        
        # Channel length and area
        self.channel_length      = self.resolution * self.standard_LUE.one()
        self.channel_area        = self.channel_width * self.channel_length
        
        # Kinematic Surface Water Routing Constants
        self.alpha       = 1.5
        self.beta        = 0.6
        self.c           = 5/3
        
        # Load initial groundWaterStorage, if no raster is supplied, use the waterBelowDEM in combination with DEM to create a initialGroundWaterStorage layer.
        try:
            self.gw_s   = lfr.from_gdal(self.input_dir + self.configuration.dataSettings['iniGroundWaterStorage'], partition_shape)
        except:
            print("Did not find a initial groundwater height file, looked at: {}".format(self.configuration.dataSettings['iniGroundWaterStorage']))
            self.gw_s   = lfr.where(self.dem > (self.gw_base + self.water_below_dem),
                                        ((self.dem - self.water_below_dem)-self.imperm_lay_height) * self.cell_area,
                                        (self.gw_base - self.imperm_lay_height) * self.cell_area)
            self.gw_s   = lfr.where(self.gw_s > self.max_gw_s, self.max_gw_s, self.gw_s)
        
        # Load initial discharge, if no raster is supplied, set to zero.
        try:
            self.height              = lfr.from_gdal(self.input_dir + self.configuration.dataSettings['iniWaterHeight'], partition_shape)
        except:
            print("Did not find a initial discharge file, looked at: {}".format(self.configuration.dataSettings['iniWaterHeight']))
            self.height              = lfr.where(self.dem > 0, self.standard_LUE.zero(), self.standard_LUE.zero())
            
        try:
            self.discharge           = lfr.from_gdal(self.input_dir + self.configuration.dataSettings['iniDischarge'], partition_shape)
        except:
            print("Did not find a initial discharge file, looked at: {}".format(self.configuration.dataSettings['iniDischarge']))
            self.discharge           = lfr.where(self.dem > 0, self.standard_LUE.zero(), self.standard_LUE.zero())
        
        # Initial InterceptionStorage and groundWaterStorage
        try:
            self.int_s = lfr.from_gdal(self.input_dir + self.configuration.dataSettings['iniInterceptionStorage'], partition_shape)
        except:
            print("Did not find a initial interception storage file, looked at: {}".format(self.configuration.dataSettings['iniInterceptionStorage']))
            self.int_s = lfr.where(self.dem > 0, self.standard_LUE.zero(), self.standard_LUE.zero())

        # REPORTING INITIAL CONDITIONS
        variables = {"ini_gw_s": self.gw_s, "ini_int_s": self.int_s, "ini_sur_h": self.height,
                             }
        
        self.report.initial(variables)

    
    def update_all_fluxes(self, gw_s, date, update_timestep):
        # Load flux and storage values
        precipitation = self.retrieve_data.csv_timeseries_to_flux(self.configuration.generalSettings['inputDir'] +
                                                                  self.configuration.dataSettings['precipitationData'],
                                                                  self.refactor, date) # m3/s into a cell (refactor incorporates cell area)
                
        ref_evaporation = self.retrieve_data.csv_timeseries_to_flux(self.configuration.generalSettings['inputDir'] +
                                                                    self.configuration.dataSettings['evapotranspirationData'],
                                                                    self.refactor, date) # m3/s into a cell (refactor incorporates cell area)
        
        self.int_s, precipitation, evapotranspiration_surface = self.calculate_flux.interception(self.int_s,
                                                                                            self.max_int_s,
                                                                                            precipitation,
                                                                                            ref_evaporation,
                                                                                            self.throughfall_frac)
        
        evapotranspiration_surface, evapotranspiration_soil = self.calculate_flux.evapotranspiration(precipitation,
                                                                                                        evapotranspiration_surface)
        
        """ direct_infiltration = self.calculate_flux.infiltration(gw_s,
                                                                self.max_gw_s,
                                                                self.Ks,
                                                                self.permeability,
                                                                self.porosity,
                                                                precipitation,
                                                                evapotranspiration_surface) """
        
        direct_infiltration, pot_reinfiltration = self.calculate_flux.adjusted_Horton(gw_s,
                                                                self.max_gw_s,
                                                                self.Ks,
                                                                self.Ki,
                                                                self.permeability,
                                                                self.porosity,
                                                                precipitation,
                                                                evapotranspiration_surface)
        
        # Groundwater LDD, gradient and flow flux
        gw_height       = self.imperm_lay_height + self.gw_s / self.cell_area
        gw_ldd          = lfr.d8_flow_direction(gw_height)
        del_h_gw        = gw_height - lfr.downstream(gw_ldd, gw_height)
        gw_grad         = (del_h_gw) / self.resolution
        gw_flow         = self.Ks * gw_grad * self.timestep * (gw_height - self.imperm_lay_height) * self.resolution                       # Groundwater velocity in m2/s
        
        # If the groundwater flow because of the impermeable layer is larger than the amount of water available, than it should be set so only the stored water will move.
        gw_flow         = lfr.where(gw_flow * update_timestep > gw_s - self.min_gw_s, (gw_s - self.min_gw_s)/update_timestep, gw_flow)
        gw_flow         = lfr.where(gw_s < self.min_gw_s, self.stdmin, gw_flow)
        gw_flow         = lfr.where(gw_ldd == 5, self.standard_LUE.zero(), gw_flow)

        # Add all vertical processes for the surfacewater and all processes groundwater
        gw_flux      =  direct_infiltration - evapotranspiration_soil + lfr.upstream(gw_ldd, gw_flow) - gw_flow
        sw_flux      =  precipitation - evapotranspiration_surface - direct_infiltration                                         # Is now in cubic meters
        
        return gw_flux, sw_flux, pot_reinfiltration

    def update_surface_fluxes(self, date):
        # Load flux and storage values
        gw_s = self.standard_LUE.one() * 1
        
        precipitation = self.retrieve_data.csv_timeseries_to_flux(self.configuration.generalSettings['inputDir'] +
                                                                  self.configuration.dataSettings['precipitationData'],
                                                                  self.refactor, date) # m3/s into a cell (refactor incorporates cell area)
                
        ref_evaporation = self.retrieve_data.csv_timeseries_to_flux(self.configuration.generalSettings['inputDir'] +
                                                                    self.configuration.dataSettings['evapotranspirationData'],
                                                                    self.refactor, date) # m3/s into a cell (refactor incorporates cell area)
        
        self.int_s, precipitation, evapotranspiration_surface = self.calculate_flux.interception(self.int_s,
                                                                                            self.max_int_s,
                                                                                            precipitation,
                                                                                            ref_evaporation,
                                                                                            self.throughfall_frac)
        
        evapotranspiration_surface, evapotranspiration_soil = self.calculate_flux.evapotranspiration(precipitation,
                                                                                                        evapotranspiration_surface)
        
        direct_infiltration = self.calculate_flux.infiltration(gw_s,
                                                                self.max_gw_s,
                                                                self.Ks,
                                                                self.permeability,
                                                                self.porosity,
                                                                precipitation,
                                                                evapotranspiration_surface)
        
        sw_flux      =  precipitation - evapotranspiration_surface - direct_infiltration                                         # Is now in cubic meters
        
        return sw_flux

    def update_subsurface_fluxes(self, gw_s, date, update_timestep):
        # Load flux and storage values
        precipitation = self.retrieve_data.csv_timeseries_to_flux(self.configuration.generalSettings['inputDir'] +
                                                                  self.configuration.dataSettings['precipitationData'],
                                                                  self.refactor, date) # m3/s into a cell (refactor incorporates cell area)
                
        ref_evaporation = self.retrieve_data.csv_timeseries_to_flux(self.configuration.generalSettings['inputDir'] +
                                                                    self.configuration.dataSettings['evapotranspirationData'],
                                                                    self.refactor, date) # m3/s into a cell (refactor incorporates cell area)
        
        self.int_s, precipitation, evapotranspiration_surface = self.calculate_flux.interception(self.int_s,
                                                                                            self.max_int_s,
                                                                                            precipitation,
                                                                                            ref_evaporation,
                                                                                            self.throughfall_frac)
        
        evapotranspiration_surface, evapotranspiration_soil = self.calculate_flux.evapotranspiration(precipitation,
                                                                                                        evapotranspiration_surface)
        
        direct_infiltration = self.calculate_flux.infiltration(gw_s,
                                                                self.max_gw_s,
                                                                self.Ks,
                                                                self.permeability,
                                                                self.porosity,
                                                                precipitation,
                                                                evapotranspiration_surface)
        
        # Groundwater LDD, gradient and flow flux
        gw_height       = self.imperm_lay_height + self.gw_s / self.cell_area
        gw_ldd          = lfr.d8_flow_direction(gw_height)
        del_h_gw        = gw_height - lfr.downstream(gw_ldd, gw_height)
        gw_grad         = (del_h_gw) / self.resolution
        gw_flow         = self.Ks * gw_grad * self.timestep * (gw_height - self.imperm_lay_height) * self.resolution                       # Groundwater velocity in m2/s
        
        # If the groundwater flow because of the impermeable layer is larger than the amount of water available, than it should be set so only the stored water will move.
        gw_flow         = lfr.where(gw_flow * update_timestep > gw_s - self.min_gw_s, (gw_s - self.min_gw_s)/update_timestep, gw_flow)
        gw_flow         = lfr.where(gw_s < self.min_gw_s, self.stdmin, gw_flow)
        gw_flow         = lfr.where(gw_ldd == 5, self.standard_LUE.zero(), gw_flow)

        # Add all vertical processes for the surfacewater and all processes groundwater
        gw_flux      =  direct_infiltration - evapotranspiration_soil + lfr.upstream(gw_ldd, gw_flow) - gw_flow
        
        return gw_flux

    def route_both(self, gw_s, gw_flux, discharge, sw_flux, pot_reinfiltration):
        """Route surface and subsurface.
        
        We use the available information in regards to surface and subsurface
        layers in order to determine the values of each variable in the next
        state. Both layers are updated and routed.
        
        Args:
            gw_s (LUE partitioned array): the groundwater storage in m3 of soil filled with water (m3).
            gw_flux (LUE partitioned array): the groundwater flux based on input and flow (m3/s).
            discharge (LUE partitioned array): the flow volume that travels through each cell (m3/s).
            sw_flux (LUE partitioned array): the surface water flux (m3/s).
            pot_reinfiltration (LUE partitioned array): the potential reinfiltration of runoff (m3/s).

        Returns:
            height (LUE partitioned array): the water height in each cell (m).
            updated_discharge (LUE partitioned array): the updated flow volume that travels through each cell (m3/s).
            updated_gw_s (LUE partitioned array): the updated groundwater storage in m3 of soil filled with water. (m3)
            seepage (LUE partitioned array): the amount of water leaving the soil (m3/s).
            reinfiltration (LUE partitioned array): the amount of runoff that infiltrates the soil (m3/s).
        """
        # The groundwater is adjusted by the fluxes
        height = lfr.pow(self.coefficient*discharge, 0.6)
        reinfiltration = lfr.where(height > pot_reinfiltration, pot_reinfiltration, height)
        gw_s         = gw_s + gw_flux + (reinfiltration / self.porosity)                                                
        
        # If the groundwater table surpases the digital elevation map, groundwater is turned into runoff.
        seepage     = lfr.where(gw_s > self.max_gw_s, (gw_s - self.max_gw_s)*self.porosity, 0)

        height    = height - reinfiltration
        discharge = lfr.pow(height, self.c) / self.coefficient
        
        # Because the kinematic wave has difficulties working with zero's, we have opted for a very small value. This will impact model results.
        discharge   = lfr.where(discharge < self.stdmin, self.stdmin, discharge)
        inflow      = (sw_flux + seepage)/self.channel_area
        inflow      = lfr.where(inflow < self.stdmin, self.stdmin, inflow)

        # Water routing based on the kinematic wave function, currently alpha is a float. Hopefully mannings raster can be used in the future.
        discharge           = lfr.kinematic_wave(self.ldd, discharge, inflow,\
                                    self.alpha, self.beta, self.timestep,\
                                    self.channel_length)
        
        height = lfr.pow(self.coefficient*discharge, 0.6)
        
        # Any water that is moved from groundwater to discharge has to be removed from the groundwaterStorage
        gw_s        = gw_s - (seepage / self.porosity)
        
        return height, discharge, gw_s, seepage, reinfiltration
    
    def route_surface(self, height, sw_flux):
        """Route only the surface

        The height and surface water flux are used to calculate a new 
        water table height, and route water in downstream direction.

        Args:
            height (LUE partitioned array): the old water height (m).
            sw_flux (LUE partitioned array): the surface water flux (m3/s).

        Returns:
            updated_height (LUE partitioned array): the water height in each cell (m).
            discharge (LUE partitioned array): the flow volume that travels through each cell (m3/s).
        """
        # Discharge is affected by the surfacewater fluxes, and seepage is added
        new_height   = height + ((sw_flux)/self.channel_area)            #- channel_infiltation
        
        discharge = lfr.pow(new_height, self.c) / self.coefficient

        # Because the kinematic wave has difficulties working with zero's, we have opted for a very small value. This will impact model results.
        discharge   = lfr.where(discharge < self.stdmin, self.stdmin, discharge)

        # Water routing based on the kinematic wave function, currently alpha is a float. Hopefully mannings raster can be used in the future.
        discharge           = lfr.kinematic_wave(self.ldd, discharge, self.stdmin,\
                                    self.alpha, self.beta, self.timestep,\
                                    self.channel_length)
        
        updated_height = lfr.pow(self.coefficient*discharge, 0.6)

        return updated_height, discharge
    
    def route_subsurface(self, gw_s, gw_flux):
        """Route only the subsurface
        
        The subsurface is routed based on the supplied groundwater storage,
        and the given groundwater flux. The new groundwater storage and the
        seepage are returned by the function.

        Args:
            gw_s (LUE partitioned array): The groundwater storage in m3 of soil filled with water.
            gw_flux (LUE partitioned array): The groundwater flux based on input and flow (m3/s).

        Returns:
            updated_gw_s (LUE partitioned array): The updated groundwater storage in m3 of soil filled with water.
            seepage (LUE partitioned array): The amount of water leaving the soil (m3/s).
        """
        # The groundwater is adjusted by the fluxes
        new_gw_s         = gw_s + gw_flux                                                             
        
        # If the groundwater table surpases the digital elevation map, groundwater is turned into runoff.
        seepage     = lfr.where(new_gw_s > self.max_gw_s, (new_gw_s - self.max_gw_s)*self.porosity, 0)
        
        updated_gw_s        = new_gw_s - (seepage / self.porosity)
        
        return updated_gw_s, seepage
    
    def simulate(self, time_step):
        date = self.start_date + datetime.timedelta(seconds=time_step)
        
        # Update all variables and route.
        if self.routing == "both":
            # Get new fluxes
            if time_step % self.update_timestep:
                self.gw_flux, self.sw_flux, self.pot_reinfiltration = self.update_all_fluxes(self.gw_s, date, self.update_timestep)
            
            
            self.height, self.discharge, self.gw_s, self.seepage, self.reinfiltration = self.route_both(
                                                                                self.gw_s,
                                                                                self.gw_flux,
                                                                                self.discharge,
                                                                                self.sw_flux,
                                                                                self.pot_reinfiltration
                                                                                )
            
            variables = {"discharge": self.discharge, "int_s": self.int_s, 
                         "gw_s": self.gw_s, "seepage": self.seepage, 
                         "reinfiltration": self.reinfiltration
                            }
        
        elif self.routing == 'surface':
            if time_step % self.update_timestep:
                self.gw_flux, self.sw_flux, self.pot_reinfiltration = self.update_surface_fluxes(date)
            
            self.height, self.discharge = self.update_and_route_surface(self.height,
                                                                        self.sw_flux
                                                                        )
            
            variables = {"discharge": self.discharge, "int_s": self.int_s,
                         "height": self.height
                            }
            
            
        elif self.routing == 'subsurface':
            if time_step % self.update_timestep:
                self.gw_flux, self.sw_flux, self.pot_reinfiltration = self.update_subsurface_fluxes(self.gw_s, date, self.update_timestep)
                
            self.gw_s, self.seepage = self.update_and_route_subsurface(self.gw_s, self.gw_flux)
            
            variables = {"gw_s": self.gw_s, "seepage": self.seepage, 
                         "groundwater_height": self.gw_height,
                            }
        
        else:
            raise("Unknown 'routing' mode, change configuration setting to one of the following: both, surface, subsurface")
        
        # Report output
        if time_step % self.report_timestep == 0:
            self.report.dynamic(date, variables)
            
        
    
    def finalize(self):
        # Create the balance report
        self.report.balance_report(self.configuration)  
    
        # Process the results into a gif
        if self.configuration.generalSettings['makeGIF'] == 'True':
            print(f"Creating a GIF for: {self.configuration.gifSettings['variables']}.")
            #tools.MakeGIF.run(self.configuration)

class MyProgressor(lfr.Progressor):
    def __init__(self):
        lfr.Progressor.__init__(self)
        

    def initialize(self):
        sys.stdout.write("[")
        sys.stdout.flush()

    def simulate(self, time_step):
        sys.stdout.write(".")
        sys.stdout.flush()

    def finalize(self):
        sys.stdout.write("]\n")
        sys.stdout.flush()

def calculate_timesteps(start_date: str, end_date: str, sep, timestep_size):
    sd  = utilityFunctions.string_to_datetime(start_date, sep)
    ed  = utilityFunctions.string_to_datetime(end_date, sep)
    timesteps   = int((ed - sd).total_seconds() / int(timestep_size))
    return timesteps

@lfr.runtime_scope
def main():
    if len(sys.argv) == 1:
        sys.exit("Missing configuration file.")
    configuration = Configuration(sys.argv[1])
    report        = Report(configuration)
    model = MyModel(configuration, report)
    progressor = MyProgressor()
    nr_time_steps = calculate_timesteps(configuration.modelSettings['startDate'],
                                        configuration.modelSettings['endDate'],
                                        ", ",
                                        configuration.modelSettings['timestep']
                                        )
    lfr.run_deterministic(model, progressor, nr_time_steps, rate_limit=15)

if __name__ == "__main__":
    main()
    sys.exit(print("--- %s seconds ---" % (time.time() - start_time)))


