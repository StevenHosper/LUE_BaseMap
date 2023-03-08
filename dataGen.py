# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 10:03:54 2023

@author: steven.hosper
"""

import lue.framework as lfr
import numpy as np
import requests
import datetime
import math as math
import configuration as config

class generate(): 
    def __init__(self):
        """
        The initial calculation phase. Get API access for this session and interpret command line arguments.
        
        """
        
        # Get access to the APIs for this session
        username = '__key__'
        password = 'Cy0BNm8p.vpytC2vYPT9g7OKdgxvqggyV0k9zzJVy'
        
        self.s = requests.Session()
        self.s.headers = {'username': username,
                          'password': password,
                          'Content-Type': 'application/json',
                          }
        
        
        print(f'partition: {config.partitionShape}', f'array: {config.arrayShape}')
        
        self.path = f'{config.path}/data/generated/{config.arrayExtent}/'
        
        print("__init__ done")
    
    def lue_zero():
        return lfr.create_array(config.arrayShape,
                                config.partitionShape,
                                dtype = np.dtype(np.float32),
                                fill_value = 0,
                                )
    
    def lue_one():
        return lfr.create_array(config.arrayShape,
                                config.partitionShape,
                                dtype = np.dtype(np.float32),
                                fill_value = 1,
                            )

    def lue_sink():
        return lfr.create_array(config.arrayShape,
                                config.partitionShape,
                                dtype = np.dtype(np.uint8),
                                fill_value = 5,
                                )
    
    def dem(self, output_name):
        """
        Create a random uniform map, able to function as a Digital Elevation Map and sade it in the \
        generated data directory.
        """
        dem = lfr.uniform(self.lue_zero(), np.float32, 0, 1)
        lfr.to_gdal(dem, f'{self.path}{output_name}_{config.arrayExtent}.tiff')
        return dem
    
    def ldd(self, dem, output_name):
        """
        Create a local drain direction map based on the dem given and safe it to a tiff file in the generated \
        data directory.
        """
        ldd = lfr.d8_flow_direction(dem)
        lfr.to_gdal(ldd, f'{self.path}{output_name}_{config.arrayExtent}.tiff')
        return ldd 
     
    def initialize_raincells(self):
        """
        Assign $% of the cells to become rain producing cells.
        """
        fraction_raincells = 0.05                                           # Determine the percentage of raincells in the array
        raincells = lfr.uniform(self.lue_zero(), np.float32, 0, 1) <= fraction_raincells
        return raincells
    
    def precipitation(self, raincells):
        """
        Simulate precipitation based on kernels and random assigned cells. For now \
        this functions as a simple method to have varied rainclouds that produce rain.
        """
        kernel = np.array(
            [
                [1,1,1],
                [1,1,1],
                [1,1,1],
            ],
            dtype=np.uint8,
            )
        
        nr_raincells_nearby = lfr.focal_sum(raincells, kernel)
               
        sRain = nr_raincells_nearby == 1
        mRain = nr_raincells_nearby == 2 | 3
        lRain = nr_raincells_nearby > 3
        
        random = max(0, np.random.randint(0, 20) - 10)
        rain_value = lfr.uniform(self.lue_zero(), np.float32, 0, 0.3)
        
        rain = lfr.where(sRain, rain_value, 0)
        rain = lfr.where(mRain, rain_value*2, rain)
        rain = lfr.where(lRain, rain_value*4, rain)
        
        return rain 
        
    def temp(self, temperature):
        """
        Assign a random value between -5 and 35 to the region as a temperature indication \
        if no temperature was assigned prior. Otherwise change the temperature with a tenth \
        of the assigned value to allow for slight changes in the system.
        """
        if temperature == None:
            temperature = np.random.randint(0,15)
        else:
            temperature = temperature + np.random.randint(-5,5)/10
        
        return temperature
        
    def pet_Hamon(self, temperature, date):
        """
        Simulates the potential evaporation based on parameters supplied.
        As of yet only takes the temperature into account, makes sure that any temperature \
        below zero has no evporation, the rest is exponentially scaled with temperature. 
        """
        latitude = 0.90                                                     # in radians
        day = int((date - self.start_of_year).days) + 1                     # calculate the day of the year
        declination = 1 + 0.033*math.cos((2*math.pi*day)/365)               # calculate the declination
        sunset_angle = math.acos(math.radians(-math.tan(math.radians(declination))*math.tan(math.radians(latitude))))
        N = (24/math.pi)*sunset_angle                                       # calculate the amount of daylight hours per 12 hours
        K = 273.3                                                           # Kelvin to degrees conversion
        k = 1                                                               # Proportionality constant
        
        saturated_vapor_pressure = 6.108*math.e**((17.27*temperature)/(temperature + K))        # calculate saturated vapor pressure
        evaporation = k * 0.165 * 216.7 * N * (saturated_vapor_pressure / (temperature + K))    # calculate the evaporation rate in mm / day
        evaporation = evaporation / 1000                                                        # convert to m / day
        evaporation = lfr.create_array(config.arrayShape,
                                       config.partitionShape,
                                       dtype = np.float32,
                                       fill_value = np.float32(evaporation),
                                       )
        return evaporation
    
     
    @lfr.runtime_scope 
    def simulate(self):
        if config.generateDEM:
            dem = self.dem("dem")
            ldd = self.ldd(dem, "ldd")
        
        raincells = self.initialize_raincells()                         # Assigns random cells that will produce rain
        
        temperature = None                                              # Initialize temperature
        
        
        # Calculate the amount of rainfall and precipitation for each of the days in the timeperiod
        for i in range(int((config.endDate - config.startDate).days)):
            print(f'Generating data for day {i + 1}.')
            # Rainfall generation per day
            if config.generatePrecip:
                if i < 5:
                    rain = self.precipitation(raincells)
                    lfr.to_gdal(rain, f'{self.path}precipitation_{config.arrayExtent}_{date}.tiff')
                elif 12 > i > 8:
                    rain = self.precipitation(raincells)
                    lfr.to_gdal(rain, f'{self.path}precipitation_{config.arrayExtent}_{date}.tiff')
                elif 23 > i > 18:
                    rain = self.precipitation(raincells)
                    lfr.to_gdal(rain, f'{self.path}precipitation_{config.arrayExtent}_{date}.tiff')
                else:
                    rain = self.lue_zero()
                    lfr.to_gdal(rain, f'{self.path}precipitation_{config.arrayExtent}_{date}.tiff')
            
            if config.generateEvapo:
                # First generate the temperature
                temperature = self.temp(temperature)
                temperature_2d = lfr.create_array(config.arrayShape,
                                           config.partitionShape,
                                           dtype = np.float32,
                                           fill_value = temperature)
                
                if config.saveTemp:             # IF the temperature should be saved, save the temperature.
                    lfr.to_gdal(temperature_2d, f'{self.path}temperature_{config.arrayExtent}_{date}.tiff')
                
                # Generature the evaporation based on the temperature and save it.
                evaporation = self.pet_Hamon(temperature, date)
                lfr.to_gdal(evaporation, f'{self.path}potential_evaporation_{config.arrayExtent}_{date}.tiff')
            
            # To the next day!
            date = date + datetime.timedelta(days=1)
        
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
if __name__ == "__main__":
    main = generate()
    main.simulate()