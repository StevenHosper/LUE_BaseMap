# -*- coding: utf-8 -*-
"""
Created on Mon Mar 6 2023

@author: steven.hosper
"""
# Configuration for the HydrologicBaseModel
import datetime
import os

### GENERAL SETTINGS ###
network      = False                                                    # use network disk of NS [True] or local disk [False]
useAPI       = False                                                    # use API functionality to access data
scenario     = "De Hupsel"                                              # scenarios : ["generated", "De Hupsel", "De Tol", "De Hupsel10", "unitTest"]
startDate    = datetime.date(year = 2023, month = 2, day = 23)
endDate      = datetime.date(year = 2023, month = 2, day = 24)
dt           = 30                                                      # The time difference between reported data in seconds (for v2)
dT           = 2                                                        # The times new data will be loaded in.
v2 = True
unitTest = False


### API SETTINGS ###
# Username / Password for API session
username = '__key__'
password = 'Cy0BNm8p.vpytC2vYPT9g7OKdgxvqggyV0k9zzJVy'

# API URLs
precipAPI = "730d6675-35dd-4a35-aa9b-bfb8155f9ca7"
evapAPI   = "e262dc03-f12b-4082-a4f4-7d0534e31fa4"
demAPI    = "a60ad336-c95b-4fb6-b852-96fc352ee808"


### MODEL SETTINGS ###
# Include processes
includePrecipitation = True
includeEvaporation   = True
includeInfiltration  = False
includeInterception  = False
includePercolation   = False
includeGroundFlow    = False                                              # Currently not working

# Report
variables = ['waterheight', 'groundWaterHeight']
fluxes    = []

# Set values
initialWaterTable = 1.30
waterBelowDEM     = 0.05
imperviousLayer   = 2.00

if scenario == "De Hupsel" or scenario == "De Tol":
    arrayExtent     = 5000
    partitionExtent = 1000
    arrayShape      = 2 * (arrayExtent,)
    partitionShape  = 2 * (partitionExtent,)
    resolution       = 1
elif scenario == "generated":
    arrayExtent     = 1000
    partitionExtent = 1000
    arrayShape      = 2 * (arrayExtent,)
    partitionShape  = 2 * (partitionExtent,)
    resolution       = 1
elif scenario == "De Hupsel10":
    arrayExtent     = 2000
    partitionExtent = 1000
    arrayShape      = 2 * (arrayExtent,)
    partitionShape  = 2 * (partitionExtent,)
    resolution      = 1
elif scenario == "unitTest":
    arrayExtent = 10
    partitionExtent = 10
    arrayShape      = 2 * (arrayExtent,)
    partitionShape  = 2 * (partitionExtent,)
    resolution      = 1
else:
    raise("Invalid scenario")

# Root path
root_path = f"{os.path.dirname(__file__)}" 
if not useAPI:
    if network:
        path = f"{root_path}/"
    else:
        path = "C:/Users/steven.hosper/Desktop/Mapje Stage/"      # Local directory
else:
    path = "C:/Users/steven.hosper/Desktop/Mapje Stage/"
    
output_path = path + f'output/{scenario}/'

# Create ID variables
concrete       = [2, 4, 6, 8, 10, 13, 14, 15, 16, 25, 28, 35, 166, 253]
green          = [40, 41, 42, 43, 44, 112, 157]
water          = [51, 254]
compacted      = [18]
other_road     = [29]
total          = [2, 4, 6, 8, 10, 13, 14, 15, 16, 25, 28, 35, 166, 253, 40, 41, 42, 43, 44, 112, 157, 51, 254, 18, 29]

### GENERATE DATA SETTINGS ###
generateDEM     = False                                          # Simulate a new digital elevation map
generatePrecip  = True                                          # Simulate new precipitation data
generateEvapo   = True                                          # Simulate new evaporation data
saveTemp        = False                                          # Save the temperature GeoTiff, this determines the amount of variables and thus speed.