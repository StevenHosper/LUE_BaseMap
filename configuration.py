# -*- coding: utf-8 -*-
"""
Created on Mon Mar 6 2023

@author: steven.hosper
"""
# Configuration for the HydrologicBaseModel
import datetime
import os

### GENERAL SETTINGS ###
network = False                                                     # use network disk of NS [True] or local disk [False]
useAPI = False                                                      # use API functionality to access data
generateData = False                                                # generate new data to be used
scenario = "generated"                                              # scenarios : ["generated", "De Wupsel", "De Tol"]
startDate = datetime.date(year = 2023, month = 2, day = 23)
endDate   = datetime.date(year = 2023, month = 3, day = 25)


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
includeInfiltration  = True
includePercolation   = True
includeGroundFlow    = True

# Report
variables = ['waterheight', 'groundWaterHeight']
fluxes    = []

# Set values
initialWaterTable = 1.30
resolution       = 1

if scenario == "De Wupsel" or scenario == "De Tol":
    arrayExtent     = 5000
    assert arrayExtent > 0
    arrayShape      = 2 * (arrayExtent,)
    partitionExtent = 1000
    assert partitionExtent > 0
    partitionShape  = 2 * (partitionExtent,)
elif scenario == "generated":
    arrayExtent     = 1000
    assert arrayExtent > 0
    arrayShape      = 2 * (arrayExtent,)
    partitionExtent = 1000
    assert partitionExtent > 0
    partitionShape  = 2 * (partitionExtent,)
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
    path = root_path
    
output_path = path + f'output/{scenario}/'

# Create ID variables
concrete       = [2, 4, 6, 8, 10, 14, 15, 16, 25, 28, 35, 166, 253]
green          = [40, 41, 43, 44, 112, 157]
water          = [51, 254]
compacted      = [18]
other_road     = [29]