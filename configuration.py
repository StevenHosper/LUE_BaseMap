# -*- coding: utf-8 -*-
"""
Created on Mon Mar 6 2023

@author: steven.hosper
"""
import datetime
import os

# Use local disk or network disk
network = False

# Use API
useAPI = False

# Configuration file
# scenarios : ["generated", "De Wupsel", "De Tol"]
scenario = "generated"

# Dates
startDate = datetime.date(year = 2023, month = 2, day = 23)
endDate   = datetime.date(year = 2023, month = 3, day = 25)

# Include processes
variables = ['precipitation', 'evaporation', 'infiltration']
includeEvaporation   = True
includePrecipitation = True
includeInfiltration  = True
includePercolation   = True
includeGroundFlow    = True

# Set values
groundWaterTable = -1.40
resolution       = 1

# Username / Password for API session
username = '__key__'
password = 'Cy0BNm8p.vpytC2vYPT9g7OKdgxvqggyV0k9zzJVy'

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
    
# API URLs
precipAPI = "730d6675-35dd-4a35-aa9b-bfb8155f9ca7"
evapAPI   = "e262dc03-f12b-4082-a4f4-7d0534e31fa4"
demAPI    = "a60ad336-c95b-4fb6-b852-96fc352ee808"