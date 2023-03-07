# -*- coding: utf-8 -*-
"""
Created on Mon Mar 6 2023

@author: steven.hosper
"""
import datetime
import os

# Username / Password for API session
username = '__key__'
password = 'Cy0BNm8p.vpytC2vYPT9g7OKdgxvqggyV0k9zzJVy'

# Configuration file
arrayExtent     = 5000
partitionExtent = 1000
resolution      = 1

# Root path
path = f"{os.path.dirname(__file__)}" 

# Dates
startDate = datetime.date(year = 2023, month = 2, day = 23)
endDate = datetime.date(year = 2023, month = 3, day = 25)

# Data source
useApi = False

# Include processes
includeEvaporation = True
includePrecipitation = True
includeInfiltration = True
includePercolation = True
includeGroundFlow = True

# Set values
groundWaterTable = -1.40

# Use API
useAPI = False

# API URLs
precipAPI = "730d6675-35dd-4a35-aa9b-bfb8155f9ca7"
evapAPI   = "e262dc03-f12b-4082-a4f4-7d0534e31fa4"
demAPI    = "a60ad336-c95b-4fb6-b852-96fc352ee808"