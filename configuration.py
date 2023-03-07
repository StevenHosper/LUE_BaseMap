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
arrayExtent = 5000
partitionExtent = 1000

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
