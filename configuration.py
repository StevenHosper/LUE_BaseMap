# -*- coding: utf-8 -*-
"""
Created on Mon Mar 6 2023

@author: steven.hosper
"""
import datetime

# Configuration file
arrayExtent = 5000
partitionExtent = 1000

# Data and file path
path = "F:/Projecten intern (2023)/Stage Steven Hosper/Model" 

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
