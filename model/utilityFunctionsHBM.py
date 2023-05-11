# -*- coding: utf-8 -*-
"""
Created on 8 May 2023

@author: steven.hosper
"""

import lue.framework as lfr

class utilityFunctions:
    def calculate_sqrd_slope(slope, max: float, min: float):
        slope_sqrd  = lfr.sqrt(slope)
        slope_sqrd  = lfr.where(slope_sqrd < min, min, slope_sqrd)
        slope_sqrd  = lfr.where(slope_sqrd > max, max, slope_sqrd)
        
        return slope_sqrd