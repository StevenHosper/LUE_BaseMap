# -*- coding: utf-8 -*-
"""
Created on 8 May 2023

@author: steven.hosper
"""

import lue.framework as lfr
import datetime

class utilityFunctions:
    def calculate_sqrd_slope(slope, max: float, min: float):
        slope_sqrd  = lfr.sqrt(slope)
        slope_sqrd  = lfr.where(slope_sqrd < min, min, slope_sqrd)
        slope_sqrd  = lfr.where(slope_sqrd > max, max, slope_sqrd)
        
        return slope_sqrd
    
    def string_to_datetime(date_string: str, seperator: str):
        date_int_list = list(map(int, date_string.split(seperator)))
        datetime_date = datetime.datetime(date_int_list[0],
                                          date_int_list[1],
                                          date_int_list[2],
                                          date_int_list[3],
                                          date_int_list[4],
                                          date_int_list[5])
        return datetime_date