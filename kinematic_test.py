# -*- coding: utf-8 -*-
"""
Created on 7 Mar 2023

@author: steven.hosper
"""
import configuration as config
import lue.framework as lfr
import numpy as np

def kinematic_test(self, ldd, head, inflow):
    """
    Summary:
        Calculates the flux of material through each cell based on the ldd, inflow, velocity, channel_length and coefficients alpha & beta.
        
    Input:
        ldd: local drain direction map based on the Digital Elevation Map / current map height
        head: head height of the current situation
        inflow: currently only the precipitation
    
    Returns:
        kinematic [unit]: the kinematic wave flux for each cell
    """
    velocity = lfr.atan((head - lfr.downstream(ldd, head))/1)
        
    channel = lfr.create_array(config.arrayShape,
                                config.partitionShape,
                                dtype = np.float32,
                                fill_value = 1,
                                )
    
    kinematic = lfr.kinematic_wave(
            ldd,
            velocity,
            inflow,
            1.5,
            0.6,
            1.0,
            channel,
        )
    return kinematic