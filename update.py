# -*- coding: utf-8 -*-
"""
Created on 7 Mar 2023

@author: steven.hosper
"""
import lue.framework as lfr
import configuration as config

class update():
    def runoff(precipitation, evaporation,\
                    infiltration):
        runoff = precipitation - evaporation - infiltration
        return runoff
    
    def groundWaterHeight(dem, Ks, waterHeight, groundWaterHeight,\
                          infiltration, percolation, zero):
        # Determine the way water flows.
        # The hydraulic head is the groundWaterHeight plus the surfaceWaterHeight.
        groundhead = waterHeight + groundWaterHeight
        groundWaterLDD = lfr.d8_flow_direction(groundhead)
        
        # Calculate the groundwater flow
        if config.includeGroundFlow:
            # Currently not working
            groundflow = None                              # Q = k * i * A
        else:
            groundflow = zero
        
        # Add the groundwaterflow to each cell
        groundWaterHeight = groundWaterHeight + lfr.upstream(groundWaterLDD, groundflow) - groundflow
        
        # If more water flows to a cell than fits, it seeps out at the surface.
        seepage = groundWaterHeight - dem
        
        groundWaterHeight = groundWaterHeight + infiltration - percolation
    
        return groundWaterHeight
    
    def waterRouting():
        pass