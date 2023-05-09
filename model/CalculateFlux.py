# -*- coding: utf-8 -*-
"""
Created on 8 May 2023

@author: steven.hosper
"""

import lue.framework as lfr
import numpy as np
import datetime
import pandas as pd
from StandardArraysLUE import StandardArraysLUE

class CalculateFlux():
    def __init__(self, configuration):
        """
        Initialize the class. 
        1) Initialize the standard arrays
        2) Set the processes that should be included.
        """
        self.std_arr_lue        = StandardArraysLUE(configuration)
        
        self.iterations = int(configuration.modelSettings['iterationsBeforeReport'])
        self.includePrecipitation        = configuration.generalSettings['includePrecipitation']
        self.includeEvapotranspiration   = configuration.generalSettings['includeEvapotranspiration']
        self.includeInfiltration         = configuration.generalSettings['includeInfiltration']
        self.includeInterception         = configuration.generalSettings['includeInterception']
        self.includePercolation          = configuration.generalSettings['includePercolation']
        self.cA                    = int(configuration.modelSettings['resolution'])**2
        
    
    def interception(self, intStor, intStorMax, pre, rev, thF):
        """Calculate the interceptionStorage, evapotranspiration of the surface and precipitation that falls through the canopy.
        
        Uses information known about the storage capacity and throughfall fraction of the vegetation to determine how much 
        precipitation is turned into interception and how much falls through. Furthermore determines the evapotranspirative
        potential still left for the soil surface. If the vegetation has a lot of water stored, no evaporation will happen
        at the surface, because the evaporative potential is already met by the vegetation.
        
        Args:
            intStor (lpa*):     Current water stored from interception in the vegetation
            intStorMax (lpa*):  Maximum amount of interception stored
            pre (lpa*):         Precipitation rate
            rev (lpa*):         Reference evapotranspiration rate
            thF (lpa*):         Throughfall fraction
        
        Returns:
            intStor (lpa*):                     Updated water stored from interception in the vegetation
            precipitation (lpa*):               Precipitation rate that falls through the vegetation towards the soil
            evapotranspirationSurface (lpa*):   Evaporative potential left after interception evaporation
            
        lpa*: lue partitioned array
        """      
        # Determine the interception rate
        interception    = (self.std_arr_lue.one() - thF) * pre
        
        # Determine if there is enough water in the canopy to meet all evaporative potential
        enoughWaterInt  = (interception + intStor/self.iterations) > rev  

        # If there is enough water, update the interception storage, otherwise set at zero
        intStor         = lfr.where(enoughWaterInt, intStor + (interception - rev) * self.iterations, 0)
        
        # If the store is exceeded, add to variable called interception storage surplus
        intStorSurplus  = lfr.where(intStor > intStorMax, intStor - intStorMax, 0)
        
        # Remove the surplus from the interception storage
        intStor         = intStor - intStorSurplus
        
        # Add the interception storage surplus back to the precipitation
        precipitation   = pre - (interception - intStorSurplus/self.iterations)
        
        evapotranspirationSurface = lfr.where(enoughWaterInt, 0, rev - (interception + intStor/self.iterations))
        return intStor, precipitation, evapotranspirationSurface

    def evapotranspiration(self, pre, evS):
        """Determine the actual surface and soil evapotranspiration
        
        Args:
            pre (lpa*): Precipitation rate
            evS (lpa*): Evaporative potential left for the surface
        
        Returns:
            evapotranspirationSurface (lpa*):   Actual surface evapotranspiration rate
            evapotranspirationSoil (lpa*):      Actual soil evaporation rate
        
        lpa*: lue partitioned array
        """
        enoughWaterE = pre > evS
        evapotranspirationSoil = lfr.where(enoughWaterE, 0, evS - pre)
        evapotranspirationSurface = lfr.where(enoughWaterE, evS, pre)
        
        return evapotranspirationSurface, evapotranspirationSoil
    
    def infiltration(self, Sgw, MaxSgw, Ks, prm, por, pre, evS):
        """Determine the direct and potential channel infiltration
        
        First determines the potential infiltration rate. Then determines if plenty of water
        is available or not. This determines the actual direct infiltration and the potential
        channel infiltration. Any potential infiltration that is not met with direct infiltration
        is turned into potential channel infiltration. 
        
        Args:
            Sgw (lpa*):     Groundwater storage
            MaxSgw (lpa*):  Maximum groundwater storage
            Ks (lpa*):      Hydraulic conductivity
            prm (lpa*):     Permeability
            por (lpa*):     Porosity
            pre (lpa*):     Precipitation rate reaching the soil
            evS (lpa*):     Actual surface evapotranspiration
        
        Returns:
            direct_infiltration (lpa*):     precipitation that directly infiltrates the soil
            potInfiltrationChannel (lpa*):  if the potential infiltration is not yet reached,
                                            water flowing into the cells could still infiltrate.
        
        lpa*: lue partitioned array
        """
        potInfiltration = Ks * prm                                                              # meters that can infiltrate the soil
        potInfiltration = lfr.where(((MaxSgw - Sgw)/self.cA)*por < potInfiltration, \
                                    ((MaxSgw - Sgw)/self.cA)*por, potInfiltration) * self.cA    # The amount that can infiltrate because of capacity times the area
        enoughWaterInf   = pre - evS > potInfiltration                                          # If there is more water on the surface available than can infiltrate
        direct_infiltration = lfr.where(enoughWaterInf, potInfiltration, pre - evS)             # Either the potInfiltration will fully infiltrate
        potInfiltrationChannel = lfr.where(enoughWaterInf, 0, 
                                           potInfiltration - direct_infiltration) / self.cA
        return direct_infiltration, potInfiltrationChannel                                      # or the available water at the surface will.
        