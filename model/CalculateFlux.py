# -*- coding: utf-8 -*-
"""
Created on 8 May 2023

@author: steven.hosper
"""

import lue.framework as lfr
from StandardArraysLUE import StandardArraysLUE

class CalculateFlux:
    def __init__(self, configuration):
        """
        Initialize the class. 
        1) Initialize the standard arrays
        2) Set the processes that should be included.
        """
        self.std_arr_lue        = StandardArraysLUE(configuration)
        self.iterationsData     = int(configuration.modelSettings['iterationsBeforeData'])
        self.timestep           = int(configuration.modelSettings['timestep'])
        self.include_precipitation        = configuration.generalSettings['includePrecipitation']
        self.include_evapotranspiration   = configuration.generalSettings['includeEvapotranspiration']
        self.include_infiltration         = configuration.generalSettings['includeInfiltration']
        self.include_interception         = configuration.generalSettings['includeInterception']
        self.include_percolation          = configuration.generalSettings['includePercolation']
        self.ca                           = int(configuration.modelSettings['resolution'])**2
        self.T                            = self.std_arr_lue.zero()
        self.k                            = 1.6
        
    
    def interception(self, int_stor, int_stor_max, pre, rev, th_f):
        """Calculate the interceptionStorage, evapotranspiration of the surface and precipitation that falls through the canopy.
        
        Uses information known about the storage capacity and throughfall fraction of the vegetation to determine how much 
        precipitation is turned into interception and how much falls through. Furthermore determines the evapotranspirative
        potential still left for the soil surface. If the vegetation has a lot of water stored, no evaporation will happen
        at the surface, because the evaporative potential is already met by the vegetation.
        
        Args:
            int_stor (lpa*):     Current water stored from interception in the vegetation (m3)
            int_stor_max (lpa*): Maximum amount of interception stored (m3)
            pre (lpa*):         Precipitation rate (m3/s)
            rev (lpa*):         Reference evapotranspiration rate (m3/s)
            th_f (lpa*):         Throughfall fraction (-)
        
        Returns:
            int_stor (lpa*):                    Updated water stored from interception in the vegetation
            precipitation (lpa*):               Precipitation rate that falls through the vegetation towards the soil
            evapotranspiration_surface (lpa*):  Evaporative potential left after interception evaporation
            
        lpa*: lue partitioned array
        """      
        # Determine the interception rate
        interception    = (self.std_arr_lue.one() - th_f) * pre
        
        # Determine if there is enough water in the canopy to meet all evaporative potential
        enough_water_int  = (interception + int_stor/self.iterationsData) > (rev * (self.std_arr_lue.one() - th_f)) 

        # If there is enough water, update the interception storage, otherwise set at zero
        int_stor         = lfr.where(enough_water_int, int_stor + (interception - (rev * (self.std_arr_lue.one() - th_f))) * self.iterationsData, 0)
        
        # If the store is exceeded, add to variable called interception storage surplus
        int_storSurplus  = lfr.where(int_stor > int_stor_max, int_stor - int_stor_max, 0)
        
        # Remove the surplus from the interception storage
        int_stor         = int_stor - int_storSurplus
        
        # Add the interception storage surplus back to the precipitation
        precipitation   = pre - (interception - int_storSurplus/self.iterationsData)
        
        evapotranspirationSurface = lfr.where(enough_water_int, (rev * th_f), rev  - (interception + int_stor/self.iterationsData) * th_f)
        return int_stor, precipitation, evapotranspirationSurface

    def evapotranspiration(self, pre, ev_s):
        """Determine the actual surface and soil evapotranspiration
        
        Args:
            pre (lpa*): Precipitation rate (m3/s)
            ev_s (lpa*): Evaporative potential left for the surface (m3/s)
        
        Returns:
            evapotranspiration_surface (lpa*):   Actual surface evapotranspiration rate (m3/s)
            evapotranspiration_soil (lpa*):      Actual soil evaporation rate (m3/s)
        
        lpa*: lue partitioned array
        """
        enough_water_e = pre > ev_s
        evapotranspiration_soil = lfr.where(enough_water_e, 0, ev_s - pre)
        evapotranspiration_surface = lfr.where(enough_water_e, ev_s, pre)
        
        return evapotranspiration_surface, evapotranspiration_soil
    
    def infiltration(self, sgw, max_sgw, Ks, prm, por, pre, ev_s):
        """Determine the direct and potential channel infiltration
        
        First determines the potential infiltration rate. Then determines if plenty of water
        is available or not. This determines the actual direct infiltration and the potential
        channel infiltration. Any potential infiltration that is not met with direct infiltration
        is turned into potential channel infiltration. 
        
        Args:
            sgw (lpa*):     Groundwater storage in m3
            max_sgw (lpa*):  Maximum groundwater storage (m3)
            Ks (lpa*):      Hydraulic conductivity (m/s)
            prm (lpa*):     Permeability coefficient (-)
            por (lpa*):     Porosity (-)
            pre (lpa*):     Precipitation rate reaching the soil (m3/s)
            ev_s (lpa*):    Actual surface evapotranspiration (m3/s)
        
        Returns:
            direct_infiltration (lpa*):     precipitation that directly infiltrates the soil
            pot_infil_channel (lpa*):  if the potential infiltration is not yet reached,
                                            water flowing into the cells could still infiltrate.
        
        lpa*: lue partitioned array
        """
        pot_infiltration = Ks * prm * self.ca                                                              # m3 that can infiltrate the soil
        pot_infiltration = lfr.where((max_sgw - sgw)*por < pot_infiltration, \
                                    (max_sgw - sgw)*por, pot_infiltration)                          # The amount that can infiltrate because of capacity times the area
        pot_infiltration = lfr.where(pot_infiltration < 0, 0, pot_infiltration)
        enough_water_inf   = pre - ev_s > pot_infiltration                                          # If there is more water on the surface available than can infiltrate
        direct_infiltration = lfr.where(enough_water_inf, pot_infiltration, pre - ev_s)             # Either the pot_infiltration will fully infiltrate
        return direct_infiltration                               # or the available water at the surface will.
    
    def adjusted_Horton(self, sgw, max_sgw, Ks, Ki, prm, por, pre, ev_s):
        """Created a version of the Horton infiltration that deals with dry periods.
        It returns to high infiltration rates over time with a relatively simple function.
        One of the drawbacks is that by using 'time' the intensity does not matter.
        So a very light rain could affect the high infiltration rate, which would not be used at all.
        It should however still significantly improve the model from the function used prior.

        Args:
            sgw (lpa*):     Groundwater storage in m3
            max_sgw (lpa*): Maximum groundwater storage (m3)
            Ks  (lpa*):     Saturated hydraulic conductivity (m/s)
            Ki  (lpa*):     Intial hydraulic conductivity (m/s)
            k   (lpa*):     Decay factor (-)
            prm (lpa*):     Permeability coefficient (-)
            por (lpa*):     Porosity (-)
            pre (lpa*):     Precipitation rate reaching the soil (m3/s)
            ev_s(lpa*):     Actual surface evapotranspiration (m3/s)

        Returns:
            infiltration (lpa*): the infiltration rate based on horton (m3/s)
        """
        
        wet = pre > ev_s
        
        # Horton infiltration
        pot_infiltration_hort = Ks + (Ki - Ks)*lfr.pow(self.std_arr_lue.one()*2.71828,(-1*self.k*self.T))
        pot_infiltration = pot_infiltration_hort * prm * self.ca     # m3 that can infiltrate the soil
        pot_infiltration = lfr.where((max_sgw -sgw)*por / (self.timestep * self.iterationsData) > pot_infiltration,
                                        pot_infiltration,
                                        (max_sgw -sgw)*por / (self.timestep * self.iterationsData)) 
        
        enough_rain = pre - ev_s > pot_infiltration
        
        self.T = lfr.where(wet, self.T + (self.iterationsData * self.timestep / 3600),
                           self.T - (self.iterationsData * self.timestep / 3600))
        
        self.T = lfr.where(self.T < 0, 0, self.T)
        
        infiltration1 = lfr.where(enough_rain, pot_infiltration, pre - ev_s)  
        infiltration2 = lfr.where(infiltration1 < 0, 0, infiltration1)
            
             
            
        pot_reinfiltration = pot_infiltration - infiltration2
            
        return infiltration2, pot_reinfiltration