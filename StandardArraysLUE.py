# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 10:03:54 2023

@author: steven.hosper
"""

import lue.framework as lfr
import numpy as np
import requests
import datetime
import math as math
import configuration as config
import configuration_v2

class StandardArraysLUE(): 
    def __init__(self, configuration):
        """
        Initialize the class. 
        1) Set extent, input and output dir.
        """
        
        self.arrayExtent    = int(configuration.modelSettings['arrayExtent'])
        self.partitionExtent= int(configuration.modelSettings['partitionExtent'])
        self.outputDir      = configuration.generalSettings['outputDir']
        self.inputDir       = configuration.generalSettings['inputDir']
        
    
    def boundary_cell(self):
        array = np.ones((self.partitionExtent - 2, self.arrayExtent - 2), dtype= np.uint8)
        boundaryCell = np.pad(array, pad_width=1, mode='constant', constant_values=0)
        boundaryCell = lfr.from_numpy(boundaryCell, 2*(self.partitionExtent,))
        lfr.to_gdal(boundaryCell, self.outputDir + '/boundaryCell.tiff')
        return boundaryCell
    
    def zero(self):
        return lfr.create_array(2*(self.arrayExtent,),
                                2*(self.partitionExtent,),
                                dtype = np.dtype(np.float64),
                                fill_value = 0,
                                )
    
    def one(self):
        return lfr.create_array(2*(self.arrayExtent,),
                                2*(self.partitionExtent,),
                                dtype = np.dtype(np.float64),
                                fill_value = 1,
                                )

    def one_int(self):
        return lfr.create_array(2*(self.arrayExtent,),
                                2*(self.partitionExtent,),
                                dtype = np.dtype(np.uint8),
                                fill_value = 1,
                                )

    def ldd_sink(self):
        return lfr.create_array(2*(self.arrayExtent,),
                                2*(self.partitionExtent,),
                                dtype = np.dtype(np.uint8),
                                fill_value = 5,
                                )