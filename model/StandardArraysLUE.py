# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 10:03:54 2023

@author: steven.hosper
"""

import lue.framework as lfr
import numpy as np
import math as math

class StandardArraysLUE: 
    def __init__(self, configuration):
        """
        Initialize the class. 
        1) Set extent, input and output dir.
        """
        
        self.array_extent    = int(configuration.modelSettings['arrayExtent'])
        self.partition_extent= int(configuration.modelSettings['partitionExtent'])
        self.output_dir      = configuration.generalSettings['outputDir']
        
    
    def boundary_cell(self):
        array = np.ones((self.partition_extent - 2, self.array_extent - 2), dtype= np.uint8)
        boundary_cell = np.pad(array, pad_width=1, mode='constant', constant_values=0)
        boundary_cell = lfr.from_numpy(boundary_cell, 2*(self.partition_extent,))
        lfr.to_gdal(boundary_cell, self.output_dir + '/boundary_cell.tiff')
        return boundary_cell
    
    def zero(self):
        return lfr.create_array(2*(self.array_extent,),
                                2*(self.partition_extent,),
                                dtype = np.dtype(np.float64),
                                fill_value = 0,
                                )
    
    def one(self):
        return lfr.create_array(2*(self.array_extent,),
                                2*(self.partition_extent,),
                                dtype = np.dtype(np.float64),
                                fill_value = 1,
                                )

    def one_int(self):
        return lfr.create_array(2*(self.array_extent,),
                                2*(self.partition_extent,),
                                dtype = np.dtype(np.uint8),
                                fill_value = 1,
                                )

    def ldd_sink(self):
        return lfr.create_array(2*(self.array_extent,),
                                2*(self.partition_extent,),
                                dtype = np.dtype(np.uint8),
                                fill_value = 5,
                                )