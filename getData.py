# -*- coding: utf-8 -*-
"""
Created on 7 Mar 2023

@author: steven.hosper
"""

import lue.framework as lfr
import configuration
import numpy as np

class getData():
    # Create useful single value arrays
    # Zero, for empty cells or unincluded variables
    zero = lfr.create_array(configuration.arrayShape,
                            configuration.partitionShape,
                            dtype = np.dtype(np.float32),
                            fill_value = 0,
                            )
    # One, for additions
    ones = lfr.create_array(configuration.arrayShape,
                            configuration.partitionShape,
                            dtype=np.float32,
                            fill_value=1,
                            )
        
        # For the ldd direction being towards the cell itself (pit)
    sink = lfr.create_array(configuration.arrayShape,
                            configuration.partitionShape,
                            dtype = np.dtype(np.uint8),
                            fill_value = 5,
                            )