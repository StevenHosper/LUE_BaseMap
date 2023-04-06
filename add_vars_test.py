# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 2023

@author: steven.hosper
"""
# The HydrologicBaseModel
import lue.framework as lfr
import math as math
import os
import sys
import time
import numpy as np

# Own scripts
import configuration as config

# Timer to add some measure of functionality to the program
start_time = time.time()

usage = """\
Run the main model of the hydrologic base model.

Usage:
    {command}

Options:
    {command} : --hpx:thread = integer;
                The integer is the amount of cores used during the model run.
""".format(
    command=os.path.basename(sys.argv[0])
)


class mainModel():
    def dynamicModel(self):
        a = lfr.create_array(config.arrayShape,
                             config.partitionShape,
                             dtype = np.float32,
                             fill_value = 2.0)
        
        b = lfr.create_array(config.arrayShape,
                             config.partitionShape,
                             dtype = np.float32,
                             fill_value = 10.0)
        
        c = lfr.create_array(config.arrayShape,
                             config.partitionShape,
                             dtype = np.float32,
                             fill_value = 30.0)
        
        height = lfr.create_array(config.arrayShape,
                             config.partitionShape,
                             dtype = np.float32,
                             fill_value = 0.0)
        
        for i in range (5):
            for j in range(60):
                height = height + a - b - c
                
            input("Press Enter to continue...")
                
# Initialize HPX runtime and run model, on the root locality -------------------
# General configuration options, which are valid on all
# platforms. Platform-specific options can be passed on the command line.
cfg = [
    # Make sure hpx_main is always executed
    "hpx.run_hpx_main!=1",
    # Allow for unknown command line options
    "hpx.commandline.allow_unknown!=1",
    # Disable HPX' short options
    "hpx.commandline.aliasing!=0",
    # Don't print diagnostics during forced terminate
    "hpx.diagnostics_on_terminate!=1",
    # Make AGAS clean up resources faster than by default
    "hpx.agas.max_pending_refcnt_requests!=50",
]

lfr.start_hpx_runtime(cfg)

# The root locality will distribute the work over all other
# localities. Never perform Python code on the other localities than the
# root locality unless you know what you are doing.
if lfr.on_root_locality():
    # Run the main model
    main = mainModel()
    main.dynamicModel()
        
        