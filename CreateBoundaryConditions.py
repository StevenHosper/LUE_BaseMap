import lue.framework as lfr
import pcraster as pcr
import configuration as config
import numpy as np
import os
import sys

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

@lfr.runtime_scope
def main():
    # Load raster
    dem = pcr.readmap("C:/Users/steven.hosper/Desktop/Mapje Stage/data/De Hupsel5/v2/demPCR.map")
    
    # Set all values to 1
    dem = (dem * 0) + 1
    demSum = pcr.windowtotal(dem, 9)
    pcr.report(dem, "C:/Users/steven.hosper/Desktop/Mapje Stage/data/De Hupsel5/v2/test.tiff")
    pcr.report(demSum, "C:/Users/steven.hosper/Desktop/Mapje Stage/data/De Hupsel5/v2/test2.tiff")

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

if __name__ == "__main__":
    main()