#!/usr/bin/env python
import docopt
import imageio.v2 as iio
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import rasterio.plot
import io
import osgeo.gdal as gdal
import os.path
import sys
import datetime
import configuration as config
from configuration_v2 import Configuration

class makeGIF():
    def __init__(self):
        pass
    
    def slice_pathname(pathname, idx, date):
        date = date + datetime.timedelta(minutes=idx)
        dateTime = date.strftime("%Y-%m-%d-%H%M")
        return "{}_{}.tiff".format(pathname, dateTime)


    def read_raster(raster_pathname, idx, dateTime):
        dataset = gdal.Open(makeGIF.slice_pathname(raster_pathname, idx, dateTime))
        return np.array(dataset.GetRasterBand(1).ReadAsArray())


    def create_animation(raster_pathname, nr_rasters, animation_pathname, vmin, vmax, date, FPS):
        with iio.get_writer(animation_pathname, mode="i", fps = FPS) as writer:
            for i in range(nr_rasters + 1):
                figure, axis = plt.subplots(figsize=(10, 10))
                axis.set_axis_off()
                data = makeGIF.read_raster(raster_pathname, i, date)
                image = rasterio.plot.show(
                    data,
                    ax=axis,
                    cmap="magma",
                    norm = colors.LogNorm(vmin, vmax),
                )

                with io.BytesIO() as buffer:
                    figure.savefig(buffer, format="raw")  # , bbox_inches="tight")
                    buffer.seek(0)
                    data = np.frombuffer(buffer.getvalue(), dtype=np.uint8)
                    nr_cols, nr_rows = figure.canvas.get_width_height()
                    image = data.reshape(nr_rows, nr_cols, -1)

                writer.append_data(image)
                plt.close()
                print(f'Processed {i} out of {nr_rasters} rasters.')


usage = """\
Create animation given a set of rasters stored in Geotiffs

Usage:
    {command} 
    
Options:

""".format(
    command=os.path.basename(sys.argv[0])
)


def run(configuration):
    argv = [arg for arg in sys.argv[1:] if not arg.startswith("--hpx")]
    arguments = docopt.docopt(usage, argv)
    
    # Use given path
    path = configuration.generalSettings['outputDir'] + configuration.generalSettings['scenario']
    
    # Determine variables and values
    variables   = configuration.gifSettings['variables'].split(", ")
    vminList    = configuration.gifSettings['vmin'].split(", ")
    vmaxList    = configuration.gifSettings['vmax'].split(", ")
    
    # Combine
    vminDict    = dict(zip(variables, vminList))
    vmaxDict    = dict(zip(variables, vmaxList))
    
    # Set amount of rasters and timestep
    sD = list(map(int, configuration.modelSettings['startDate'].split(", ")))
    startDate   = datetime.datetime(sD[0], sD[1], sD[2], sD[3], sD[4], sD[5])
    nr_rasters  = int(configuration.gifSettings['nrRasters'])
    timestep    = int(configuration.modelSettings['timestep'])
    fps         = int(configuration.gifSettings['fps'])
    assert nr_rasters >= 0
    
    # Create animations
    for var in variables:
        raster_pathname = f'{path}/{timestep}_{var}'
        animation_pathname = f'{path}/{var}.gif'
        assert not os.path.splitext(raster_pathname)[1]
        
        vmin = vminDict[var]
        vmax = vmaxDict[var]
        
        makeGIF.create_animation(raster_pathname, nr_rasters, animation_pathname, vmin, vmax, startDate, fps)


if __name__ == "__main__":
    configuration = Configuration("F:/Projecten intern (2023)/Stage Steven Hosper/Model/v1/config.ini")
    run(configuration)