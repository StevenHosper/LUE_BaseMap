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
    
    def slice_pathname(pathname, idx):
        date = datetime.date(year = 2023, month = 2, day = 24)
        if config.v2:
            second = (idx) * config.timestep
            date = f'{date + datetime.timedelta(seconds=idx)}_{second}'
        else:
            date = date + datetime.timedelta(idx)
        return "{}_{}.tiff".format(pathname, date)


    def read_raster(raster_pathname, idx):
        dataset = gdal.Open(makeGIF.slice_pathname(raster_pathname, idx))
        return np.array(dataset.GetRasterBand(1).ReadAsArray())


    def create_animation(raster_pathname, nr_rasters, animation_pathname, vmin, vmax):
        with iio.get_writer(animation_pathname, mode="i", fps=10) as writer:
            for i in range(nr_rasters + 1):
                figure, axis = plt.subplots(figsize=(10, 10))
                axis.set_axis_off()
                data = makeGIF.read_raster(raster_pathname, i)
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
                
    def create_animation2(raster_pathname, nr_rasters, animation_pathname, vmin, vmax):
        with iio.get_writer(animation_pathname, mode="i", fps=10) as writer:
            for i in range(nr_rasters + 1):
                figure, axis = plt.subplots(figsize=(10, 10))
                axis.set_axis_off()
                data = makeGIF.read_raster(raster_pathname, i)
                image = rasterio.plot.show(
                    data,
                    ax=axis,
                    cmap="magma",
                    vmin = vmin,
                    vmax = vmax,
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
    nr_rasters  = int(configuration.gifSettings['nrRasters'])
    timestep    = int(configuration.modelSettings['timestep'])
    assert nr_rasters >= 0
    
    # Create animations
    for var in variables:
        raster_pathname = f'{path}/{timestep}_{var}'
        animation_pathname = f'{path}/{var}.gif'
        print(raster_pathname)
        assert not os.path.splitext(raster_pathname)[1]
        
        vmin = vminDict[var]
        vmax = vmaxDict[var]
        
        makeGIF.create_animation(raster_pathname, nr_rasters, animation_pathname, vmin, vmax)


if __name__ == "__main__":
    configuration = Configuration("F:/Projecten intern (2023)/Stage Steven Hosper/Model/v1/config.ini")
    run(configuration)