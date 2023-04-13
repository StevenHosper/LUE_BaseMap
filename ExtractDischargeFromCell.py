from osgeo import gdal
import configuration as config

# Open the GeoTIFF file
dataset = gdal.Open(config.path + f'output/{config.scenario}/1_discharge_2023-02-24_50.tiff')

# Get the raster band
band = dataset.GetRasterBand(1)

# Get the cell value at the given coordinates
x = 58
y = 272
cols = dataset.RasterXSize
rows = dataset.RasterYSize
transform = dataset.GetGeoTransform()
xOrigin = transform[0]
yOrigin = transform[3]
pixelWidth = transform[1]
pixelHeight = -transform[5]
col = int((x - xOrigin) / pixelWidth)
row = int((yOrigin - y ) / pixelHeight)
value = band.ReadAsArray(col,row,1,1)

# Print the cell value
print(value)