import sys
import lue.framework as lfr

sys.path.insert(0, 'F:/Projecten intern (2023)/Stage Steven Hosper/Model/v1')
import DataGeneration as dg

dem = lfr.from_gdal("F:/Projecten intern (2023)/Stage Steven Hosper/Model/v1/data/dem/dem_1000x1000.tiff", (1000, 1000))
dg.dataGeneration.ldd_simulation(dg, dem, "ldd_1000x1000")