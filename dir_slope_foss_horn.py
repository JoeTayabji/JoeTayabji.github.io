import time

start_time = time.time()
'''
JOSEPH TAYABJI GEOG 6150

This script is designed to taken two inputs: a DEM raster and a line. The output is a 
line. The input line is resampled at a specified interval (default = 0.5 m), and the 
directional slope is calculated and added as a field at each segment. 

For standalone use the interpreter should be (or similar to): 
C:\Program Files\QGIS 3.32.1\apps\Python39\python3.exe
User may need to edit system PATH and USER variables, instructions can be found at 
(first 8 steps):
https://guides.lib.utexas.edu/gis/python-qgis-scripting#:~:text=Identify%20the%20path%20to%20the,of%20the%20python.exe%20file.
'''

import os
from math import *

# import necessary libraries and modules
# note: some imports will not be necessary inside the Q gui
import numpy as np
from osgeo import gdal, ogr
from shapely import wkt
from shapely.geometry import MultiLineString, Point

gdal.UseExceptions() # other option is ogr.DontUseExceptions()

# load DEM and necessary metadata

'''
YOUR DEM TIF PATH SHOULD BE ASSIGNED HERE
'''
dem_path = r"<YOUR DEM PATH>"
dem_raster = gdal.Open(dem_path)
# get geotransform coefficients of raster to allow translation between raster and pixel coordinates
dem_gt = dem_raster.GetGeoTransform()
# inverse geotransform will be used later to extract raster values at the resampled points
dem_reverse_gt = gdal.InvGeoTransform(dem_gt)
if dem_raster is not None:
    print("DEM was loaded.")
else:
    print("DEM was not loaded.")
cell_size = dem_gt[1]   
print("The DEM's spatial resolution is {0} m".format(cell_size))

# set file directories for line shapefiles
'''
THE SAVE PATH SHOULD BE SET HERE
'''
out_dir = r"<FOLDER TO SAVE TO>"
os.makedirs(out_dir, exist_ok=True) # Create the output directory if it doesn't exist
'''
THE FULL PATH (INCLUDING .shp FILE) SHOULD BE SET HERE
'''
line_in_path = r"<INPUT FILE>"
'''
DESIRED FILE NAME SHOULD BE SET HERE
'''
shp_out_name = "output_dir_slope_line.shp" # output shapefile name
line_out_path = os.path.join(out_dir, shp_out_name)
line_in_ds = ogr.Open(line_in_path, 0) # 0 is read only, 1 is write
if line_in_ds is not None:
    print("Input line was loaded.")
else:
    print("Input line was not loaded")

# get required line metadata
line_in_layer = line_in_ds.GetLayer()
line_in_spref = line_in_layer.GetSpatialRef()
line_in_layer_def = line_in_layer.GetLayerDefn() # schema of input line from the .dbf
num_fields = line_in_layer_def.GetFieldCount()
field_name = [line_in_layer_def.GetFieldDefn(i).GetName() for i in range(num_fields)]

# create the directional slope line (output)
line_out_ds = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(line_out_path)
line_out_layer = line_out_ds.CreateLayer("Dir Slope Line", geom_type = ogr.wkbLineString, srs = line_in_spref)

# get input field names and copy them to the output line
for i in range(num_fields):
    field_defn = line_in_layer_def.GetFieldDefn(i)
    new_field = ogr.FieldDefn(field_defn.GetName(), field_defn.GetType())
    line_out_layer.CreateField(new_field)
# add a bearing field
bearing_field = ogr.FieldDefn("Bearing", ogr.OFTReal)
line_out_layer.CreateField(bearing_field)
# add a directional slope field
dir_slope_field = ogr.FieldDefn("Dir_Slope", ogr.OFTReal)
line_out_layer.CreateField(dir_slope_field)

def get_horn(px, py):
    '''
    The get_horn function is used to calculate slope at one point's location.
    Horn's method (Horn, 1981) is used, which is preferrential for variable terrain.
    Horn's method is also the default in ESRI products.
    The variables a through i reference cells in a 3x3 grid centered around e.
    This convention is used in ESRI's Slope documentation, found at:

    Parameters:
    px is the x coordinate of the input point, using pixel coordinates.
    py is the y coordinate of the input point, using pixel coordinates.

    Returns:
    Two values are returned: dz_dx, and dz_dy. These can be used with a bearing
    value to calculate directional slope.
    '''

    # Get a 3x3 NumPy array with DEM values centered at location of px, py
    local_dem = dem_raster.ReadAsArray(px-1, py-1, 3, 3)
    # top left = a
    a = local_dem[0, 0]
    # top center = b
    b = local_dem[0, 1]
    # top right = c
    c = local_dem[0, 2]
    # left = d
    d = local_dem[1, 0]
    # the pixel = e
    e = local_dem[1, 1]
    # right = f
    f = local_dem[1, 2]
    # bottom left = g
    g = local_dem[2, 0]
    # bottom center = h
    h = local_dem[2, 1]
    # bottom right = i
    i = local_dem[2, 2]

    # right column - left column
    dz_dx = ((c+(2*f)+i)-(a+(2*d)+g))/(8*cell_size)
    # top row - bottom row
    dz_dy = ((c+(2*b)+a)-(i+(2*h)+g))/(8*cell_size)

    return dz_dx, dz_dy

def dir_slope(x, y, previous_x, previous_y):
    '''
    The dir_slope function calculates the directional slope and adds it to the new line.
    It requires the x and y coordinates of two points (must be in a project coordinate system).
    The two points should be consecutive points along the input point, separated by the specified distance.
    First the coordinates of the current point are transformed to pixel coordinates and the two slope
    values, dz_dx and dz_dy, are retrieved from the get_horn function.
    Then, the direction (bearing) is calculated using the two points.
    Using this these three values, the directional slope can be calculated (Neteler and Mitášová, 2008).
    Negative values indicate uphill, positive values indicate downhill.
    The bearing and directional slope values are added to the associated fields, and the
    remaining fields are copied from the input line. Finally, this new line segmented is added to the layer.

    Parameters: 
    x, y, previous_x, and previous_y are the x and y coordinates of the current and previous point.
    '''

    # transform x, y coordinates to pixel coordinates
    px, py = gdal.ApplyGeoTransform(dem_reverse_gt, x, y)
    px = floor(px)
    py = floor(py)
    # extract the dzdx and dzdy values
    dx, dy = get_horn(px, py)

    # calculate bearing in radians, must be a projected coordinate system
    delta_x = x - previous_x
    delta_y = y - previous_y
    bearing = atan2(delta_x, delta_y) 
    bearing = (bearing + 2*pi) % (2*pi)
    tan_alpha = ((dx*sin(bearing)) + (dy * cos(bearing)))
    dir_slope = atan(tan_alpha)*(-180.0)/pi # this is the directional slope value!

    # add line to layer
    line_out_geom = ogr.Geometry(ogr.wkbLineString)
    line_out_geom.AddPoint(previous_x, previous_y)
    line_out_geom.AddPoint(x, y)
    line_out_feat = ogr.Feature(line_out_layer.GetLayerDefn())
    line_out_feat.SetGeometry(line_out_geom)

    # Copy attributes from the input line to the output line
    for fn in field_name:
        line_out_feat.SetField(fn, p.GetField(fn))
    # add bearing and directional slope fields, add line to layer
    line_out_feat.SetField("Bearing", degrees(bearing))
    line_out_feat.SetField("Dir_Slope", dir_slope)
    line_out_layer.CreateFeature(line_out_feat)

'''
Get input line length and get points along it separated by the desired interval. 
At each point, calculate the directional slope (and add line to layer).
'''
sample_interval = 0.5 # resampling distance
for p in line_in_layer:
    geom = p.GetGeometryRef().ExportToWkt()
    shapely_line = wkt.loads(geom)
    # Point at distance 0
    first_point = Point(list(shapely_line.coords)[0])
    previous_x = first_point.xy[0][0]
    previous_y = first_point.xy[1][0]

    # Go through the input line in increments equal to sample_interval
    line_length = shapely_line.length
    print("The line length is {:.2f} m.".format(line_length))
    distance = sample_interval # initialize distance
    while distance < line_length:
        # get point coordinates for current location on the input line
        point = shapely_line.interpolate(distance)
        x, y = point.xy[0][0], point.xy[1][0]

        dir_slope(x, y, previous_x, previous_y)

        # set up for next iteration
        distance += sample_interval
        previous_x = x
        previous_y = y
    
    # do the above steps one last time for the last point on the line
    last_point = Point(list(shapely_line.coords)[-1])
    x, y = last_point.xy[0][0], last_point.xy[1][0]
    dir_slope(x, y, previous_x, previous_y)
    print("Output line created.")    
line_out_ds = None

end_time = time.time()
total_time = end_time - start_time
print("Success. This script took {:.2f} seconds".format(total_time))