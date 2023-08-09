from qgis.core import QgsVectorLayer, QgsFeature, QgsGeometry, QgsPointXY
import math
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, Polygon

# Returns the centroid coords of the crater rim
def find_largest_shape_center(shapefile_path):
  
    # Read the shapefile using geopandas
    gdf = gpd.read_file(shapefile_path)
    
    # Find the largest shape by area
    largest_shape = gdf.geometry.area.idxmax()
    
    # Get the centroid of the largest shape
    center_point = gdf.geometry.centroid.iloc[largest_shape]
    radius = math.sqrt(gdf.geometry.area.iloc[largest_shape] / math.pi)
    return (int(center_point.x), int(center_point.y)), radius / 3.75

# draws the camembert shapes
def add_coordinates_to_shapefile(shapefile_path, coordinates):
  
    # Load the shapefile as a QGIS vector layer
    layer = QgsVectorLayer(shapefile_path, 'New Feature', 'ogr')
    if not layer.isValid():
        print("Error: Invalid layer!")
        return

    # Start an edit session
    layer.startEditing()

    # Create a new feature with the given coordinates
    feature = QgsFeature()
    feature.setGeometry(QgsGeometry.fromPolygonXY([[QgsPointXY(coordinates[0][0], coordinates[0][1]),
                                                    QgsPointXY(coordinates[1][0], coordinates[1][1]),
                                                    QgsPointXY(coordinates[2][0], coordinates[2][1]),
                                                    QgsPointXY(coordinates[3][0], coordinates[3][1])]]))
    feature.setFields(layer.fields())

    # Add the feature to the layer
    layer.addFeature(feature)

    # Commit the changes and save the layer
    layer.commitChanges()
    layer.updateExtents()

    print("Coordinates added successfully!")

# The shapefile containing the crater outline as a feature
shapefile_path = "/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/CSFD_regions_cp0480.shp"

# The shapefile on which to draw the camemberts
radial_regions = "/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/cp0480_camembert.shp"

center_coords, radius = find_largest_shape_center(shapefile_path)

# starting angle
theta = 0
# The angle of each region (must agree with num_regions)
alpha = 22.5
# Number of regions to be generated (must agree with alpha)
num_regions = 16
# How far out the camemberts should extend from the rim (in crater radii)
length = 2.5

for i in range(num_regions):
    c1 = center_coords
    c4 = center_coords
    for j in range(10):
        coordinates = [c1, ((radius * (j + 1)) * math.cos(math.radians(theta)) + center_coords[0], (radius * (j + 1)) * math.sin(math.radians(theta)) + center_coords[1]), ((radius * (j + 1)) * math.cos(math.radians(theta + 22.5)) + center_coords[0], (radius * (j + 1)) * math.sin(math.radians(theta + 22.5)) + center_coords[1]), c4]
        c1, c4 = coordinates[1], coordinates[2]
        if not j <= 2:
            coordinates = [c1, (((radius * (j + 1)) * length) * math.cos(math.radians(theta)) + center_coords[0], ((radius * (j + 1)) * length) * math.sin(math.radians(theta)) + center_coords[1]), (((radius * (j + 1)) * length) * math.cos(math.radians(theta + alpha)) + center_coords[0], ((radius * (j + 1)) * length) * math.sin(math.radians(theta + alpha)) + center_coords[1]), c4]
            add_coordinates_to_shapefile(radial_regions, coordinates)
            break
    theta += alpha
