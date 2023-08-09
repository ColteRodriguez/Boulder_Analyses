from qgis.core import QgsVectorLayer, QgsFeature, QgsGeometry, QgsPointXY
# from qgis.PyQt.QtCore import QVariant
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
    return (int(center_point.x), int(center_point.y)), radius

# draws the camembert and buffer shapes
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

# Shapefile on which to draw the ROIs
shapefile_path = "/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/CSFD_regions_cp1701.shp"
# Shapefile containing the crater perimeter as a feature
shapefile2_path = "/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/cpcrater1701/shp/inputs/M169764515RE-cpcrater1701-boulder-mapping.shp"

center_coords, radius = find_largest_shape_center(shapefile2_path)
print(center_coords)

# Starting angle
theta = 0
# angle of offset (must agree with num_regions)
alpha = 11.25
# Number of radian regions to be drawn (must agree with alpha)
num_regions = 32
# Length of each ROI relative to the crater radius
radius = radius/2

for i in range(num_regions):
    c1 = center_coords
    c4 = center_coords
    for j in range(7):
        coordinates = [c1, ((radius * (j + 1)) * math.cos(math.radians(theta)) + center_coords[0], (radius * (j + 1)) * math.sin(math.radians(theta)) + center_coords[1]), ((radius * (j + 1)) * math.cos(math.radians(theta + alpha)) + center_coords[0], (radius * (j + 1)) * math.sin(math.radians(theta + alpha)) + center_coords[1]), c4]
        c1, c4 = coordinates[1], coordinates[2]
        if not j <= 1:
            add_coordinates_to_shapefile(shapefile_path, coordinates)
    theta +=alpha
