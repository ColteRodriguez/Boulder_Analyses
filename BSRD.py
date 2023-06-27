import matplotlib.pyplot as plt
from qgis.core import QgsVectorLayer
from qgis.gui import QgsMapTool
from PyQt5.QtWidgets import QApplication
from qgis.utils import iface
import numpy as np


layer = iface.activeLayer()
selected_features = layer.selectedFeatures()

# Create empty lists to store attribute values
distances = []
sizes = []

def calculate_standard_error(data_set):
    n = len(data_set)
    standard_errors = []
    
    for data_point in data_set:
        standard_error = np.std(data_set) / np.sqrt(n)
        standard_errors.append(standard_error)
    
    return standard_errors

# Iterate over the features in the layer
for feature in selected_features:
    # Get the attribute values for distance and size
    if feature['Aprx_Area'] < 1000:
        distance = feature['Dist_Rim']
        size = feature['Aprx_Area']
    
        # Append the values to the respective lists
        distances.append(distance)
        sizes.append(size)

# Create a scatter plot using Matplotlib
fig, ax = plt.subplots()

ax.errorbar(distances, sizes, calculate_standard_error(sizes), fmt='o', linewidth=2, capsize=4)

plt.xlabel('Distance to crater rim (m)')
plt.ylabel('Aprox. Boulder Area (m^2)')
plt.title('Selected Size-Range Distributions')
plt.show()
