from qgis.gui import QgsMapTool
from PyQt5.QtWidgets import QApplication
from qgis.utils import iface
import matplotlib.pyplot as plt
from qgis.core import QgsVectorLayer


small_features = 0
median_features = 0
large_features = 0
areas = []


# Get the selected features
layer = iface.activeLayer()
selected_features = layer.selectedFeatures()
attribute_field = "Aprx_Area"

# Print the attributes of the selected features
for feature in selected_features:
    attribute_index = layer.fields().indexFromName(attribute_field)
    attribute_value = feature[attribute_index]
    if not attribute_value > 20000:
        if attribute_value > 25.0:
            if attribute_value > 79.0:
                areas.append(attribute_value)
            else:
                areas.append(attribute_value)
        else:
            areas.append(attribute_value)
    

bins = int(max(areas) - min(areas)) / 8
# Clear the selection
layer.removeSelection()

# Create a scatter plot using Matplotlib
fig, ax = plt.subplots()

ax.hist(areas, bins = int(bins))
plt.xlabel('Aproximate area (8m bins)')
plt.ylabel('Frequency')
plt.title('Selected Frequency Distributions')
plt.show()
