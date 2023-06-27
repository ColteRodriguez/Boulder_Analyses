from qgis.gui import QgsMapTool
from PyQt5.QtWidgets import QApplication
from qgis.utils import iface
import matplotlib.pyplot as plt
from qgis.core import QgsVectorLayer

# Distance bins will be every 57 meters, 10 in one crater radius
distances = []
densities = []
one_deg = 3033
density_window_size = 100
bounds = 1 / (one_deg / density_window_size)


# Get the selected features
layer = iface.activeLayer()
selected_features = layer.selectedFeatures()

# Bounding the selected region
lats = []
lons = []
for feature in selected_features:
    if not feature["Aprx_Area"] > 20000:
        lats.append(feature["Lat"])
        lons.append(feature["Lon"])
min_lon = min(lons) + 50
max_lon = max(lons) - 50
min_lat = min(lats) + 50
max_lat = max(lats) - 50

for feature in selected_features:
    if not feature["Aprx_Area"] > 20000:
        if min_lon < feature["Lon"] < max_lon:
            if min_lat < feature["Lat"] < max_lat:
                centeral_lat = feature["Lat"]
                centeral_lon = feature["Lon"]
                distance = feature["Dist_Rim"]
                
                count = 0
                for feature in selected_features:
                    if centeral_lat - 50 < feature["Lat"] < centeral_lat + 50:
                        if centeral_lon - 50 < feature["Lon"] < centeral_lon + 50:
                            count+=1
                distances.append(distance)
                densities.append(count)
  
  
# Create a scatter plot using Matplotlib
plt.scatter(distances, densities)
plt.xlabel('Distance to crater rim (m)')
plt.ylabel('density')
plt.title('Scatterplot of Features')
plt.show()

    