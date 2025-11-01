import ee
import geemap
import leafmap
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import os
import rasterio
from rasterio.plot import reshape_as_image
from dateutil.relativedelta import relativedelta
import glob
from urllib.parse import urlencode

import matplotlib.pyplot as plt
ee.Initialize(project='valued-proton-423000-v7')
india_states = ee.FeatureCollection('WM/geoLab/geoBoundaries/600/ADM1')
india_boundary = india_states.filter(ee.Filter.eq('shapeGroup', 'IND'))
# study_area = india_boundary.filter(
#     ee.Filter.inList('shapeName', ['Delhi', 'Noida', 'Gurugram', 'Ghaziabad', 'Greater Noida'])
# )
# study_area_geometry = study_area.geometry()

# Define bounding box coordinates (min_lon, min_lat, max_lon, max_lat) covering Delhi, Noida, Gurugram, Ghaziabad, Greater Noida
min_lon, min_lat = 76.85, 28.2
max_lon, max_lat = 77.65, 29.0

# Create a rectangle geometry for the study area
study_area_geometry = ee.Geometry.Rectangle([min_lon, min_lat, max_lon, max_lat])



def generate_html_map_from_tif(start_date, end_date, export_data):
    import folium
    from urllib.parse import parse_qs
    import sys

    # Get layers from querystring if running as a script with querystring argument
    layers_to_show = set(export_data + ['powerplants'])
    if len(sys.argv) > 1 and "layers=" in sys.argv[1]:
        qs = parse_qs(sys.argv[1].lstrip('?'))
        if "layers" in qs:
            layers_to_show = set(qs["layers"][0].split(','))

    m = folium.Map(location=[28.6139, 77.2090], zoom_start=10, tiles='OpenStreetMap')
    output_html = os.path.join("exported_data", start_date, f"exported_map.html")
    for key in export_data:
        if key not in layers_to_show:
            continue
        tif_path = os.path.join("exported_data", start_date, f"{key}_{start_date}_{end_date}.tif")
        if os.path.exists(tif_path):
            with rasterio.open(tif_path) as src:
                import json
                from shapely.geometry import shape, mapping
                from rasterio.mask import mask

                geojson = study_area_geometry.getInfo()
                shapes = [shape(geojson)]

                apply_mask = False

                if apply_mask:
                    masked_img, masked_transform = mask(src, shapes, crop=True, nodata=np.nan)
                    img = masked_img[0]
                else:
                    img = src.read(1)
                    masked_transform = src.transform

                min_val = np.nanmin(img)
                max_val = np.nanmax(img)
                norm_img = (img - min_val) / (max_val - min_val + 1e-6)
                import matplotlib
                cmap = matplotlib.cm.get_cmap('jet')
                rgba_img = cmap(norm_img)
                rgb_img = (rgba_img[:, :, :3] * 255).astype(np.uint8)

                from PIL import Image
                temp_png = tif_path.replace('.tif', '.png')
                Image.fromarray(rgb_img).save(temp_png)

                from rasterio.transform import array_bounds
                height, width = img.shape
                bounds = array_bounds(height, width, masked_transform)

                folium.raster_layers.ImageOverlay(
                    image=temp_png,
                    bounds=[[bounds[1], bounds[0]], [bounds[3], bounds[2]]],
                    name=f"{key.upper()} Mean",
                    opacity=0.7,
                    interactive=True,
                    cross_origin=False,
                    zindex=1,
                ).add_to(m)
        else:
            print(f"File not found: {tif_path}")

    powerplants_geojson = os.path.join("AminitiesGeoJSON", "PowerPlants.geojson")
    if "powerplants" in layers_to_show and os.path.exists(powerplants_geojson):
        folium.GeoJson(
            powerplants_geojson,
            name="Power Plants",
            tooltip=folium.GeoJsonTooltip(fields=["name"], aliases=["Power Plant Name"])
        ).add_to(m)
    elif "powerplants" in layers_to_show:
        print(f"GeoJSON file not found: {powerplants_geojson}")

    folium.LayerControl().add_to(m)
    m.save(output_html)
    print(f"Interactive map saved to {output_html}")

def generate_tif(start_date, end_date):
    aerosol_collection = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_AER_AI') \
                .filterDate(start_date, end_date) \
                .filterBounds(study_area_geometry) \
                .select('absorbing_aerosol_index')

    aerosol_collection_size = aerosol_collection.size().getInfo()
    if aerosol_collection_size <= 0:
        print("No aerosol data available for the specified date range and location.")

    no2_collection = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2') \
                .filterDate(start_date, end_date) \
                .filterBounds(study_area_geometry) \
                .select('NO2_column_number_density')

    no2_collection_size = no2_collection.size().getInfo()
    if no2_collection_size <= 0:
        print("No NO2 data available for the specified date range and location.")

    so2_collection = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_SO2') \
                .filterDate(start_date, end_date) \
                .filterBounds(study_area_geometry) \
                .select('SO2_column_number_density')

    so2_collection_size = so2_collection.size().getInfo()
    if so2_collection_size <= 0:
        print("No SO2 data available for the specified date range and location.")

    co_collection = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_CO') \
                .filterDate(start_date, end_date) \
                .filterBounds(study_area_geometry) \
                .select('CO_column_number_density')

    co_collection_size = co_collection.size().getInfo()
    if co_collection_size <= 0:
        print("No CO data available for the specified date range and location.")

    export_data = {
        "aerosol": aerosol_collection,
        "no2": no2_collection,
        "so2": so2_collection,
        "co": co_collection
    }
    
    output_dir = "exported_data"
    os.makedirs(output_dir, exist_ok=True)
    for key, collection in export_data.items():
        band_name = collection.first().bandNames().get(0).getInfo()
        band_dir = os.path.join(output_dir, start_date)
        os.makedirs(band_dir, exist_ok=True)
        # Export mean image clipped to study area
        mean_img = collection.mean().select(band_name).clip(study_area_geometry)
        mean_file_path = os.path.join("exported_data", start_date, f"{key}_{start_date}_{end_date}.tif")
        # mean_file_path = os.path.join(
        #     band_dir, f"{key}.tif"
        # )
        geemap.ee_export_image(mean_img, filename=mean_file_path, scale=1000, region=study_area_geometry, file_per_band=False)

    print(f"Data exported to {output_dir}. You can reload from these files without invoking EE.")

year = datetime.now().year
start = datetime(year, 1, 1)
end = datetime(year, 12, 31)

current_start = start
while current_start < end:
    current_end = current_start + relativedelta(months=1) - timedelta(days=1)
    if current_end > end:
        current_end = end
    start_date = current_start.strftime('%Y-%m-%d')
    end_date = current_end.strftime('%Y-%m-%d')
    # generate_tif(start_date, end_date)
    generate_html_map_from_tif(start_date, end_date, ['aerosol', 'no2', 'so2', 'co'])
    current_start = current_start + relativedelta(months=1)


def generate_wrapper_html():
    output_dir = "exported_data"
    layer_options = ['aerosol', 'no2', 'so2', 'co', 'powerplants']
    subfolders = sorted([f for f in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, f))])

    wrapper_html_path = os.path.join(output_dir, "wrapper.html")
    with open(wrapper_html_path, "w") as f:
        f.write("""
<!DOCTYPE html>
<html>
<head>
    <title>Exported Maps Wrapper</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .map-list { margin-bottom: 20px; }
        .map-item { margin-bottom: 10px; }
        .iframe-container { width: 100%; height: 600px; border: 1px solid #ccc; }
        .checkbox-group { margin-bottom: 20px; }
        .checkbox-group label { margin-right: 15px; }
        .slider-container { margin-bottom: 20px; }
    </style>
</head>
<body>
    <h2>Exported Maps Viewer</h2>
    <div class="checkbox-group">
        <label>Select layers to show:</label>
""")
        for layer in layer_options:
            display_name = layer.upper() if layer != "powerplants" else "Power Plants"
            f.write(f'        <label><input type="checkbox" class="layer-checkbox" value="{layer}" checked> {display_name}</label>\n')
        f.write("""
    </div>
    <div class="slider-container">
        <label for="date-slider">Select map date:</label>
        <input type="range" id="date-slider" min="0" max="{max_idx}" value="0" step="1">
        <span id="date-label"></span>
    </div>
    <div class="iframe-container">
        <iframe id="map-frame" src="" width="100%" height="600"></iframe>
    </div>
    <script>
        var subfolders = {subfolders};
        var dateSlider = document.getElementById('date-slider');
        var dateLabel = document.getElementById('date-label');
        dateSlider.max = subfolders.length - 1;
        function updateDateLabel() {
            var idx = parseInt(dateSlider.value);
            dateLabel.textContent = subfolders[idx];
        }
        function updateIframe() {
            var idx = parseInt(dateSlider.value);
            var selectedDate = subfolders[idx];
            var checkboxes = document.getElementsByClassName('layer-checkbox');
            var selectedLayers = [];
            for (var i = 0; i < checkboxes.length; i++) {
                if (checkboxes[i].checked) {
                    selectedLayers.push(checkboxes[i].value);
                }
            }
            var params = selectedLayers.length ? "?layers=" + selectedLayers.join(",") : "";
            var iframe = document.getElementById('map-frame');
            
            iframe.src = "./" + selectedDate + "/exported_map.html" + params;
            alert(iframe.src);
            updateDateLabel();
        }
        dateSlider.addEventListener('input', updateIframe);
        var checkboxes = document.getElementsByClassName('layer-checkbox');
        for (var i = 0; i < checkboxes.length; i++) {
            checkboxes[i].addEventListener('change', updateIframe);
        }
        window.onload = function() {
            updateDateLabel();
            updateIframe();
        };
    </script>
</body>
</html>
""".replace("{max_idx}", str(len(subfolders)-1)).replace("{subfolders}", str(subfolders)))
    print(f"Wrapper HTML generated at {wrapper_html_path}")

generate_wrapper_html()