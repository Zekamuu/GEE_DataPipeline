import ee
import geemap
import leafmap
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import os
import rasterio
from rasterio.plot import reshape_as_image

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

start_date = '2025-08-01'
end_date = '2025-08-15'

def generate_html_map_from_tif(export_data, standard_values, output_html="exported_map.html"):
    import folium
    m = folium.Map(location=[28.6139, 77.2090], zoom_start=8, tiles='OpenStreetMap')

    for key in export_data:
        
        print(f"Accepted standard value for {key}: {standard_values.get(key, 'N/A')}")
        tif_path = os.path.join("exported_data", start_date, f"{key}.tif")
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

                std_val = standard_values.get(key, None)
                if std_val is None:
                    print(f"No standard value for {key}, skipping color assignment.")
                    color_img = np.zeros((img.shape[0], img.shape[1], 4), dtype=np.uint8)
                else:
                    color_img = np.zeros((img.shape[0], img.shape[1], 4), dtype=np.uint8)  # RGBA

                    # Masks for buckets
                    mask_transparent = img < std_val
                    mask_blue = (img >= std_val) & (img < 1.5 * std_val)
                    mask_green = (img >= 1.5 * std_val) & (img < 2 * std_val)
                    mask_red = img >= 2 * std_val

                    # Assign colors
                    color_img[mask_transparent] = [0, 0, 0, 0]      # Transparent
                    color_img[mask_blue] = [0, 0, 255, 180]         # Blue
                    color_img[mask_green] = [0, 255, 0, 180]        # Green
                    color_img[mask_red] = [255, 0, 0, 180]          # Red

                from PIL import Image
                temp_png = tif_path.replace('.tif', '.png')
                Image.fromarray(color_img).save(temp_png)

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
    if os.path.exists(powerplants_geojson):
        folium.GeoJson(
            powerplants_geojson,
            name="Power Plants",
            tooltip=folium.GeoJsonTooltip(fields=["name"], aliases=["Power Plant Name"])
        ).add_to(m)
    else:
        print(f"GeoJSON file not found: {powerplants_geojson}")

    folium.LayerControl().add_to(m)
    m.save(output_html)
    print(f"Interactive map saved to {output_html}")

def generate_tif():
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
        mean_file_path = os.path.join(
            band_dir, f"{key}.tif"
        )
        geemap.ee_export_image(mean_img, filename=mean_file_path, scale=1000, region=study_area_geometry, file_per_band=False)

    print(f"Data exported to {output_dir}. You can reload from these files without invoking EE.")

# generate_tif()

# key can be 'aerosol', 'no2', 'so2', 'co'
# Define accepted standard values (example values, update as per latest guidelines)
standard_values = {
    "aerosol": 0.5,   # Absorbing Aerosol Index (unitless, threshold for concern)
    "no2": 0.04,      # NO2 column number density (mol/m^2, WHO annual mean guideline ~0.04)
    "so2": 0.02,      # SO2 column number density (mol/m^2, WHO guideline ~0.02)
    "co": 0.05        # CO column number density (mol/m^2, WHO guideline ~0.05)
}

standard_values = {
    "aerosol": 0.0,   # Absorbing Aerosol Index (unitless, threshold for concern)
    "no2": 0.00,      # NO2 column number density (mol/m^2, WHO annual mean guideline ~0.04)
    "so2": 0.00,      # SO2 column number density (mol/m^2, WHO guideline ~0.02)
    "co": 0.00        # CO column number density (mol/m^2, WHO guideline ~0.05)
}

generate_html_map_from_tif(['aerosol', 'no2', 'so2', 'co'], standard_values)
