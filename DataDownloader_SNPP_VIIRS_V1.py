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
import folium
from urllib.parse import parse_qs
import sys
import json
from shapely.geometry import shape, mapping
from rasterio.mask import mask
import matplotlib
from PIL import Image
from rasterio.transform import array_bounds
import json
from shapely.geometry import shape, mapping
from rasterio.mask import mask
import matplotlib.pyplot as plt

# Initialize Earth Engine
ee.Initialize(project='valued-proton-423000-v7')

# Define Region of Interest (ROI) - Punjab, India
# Load the L1 administrative boundaries and filter for Punjab, India
admin_boundaries = ee.FeatureCollection('FAO/GAUL/2015/level1')
punjab = admin_boundaries.filter(ee.Filter.eq('ADM0_NAME', 'India')).filter(ee.Filter.eq('ADM1_NAME', 'Punjab'))

# Get Punjab geometry for clipping
punjab_geometry = punjab.geometry()

def generate_html_map_from_tif(start_date, end_date, export_data):
    """
    Generate an interactive HTML map from exported TIF files
    """
    # Get layers from querystring if running as a script with querystring argument
    layers_to_show = set(export_data + ['powerplants'])
    if len(sys.argv) > 1 and "layers=" in sys.argv[1]:
        qs = parse_qs(sys.argv[1].lstrip('?'))
        if "layers" in qs:
            layers_to_show = set(qs["layers"][0].split(','))

    # Center map on Punjab
    m = folium.Map(location=[30.8143, 75.1712], zoom_start=8, tiles='OpenStreetMap')
    output_html = os.path.join("exported_data", start_date, f"exported_map.html")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_html), exist_ok=True)
    
    for key in export_data:
        if key not in layers_to_show:
            continue
        tif_path = os.path.join("exported_data", f"{key}_{start_date}.tif")
        if os.path.exists(tif_path):
            with rasterio.open(tif_path) as src:
                geojson = punjab_geometry.getInfo()
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
                
                # Create normalized image, but keep zero values as zero for transparency
                norm_img = np.zeros_like(img)
                if max_val > min_val:
                    # Only normalize non-zero values
                    mask = img > 0
                    norm_img[mask] = (img[mask] - min_val) / (max_val - min_val)
                
                # Use dark red colormap for all fires
                cmap = matplotlib.cm.get_cmap('Reds')
                rgba_img = cmap(norm_img)
                
                # Force all fire pixels to be dark red
                rgba_img[img > 0, 0] = 0.8  # Red channel
                rgba_img[img > 0, 1] = 0.0  # Green channel  
                rgba_img[img > 0, 2] = 0.0  # Blue channel
                rgba_img[img > 0, 3] = 1.0  # Alpha channel (fully opaque)
                
                # Make background transparent (where img == 0)
                rgba_img[img == 0, 3] = 0  # Set alpha to 0 for zero values
                
                # Convert to RGBA image
                rgba_img = (rgba_img * 255).astype(np.uint8)

                temp_png = tif_path.replace('.tif', '.png')
                Image.fromarray(rgba_img).save(temp_png)

                # Calculate bounds for the image overlay
                height, width = img.shape
                bounds = array_bounds(height, width, masked_transform)

                # Add the fire data as an image overlay with transparent background
                folium.raster_layers.ImageOverlay(
                    image=temp_png,
                    bounds=[[bounds[1], bounds[0]], [bounds[3], bounds[2]]],
                    name=f"{key.upper()} (Fire Radiative Power)",
                    opacity=0.8,
                    interactive=True,
                    cross_origin=False,
                    zindex=1,
                ).add_to(m)
        else:
            print(f"File not found: {tif_path}")

    # Add Punjab boundary for context
    folium.GeoJson(
        punjab_geometry.getInfo(),
        name="Punjab Boundary",
        style_function=lambda x: {'color': 'black', 'weight': 2, 'fillOpacity': 0}
    ).add_to(m)

    folium.LayerControl().add_to(m)
    m.save(output_html)
    print(f"Interactive map saved to {output_html}")

def generate_tif(start_date, end_date):
    """
    Generate TIF files for VIIRS Active Fire data
    """
    print(f"Processing VIIRS Active Fire data for {start_date} to {end_date}")
    
    # Load the VIIRS Active Fire Dataset
    # This is the specific NRT (Near Real-Time) dataset
    viirs_fires = ee.ImageCollection('NASA/LANCE/SNPP_VIIRS/C2') \
                   .filterBounds(punjab_geometry) \
                   .filterDate(start_date, end_date)
    
    # Check if data is available
    collection_size = viirs_fires.size().getInfo()
    if collection_size <= 0:
        print(f"No VIIRS fire data available for the specified date range {start_date} to {end_date}")
        return
    
    print(f"Found {collection_size} VIIRS fire images for the specified period")
    
    # Select the Fire Radiative Power band
    frp_collection = viirs_fires.select('frp')
    
    # Create a max composite to show the most intense fire at each pixel
    max_frp = frp_collection.max()
    
    # Export data
    export_data = {
        "viirs_frp": max_frp
    }
    
    output_dir = "exported_data"
    os.makedirs(output_dir, exist_ok=True)
    
    for key, image in export_data.items():
        band_dir = os.path.join(output_dir, start_date)
        os.makedirs(band_dir, exist_ok=True)
        
        # Export max FRP image clipped to Punjab
        clipped_img = image.clip(punjab_geometry)
        mean_file_path = os.path.join("exported_data", f"{key}_{start_date}.tif")
        
        # Export with appropriate scale for fire detection (375m native resolution)
        geemap.ee_export_image(
            clipped_img, 
            filename=mean_file_path, 
            scale=375,  # Native resolution of VIIRS
            region=punjab_geometry, 
            file_per_band=False
        )
        
        # Generate visualization
        with rasterio.open(mean_file_path) as src:
            geojson = punjab_geometry.getInfo()
            shapes = [shape(geojson)]

            apply_mask = False

            if apply_mask:
                masked_img, masked_transform = mask(src, shapes, crop=True, nodata=np.nan)
                img = masked_img[0]
            else:
                img = src.read(1)
                masked_transform = src.transform

            # Handle NaN values
            img = np.where(np.isnan(img), 0, img)
            
            min_val = np.nanmin(img)
            max_val = np.nanmax(img)
            
            if max_val > min_val:
                norm_img = (img - min_val) / (max_val - min_val)
            else:
                norm_img = img
            
            # Use fire-specific colormap (red to dark red)
            cmap = matplotlib.cm.get_cmap('Reds')
            rgba_img = cmap(norm_img)
            rgb_img = (rgba_img[:, :, :3] * 255).astype(np.uint8)

            temp_png = mean_file_path.replace('.tif', '.png')
            Image.fromarray(rgb_img).save(temp_png)
            
            print(f"Exported {key} data: {mean_file_path}")
            print(f"FRP range: {min_val:.2f} - {max_val:.2f} MW")

    print(f"VIIRS fire data exported to {output_dir}")

def download_monthly_data(year, start_month=10, end_month=11):
    """
    Download VIIRS fire data for stubble burning season (October-November)
    """
    print(f"Downloading VIIRS fire data for {year} stubble burning season (months {start_month}-{end_month})")
    
    for month in range(start_month, end_month + 1):
        # Calculate start and end dates for the month
        start_date = datetime(year, month, 1)
        if month == 12:
            end_date = datetime(year + 1, 1, 1) - timedelta(days=1)
        else:
            end_date = datetime(year, month + 1, 1) - timedelta(days=1)
        
        start_date_str = start_date.strftime('%Y-%m-%d')
        end_date_str = end_date.strftime('%Y-%m-%d')
        
        print(f"\nProcessing month {month} ({start_date_str} to {end_date_str})")
        generate_tif(start_date_str, end_date_str)

def download_specific_period(start_date, end_date):
    """
    Download VIIRS fire data for a specific date range
    """
    print(f"Downloading VIIRS fire data for {start_date} to {end_date}")
    generate_tif(start_date, end_date)

if __name__ == "__main__":
    # Regenerate PNG files with proper transparency
    import glob
    
    print("Regenerating PNG files with proper transparency...")
    tif_files = glob.glob("exported_data/viirs_frp_*.tif")
    
    for tif_file in sorted(tif_files):
        print(f"Processing {tif_file}")
        
        with rasterio.open(tif_file) as src:
            img = src.read(1)
            masked_transform = src.transform

            # Handle NaN values
            img = np.where(np.isnan(img), 0, img)
            
            min_val = np.nanmin(img)
            max_val = np.nanmax(img)
            
            # Create normalized image, but keep zero values as zero for transparency
            norm_img = np.zeros_like(img)
            if max_val > min_val:
                # Only normalize non-zero values
                mask = img > 0
                norm_img[mask] = (img[mask] - min_val) / (max_val - min_val)
            
            # Use fire-specific colormap (red to dark red) with transparency
            cmap = matplotlib.cm.get_cmap('Reds')
            rgba_img = cmap(norm_img)
            
            # Force all fire pixels to be dark red
            rgba_img[img > 0, 0] = 0.8  # Red channel
            rgba_img[img > 0, 1] = 0.0  # Green channel  
            rgba_img[img > 0, 2] = 0.0  # Blue channel
            rgba_img[img > 0, 3] = 1.0  # Alpha channel (fully opaque)
            
            # Make background transparent (where img == 0)
            rgba_img[img == 0, 3] = 0  # Set alpha to 0 for zero values
            
            # Convert to RGBA image
            rgba_img = (rgba_img * 255).astype(np.uint8)

            png_path = tif_file.replace('.tif', '.png')
            Image.fromarray(rgba_img).save(png_path)
            print(f"Regenerated {png_path} with RGBA transparency")
    
    print("All PNG files regenerated with proper transparency!")
