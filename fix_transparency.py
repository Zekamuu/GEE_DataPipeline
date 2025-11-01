import os
import rasterio
import numpy as np
import matplotlib
from PIL import Image

def regenerate_png_with_transparency(tif_path):
    """Regenerate PNG with proper transparency"""
    print(f"Processing {tif_path}")
    
    with rasterio.open(tif_path) as src:
        img = src.read(1)

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

        png_path = tif_path.replace('.tif', '.png')
        Image.fromarray(rgba_img).save(png_path)
        print(f"Regenerated {png_path} with RGBA transparency")

# Find all VIIRS TIF files and regenerate their PNGs
tif_files = [f for f in os.listdir('exported_data') if f.startswith('viirs_frp_') and f.endswith('.tif')]

print(f"Found {len(tif_files)} VIIRS TIF files to regenerate:")
for tif_file in sorted(tif_files):
    tif_path = os.path.join('exported_data', tif_file)
    regenerate_png_with_transparency(tif_path)
    # Export the TIF as a CSV (array rows -> CSV rows, columns -> CSV columns)
    with rasterio.open(tif_path) as src:
        arr = src.read(1).astype(float)

        # Respect nodata as NaN so it won't be mistaken for a real value
        if src.nodata is not None:
            arr[arr == src.nodata] = np.nan

        height, width = arr.shape

        # Get row/col indices for each pixel, then pixel center coordinates in source CRS
        rows, cols = np.meshgrid(np.arange(height), np.arange(width), indexing='ij')
        rows_flat = rows.ravel()
        cols_flat = cols.ravel()
        xs, ys = rasterio.transform.xy(src.transform, rows_flat, cols_flat, offset="center")

        # Transform coordinates to EPSG:4326 (lon, lat). If transform fails, fall back to source coords.
        try:
            lon, lat = rasterio.warp.transform(src.crs, 'EPSG:4326', xs, ys)
        except Exception:
            lon, lat = xs, ys

        # Prepare output: lat, lon, mean_value (for a single-band raster the pixel value is used as "mean")
        values = arr.ravel()
        out_arr = np.column_stack((lat, lon, values))

        csv_path = tif_path.replace('.tif', '.csv')
        header = 'lat,lon,mean'
        np.savetxt(csv_path, out_arr, delimiter=',', fmt='%.6f', header=header, comments='')
        print(f"Exported {csv_path} with lat, lon and mean columns")

print("All PNG files regenerated with proper transparency!")

