import os
import rasterio
import numpy as np
import matplotlib
from PIL import Image
from rasterio.transform import array_bounds
from shapely.geometry import shape

def regenerate_png_with_transparency(tif_path):
    """Regenerate PNG with proper transparency"""
    print(f"Regenerating {tif_path}")
    
    with rasterio.open(tif_path) as src:
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

        png_path = tif_path.replace('.tif', '.png')
        Image.fromarray(rgba_img).save(png_path)
        print(f"Saved {png_path} with RGBA transparency")

# Find all VIIRS TIF files and regenerate their PNGs
tif_files = [f for f in os.listdir('exported_data') if f.startswith('viirs_frp_') and f.endswith('.tif')]

print(f"Found {len(tif_files)} VIIRS TIF files to regenerate:")
for tif_file in sorted(tif_files):
    tif_path = os.path.join('exported_data', tif_file)
    regenerate_png_with_transparency(tif_path)

print("All PNG files regenerated with proper transparency!")
