from PIL import Image
import os

files = [f for f in os.listdir('exported_data') if f.startswith('viirs_frp_') and f.endswith('.png')]
print('PNG file properties:')
for f in sorted(files):
    img = Image.open(f"exported_data/{f}")
    print(f'{f}: mode={img.mode}, has_transparency={img.info.get("transparency", False)}')
    img.close()
