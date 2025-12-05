"""
Create Figure 2: Mechanism Insights (4-panel combined figure)
=============================================================
Panel A: Liver Variance vs Leverage Scatter
Panel B: Brain Leverage Bar Plot
Panel C: Spatial Regularization Low Coverage
Panel D: Alpha Sensitivity with spatial visualization

Use high-resolution PDF-based approach for publication quality.
"""

from PIL import Image
import numpy as np

# Load high-resolution images
img_a = Image.open('paper/figures/supplementary_rare_vs_abundant_final.png')
img_b = Image.open('paper/figures/supplementary_brain_leverage.png')
img_c = Image.open('paper/figures/deep_dive_low_coverage.png')
img_d = Image.open('paper/figures/supplementary_alpha_sensitivity.png')

print(f"Original sizes: A={img_a.size}, B={img_b.size}, C={img_c.size}, D={img_d.size}")

# Target width for each panel (half of final width)
target_width = 2400  # pixels

# Resize images maintaining aspect ratio
def resize_to_width(img, target_w):
    w, h = img.size
    new_h = int(h * target_w / w)
    return img.resize((target_w, new_h), Image.LANCZOS)

img_a_r = resize_to_width(img_a, target_width)
img_b_r = resize_to_width(img_b, target_width)
img_c_r = resize_to_width(img_c, target_width)
img_d_r = resize_to_width(img_d, target_width)

print(f"Resized: A={img_a_r.size}, B={img_b_r.size}, C={img_c_r.size}, D={img_d_r.size}")

# Calculate heights for each row
row1_height = max(img_a_r.size[1], img_b_r.size[1])
row2_height = max(img_c_r.size[1], img_d_r.size[1])

# Total canvas size
total_width = target_width * 2 + 20  # 20px gap
total_height = row1_height + row2_height + 30  # 30px gap between rows

# Create white canvas
canvas = Image.new('RGB', (total_width, total_height), 'white')

# Paste images
# Row 1
canvas.paste(img_a_r, (0, 0))
canvas.paste(img_b_r, (target_width + 20, 0))

# Row 2
y_offset = row1_height + 30
canvas.paste(img_c_r, (0, y_offset))
canvas.paste(img_d_r, (target_width + 20, y_offset))

# Add panel labels using PIL
from PIL import ImageDraw, ImageFont

draw = ImageDraw.Draw(canvas)
try:
    font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 72)
except:
    font = ImageFont.load_default()

# Panel labels with white background
label_positions = [
    ('a', 30, 30),
    ('b', target_width + 50, 30),
    ('c', 30, y_offset + 30),
    ('d', target_width + 50, y_offset + 30),
]

for label, x, y in label_positions:
    # White background rectangle
    bbox = draw.textbbox((x, y), label, font=font)
    draw.rectangle([bbox[0]-10, bbox[1]-5, bbox[2]+10, bbox[3]+5], fill='white')
    draw.text((x, y), label, fill='black', font=font)

# Save
canvas.save('paper/figures/figure2_mechanism.png', dpi=(300, 300))
canvas.save('paper/figures/figure2_mechanism.pdf', dpi=(300, 300))

print(f"✅ Saved: paper/figures/figure2_mechanism.png/pdf ({total_width}x{total_height}px)")
