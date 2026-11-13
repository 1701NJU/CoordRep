"""Crop standalone vector C/D PDFs after render_figure4.py; no model inference."""
from pathlib import Path
import argparse,json
from pypdf import PdfReader,PdfWriter

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--input',type=Path,default=Path(__file__).resolve().parent/'rendered')
args=p.parse_args()
root=args.input
a=json.loads((root/'FIGURE4_UPDATE_AUDIT.json').read_text())
w,h=a['source_image_size'];scale=504/w
for name,b in a['pdf_crop_bounds_pixels'].items():
    page=PdfReader(root/'qa'/'CD_vector_canvas.pdf').pages[0]
    x,y,bw,bh=b;box=[x*scale,(h-y-bh)*scale,(x+bw)*scale,(h-y)*scale]
    for attr in ['mediabox','cropbox']:
        getattr(page,attr).lower_left=(box[0],box[1]);getattr(page,attr).upper_right=(box[2],box[3])
    writer=PdfWriter();writer.add_page(page)
    with (root/'output'/'pdf'/f'Figure4{name}_updated.pdf').open('wb') as f:writer.write(f)
print('PASS: vector C, D and CD PDFs exported')
