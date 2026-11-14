"""Update Figure 4 C/D; preserve the submitted A/B image without redrawing it."""
from __future__ import annotations
import argparse, csv, hashlib, io, json, sys
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image, ImageChops
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

sys.stdout.reconfigure(encoding='utf-8')
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=HERE/'rendered')
parser.add_argument('--source',type=Path,default=HERE/'Figure4CD_source_data.csv')
parser.add_argument('--ab-image',type=Path,default=HERE/'Figure4_AB_original.png')
parser.add_argument('--font-dir',type=Path,default=None,help='Directory containing licensed Arial TTF files')
args=parser.parse_args()
OUT=args.output.resolve()
OUT.mkdir(parents=True,exist_ok=True)
PDF=OUT/'output'/'pdf'
PDF.mkdir(parents=True,exist_ok=True)
QA=OUT/'qa'
QA.mkdir(exist_ok=True)
W,H=2704,2034
KEEP_AB=1000
NATIVE_DPI=W/7.0
HEIGHT=(H+1e-6)/NATIVE_DPI
NS='http://www.w3.org/2000/svg'

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def json_write(p,obj): p.write_text(json.dumps(obj,indent=2,ensure_ascii=True),encoding='utf-8')
def panel_data(panel):
    with args.source.open(encoding='utf-8-sig',newline='') as f:
        rows=[r for r in csv.DictReader(f) if r['panel']==panel]
    models=(['gine_wide','coordrep_relation_masked','coordrep_relation_full','schnet_3d','visnet_3d']
            if panel=='C' else ['gine_wide','hybrid_mask_shape','hybrid','schnet_3d','visnet_3d'])
    assert [r['model'] for r in rows]==models
    for r in rows:
        for key in ['value','ci95_low','ci95_high','direction','direction_ci95_low','direction_ci95_high']:
            r[key]=float(r[key])
        for key in ['n_pairs','n_graph_clusters','bootstrap_replicates']:r[key]=int(r[key])
        assert r['n_pairs']==(34 if panel=='C' else 44)
        assert r['n_graph_clusters']==(34 if panel=='C' else 31)
        assert r['bootstrap_replicates']==20000
        assert 0<=r['ci95_low']<=r['value']<=r['ci95_high']
    return rows

COLORS=['#929AA2','#58A9B5','#155B86','#EAA21A','#76508F']
LABELS_C=['GINE\n(2D graph)','CoordRep hybrid\nrelations hidden','CoordRep hybrid\nfull record','SchNet\n(3D)','ViSNet\n(equivariant 3D)']
LABELS_D=['GINE\n(2D graph)','CoordRep hybrid\nshape/CShM hidden','CoordRep hybrid\nfull record','SchNet\n(3D)','ViSNet\n(equivariant 3D)']

def draw_panel(fig,rows,panel):
    left,right=(280,1312) if panel=='C' else (1701,2673)
    top,bottom=1134,1903
    ax=fig.add_axes([left/W,(H-bottom)/H,(right-left)/W,(bottom-top)/H])
    ax.set_gid('panel_'+panel)
    ax.axhspan(1.55,3.45,color='#EAF4F7',zorder=0)
    ax.axhline(1.5,color='#D5DADF',lw=.65,zorder=0)
    ax.set_xlim(0,max(r['ci95_high'] for r in rows)*1.18)
    ax.set_ylim(-.65,4.75)
    ax.grid(axis='x',color='#DDE2E6',lw=.55,zorder=0)
    artists=[]
    for y,r,c in zip([4,3,2,1,0],rows,COLORS):
        line=ax.plot([r['ci95_low'],r['ci95_high']],[y,y],color=c,
            lw=2 if y==2 else 1.45,solid_capstyle='round',zorder=2)[0]
        line.set_gid(f'{panel}_{r["model"]}_CI95')
        point=ax.scatter([r['value']],[y],s=34 if y==2 else 27,
            facecolor=c if y==2 else 'white',edgecolor=c,linewidth=1.25,zorder=3)
        point.set_gid(f'{panel}_{r["model"]}_point')
        text=ax.text(r['value'],y+.27,f'{r["value"]:.3f}',ha='center',va='bottom',fontsize=5.7)
        text.set_gid(f'{panel}_{r["model"]}_value')
        assert list(line.get_xdata())==[r['ci95_low'],r['ci95_high']]
        assert float(point.get_offsets()[0][0])==r['value']
        artists.append((r,line,point,text))
    ax.set_yticks([4,3,2,1,0],LABELS_C if panel=='C' else LABELS_D)
    ax.tick_params(axis='y',length=0,pad=5)
    ax.set_xlabel('Dipole pair-difference MAE (D)' if panel=='C' else 'Gap pair-difference MAE (eV)')
    ax.spines[['top','right','left']].set_visible(False)
    title='Donor relations recover cis/trans ordering (n = 34)' if panel=='C' else 'Shape fields recover gap ordering (n = 44)'
    ax.text(-.145 if panel=='C' else -.14,1.035,f'({panel})',transform=ax.transAxes,
        ha='left',va='bottom',fontsize=9.1,fontweight='bold',color='#17212A',clip_on=False)
    ax.text(0,1.035,title,transform=ax.transAxes,ha='left',va='bottom',
        fontsize=7.9,fontweight='bold',color='#17212A',clip_on=False)
    reduction=1-rows[2]['value']/rows[1]['value']
    x=rows[1]['value']*1.06
    ax.annotate('',xy=(x,2.04),xytext=(x,2.96),
        arrowprops=dict(arrowstyle='-|>',color=COLORS[2],lw=1.1))
    ax.text(x*1.015,2.5,f'{100*reduction:.1f}% lower MAE\n'
        f'ordering {100*rows[1]["direction"]:.1f}% \u2192 {100*rows[2]["direction"]:.1f}%',
        ha='left',va='center',fontsize=5.65,color=COLORS[2],fontweight='bold')
    return ax,artists

def svg_crop(source,dest,bounds,keep):
    tree=ET.parse(source); r=tree.getroot()
    for parent in list(r.iter()):
        for node in list(parent):
            gid=node.get('id','')
            if gid.startswith('panel_') and gid not in keep: parent.remove(node)
    scale=504/W
    x,y,w,h=bounds
    r.set('viewBox',' '.join(f'{v*scale:.9f}' for v in [x,y,w,h]))
    r.set('width',f'{w*scale:.9f}pt');r.set('height',f'{h*scale:.9f}pt')
    ET.register_namespace('',NS); ET.register_namespace('xlink','http://www.w3.org/1999/xlink')
    tree.write(dest,encoding='utf-8',xml_declaration=True)

def main():
    font_dir=args.font_dir or Path('C:/Windows/Fonts')
    font_paths=[font_dir/n for n in ['arial.ttf','arialbd.ttf','ariali.ttf','arialbi.ttf'] if (font_dir/n).is_file()]
    for p in font_paths: font_manager.fontManager.addfont(str(p))
    mpl.rcParams.update({'font.family':'Arial','font.size':6.7,'axes.labelsize':6.7,
        'xtick.labelsize':6,'ytick.labelsize':6,'axes.linewidth':.65,
        'xtick.major.width':.65,'ytick.major.width':.65,'svg.fonttype':'none',
        'pdf.fonttype':42,'ps.fonttype':42,'text.color':'#202326',
        'axes.labelcolor':'#202326','xtick.color':'#202326','ytick.color':'#202326'})
    resolved_font=font_manager.findfont(font_manager.FontProperties(family='Arial'),fallback_to_default=False)
    assert font_manager.FontProperties(fname=resolved_font).get_name()=='Arial'
    ab=Image.open(args.ab_image).convert('RGB')
    assert ab.size==(W,KEEP_AB)
    c,d=panel_data('C'),panel_data('D')
    fig=plt.figure(figsize=(7,HEIGHT),dpi=NATIVE_DPI,facecolor='white')
    ax_ab=fig.add_axes([0,(H-KEEP_AB)/H,1,KEEP_AB/H],zorder=0)
    ax_ab.set_gid('panel_AB_original')
    ax_ab.imshow(np.asarray(ab),interpolation='none',aspect='auto',extent=(0,W,KEEP_AB,0))
    ax_ab.set_xlim(0,W);ax_ab.set_ylim(KEEP_AB,0);ax_ab.axis('off')
    ax_c,ca=draw_panel(fig,c,'C');ax_d,da=draw_panel(fig,d,'D')
    fig.canvas.draw()
    fonts=sorted({t.get_fontfamily()[0] for t in fig.findobj(mpl.text.Text) if t.get_text()})
    assert fonts==['Arial'],fonts
    stem=OUT/'Figure4_updated_20261114'
    fig.savefig(stem.with_suffix('.svg'),facecolor='white')
    fig.savefig(PDF/'Figure4_updated_20261114.pdf',facecolor='white')
    buf=io.BytesIO();fig.savefig(buf,format='png',dpi=NATIVE_DPI,facecolor='white');buf.seek(0)
    native=Image.open(buf).convert('RGB')
    assert native.size==(W,H),native.size
    native.paste(ab,(0,0))
    native.save(OUT/'Figure4_updated_native.png',dpi=(NATIVE_DPI,NATIVE_DPI))
    assert ImageChops.difference(ab,native.crop((0,0,W,KEEP_AB))).getbbox() is None
    buf=io.BytesIO();fig.savefig(buf,format='png',dpi=600,facecolor='white');buf.seek(0)
    high=Image.open(buf).convert('RGB')
    new_h=round(KEEP_AB*high.width/W)
    high.paste(ab.resize((high.width,new_h),Image.Resampling.LANCZOS),(0,0))
    high.save(stem.with_suffix('.png'),dpi=(600,600))
    high.save(stem.with_suffix('.tif'),dpi=(600,600),compression='tiff_lzw')
    bounds={'C':(0,1040,1355,H-1040),'D':(1360,1040,W-1360,H-1040),'CD':(0,1040,W,H-1040)}
    for panel,b in bounds.items():
        keep={'panel_C','panel_D'} if panel=='CD' else {'panel_'+panel}
        svg_crop(stem.with_suffix('.svg'),OUT/f'Figure4{panel}_updated.svg',b,keep)
        x,y,w,h=b
        crop=native.crop((x,y,x+w,y+h))
        crop=crop.resize((round(w*600/NATIVE_DPI),round(h*600/NATIVE_DPI)),Image.Resampling.LANCZOS)
        crop.save(OUT/f'Figure4{panel}_updated.png',dpi=(600,600))
    ax_ab.set_visible(False)
    fig.savefig(QA/'CD_vector_canvas.pdf',facecolor='white')
    plt.close(fig)
    with (OUT/'Figure4CD_source_data.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(c[0]));writer.writeheader();writer.writerows(c+d)
    (OUT/'Figure4_C_D_caption_updated.txt').write_text(
        '(C) On the frozen 34-pair cis/trans cohort selected without target information, '
        'restoring donor relations lowers the CoordRep-hybrid dipole pair-difference MAE '
        'from 5.198 to 1.340 D (74.2%) and raises direction concordance from 73.5% to 94.1%. '
        'GINE, SchNet, and ViSNet give MAEs of 5.869, 0.763, and 0.530 D, respectively. '
        '(D) In 44 graph-identical, different-shape pairs with |delta E_gap| > 0.10 eV, '
        'restoring shape/CShM lowers MAE from 0.358 to 0.255 eV (28.7%) and raises '
        'concordance from 76.1% to 88.6%; GINE, SchNet, and ViSNet give 0.536, 0.351, '
        'and 0.400 eV. Error bars are pointwise 95% attributed-graph-cluster bootstrap '
        'intervals conditional on the frozen three-seed ensembles.\n',encoding='utf-8')
    audit=dict(status='NUMERIC_AND_AB_PRESERVATION_PASS',
        source_image_sha256=digest(args.ab_image),
        source_image_size=[W,H],ab_preserved_native_rectangle=[0,0,W,KEEP_AB],
        ab_native_pixels_identical=True,ab_redrawn=False,all_new_text_font_family=fonts,
        arial_font_files={p.name:digest(p) for p in font_paths},
        high_resolution_size=list(high.size),raster_dpi=600,
        exact_plot_artist_positions_and_asymmetric_CI_checked=True,
        C=c,D=d,relative_reduction_C=1-c[2]['value']/c[1]['value'],
        relative_reduction_D=1-d[2]['value']/d[1]['value'],
        source_hashes={'Figure4CD_source_data.csv':digest(args.source)},
        pdf_crop_bounds_pixels=bounds,
        graphic_scope='C/D editable vector; A/B original submitted raster retained without redraw')
    json_write(OUT/'FIGURE4_UPDATE_AUDIT.json',audit)
    print(json.dumps({'status':audit['status'],'native':native.size,'high':high.size,'fonts':fonts}))

if __name__=='__main__':main()
