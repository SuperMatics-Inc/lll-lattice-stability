"""Render the immutable certified PARTIAL map; no new spectral calculation.
Run: uv run --no-project --with-requirements requirements.txt python figures/plot_partial_map.py
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib
import importlib.util
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import Rectangle, Polygon, Patch, PathPatch
from shapely.geometry import box
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
FIG = Path(__file__).resolve().parent
plt.style.use(FIG / 'pub.mplstyle')
plt.rcParams['hatch.linewidth'] = .35
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path, default=FIG)
OUT = parser.parse_args().output_dir
OUT.mkdir(parents=True, exist_ok=True)
source = ROOT / 'gstmap/data/tree_smallq.jsonl'
assert hashlib.sha256(source.read_bytes()).hexdigest() == '647585a453c630987b8ee98a03f2af03065d6e6598636f11e4dd2f18e49012ee'
spec = importlib.util.spec_from_file_location('gen_data', ROOT / 'gstmap/gen_data.py')
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)
nodes, errors, leaves = g.load(source), [], []
for k in range(4):
    for j in range(19):
        _, part = g.build(nodes, (F(1,16)+F(k,8), F(7+j,8), F(1,16)), errors)
        leaves.extend(part)
assert not errors
counts = Counter(kind for kind, _ in leaves)
assert counts == dict(K=66, S=4178, U=25088, X=11238)
areas = {kind: str(sum((4*w*w for kd, (_,_,w) in leaves if kd==kind), F(0))) for kind in counts}
# Dissolve internal dyadic edges before vector export to avoid PDF hairline seams.
# Every coordinate is binary-exact; the certified leaf data remain unchanged.
paths = {}
for kind in ['U', 'S', 'X']:
    region = unary_union([box(float(a-w), float(b-w), float(a+w), float(b+w))
                          for kd, (a, b, w) in leaves if kd == kind])
    assert region.area == float(F(areas[kind]))
    polygons = [region] if region.geom_type == 'Polygon' else region.geoms
    vertices, codes = [], []
    for poly in polygons:
        poly = orient(poly, sign=1)
        for ring in [poly.exterior, *poly.interiors]:
            xy = list(ring.coords)
            vertices.extend(xy)
            codes.extend([MplPath.MOVETO] + [MplPath.LINETO]*(len(xy)-2) + [MplPath.CLOSEPOLY])
    paths[kind] = MplPath(vertices, codes)
fig, axes = plt.subplots(1, 2, figsize=(6.3, 3.7), gridspec_kw={'width_ratios':[0.85,1.15]}, layout='constrained')
styles = {'U':dict(facecolor='#F6EEE8', edgecolor='#A77755', hatch='/', linewidth=0),
          'S':dict(facecolor='#0072B2', edgecolor='none', linewidth=0),
          'X':dict(facecolor='#000000', edgecolor='none', linewidth=0)}
x = np.linspace(0, .5, 2001)
arc = np.sqrt(1-x*x)
for ax, label in zip(axes, ['(a)', '(b)']):
    # A common clip path restricts every exact rectangle to F+.
    clip = Polygon(np.c_[np.r_[x,.5,0],np.r_[arc,3.2,3.2]], closed=True, transform=ax.transData)
    for kind in ['U','S','X']:
        collection = PathPatch(paths[kind], **styles[kind])
        collection.set_clip_path(clip)
        ax.add_patch(collection)
        collection.set_clip_box(ax.bbox)
    ax.plot(x, arc, color='0.2', linestyle='--', linewidth=.8)
    ax.axhline(3, color='0.2', linestyle=':', linewidth=.9)
    ax.set(xlim=(0,.5), xlabel=r'$x=\operatorname{Re}\tau$', ylabel=r'$y=\operatorname{Im}\tau$')
    ax.text(.035,.965,label,transform=ax.transAxes,ha='left',va='top',fontsize=8,bbox=dict(facecolor='white',edgecolor='none',pad=1.5))
axes[0].set_ylim(.84,3.19)
axes[0].add_patch(Rectangle((.28,.86),.22,.22,fill=False,edgecolor='0.25',linestyle='-.',linewidth=.8))
axes[1].set(xlim=(.28,.5),ylim=(.86,1.08))
axes[1].set_aspect('equal', adjustable='box')
axes[1].plot([.4],[1],marker='o',color='white',markeredgecolor='black',linestyle='none')
axes[1].annotate(r'$\tau_+$',(.4,1),xytext=(5,7),textcoords='offset points',fontsize=8)
axes[1].plot([.5],[np.sqrt(3)/2],marker='^',color='white',markeredgecolor='black',linestyle='none',clip_on=False)
legend=[Patch(facecolor='#0072B2',label='S: stable'),Patch(facecolor='#F6EEE8',edgecolor='#A77755',hatch='/',label='U: unstable'),Patch(facecolor='black',label='X: undecided')]
fig.legend(handles=legend,loc='outside lower center',ncol=3)
fig.savefig(OUT/'partial-map.pdf',metadata={'CreationDate':None,'ModDate':None})
fig.savefig(OUT/'partial-map.png',dpi=200)
receipt={'claim':'CERTIFIED PARTIAL map; rendering only, no new physics', 'tree_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'leaf_counts':dict(counts),'coordinate_areas':areas,'style_sha256':hashlib.sha256((FIG/'pub.mplstyle').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pdf_sha256':hashlib.sha256((OUT/'partial-map.pdf').read_bytes()).hexdigest()}
(OUT/'partial-map-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
