#!/usr/bin/env python3
"""Build Figure 2 (178 mm wide, mSystems full width) from exported TSVs, using only the Python standard library.
Run: python3 -I manuscript_figures/figure2/build_figure2.py
Polar SVG primitives follow .codeome/context/sled_fig2/build_figure2.py.
All coordinates are mm and text sizes are points. No statistical analysis here.
"""
import csv
import re
import colorsys
import math
import os
import subprocess
import textwrap
from collections import defaultdict
from pathlib import Path
from xml.sax.saxutils import escape

OUT = Path(__file__).resolve().parent
INKSCAPE = '/Applications/Inkscape.app/Contents/MacOS/inkscape'
COL = {'Podarcis_liolepis': '#d57d2c', 'Podarcis_muralis': '#008080'}
W=178
o = ['HEADER', 'BG']
def read(name):
    with (OUT / (name + '.tsv')).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))
# Match only the two named lizard species, never taxonomy/function labels.
SPECIES_TEXT=re.compile(r'(?<![A-Za-z])(?:P\. |Podarcis )(?:liolepis|muralis)(?![A-Za-z])')
def styled_text(value):
    value=str(value); parts=[]; end=0
    for match in SPECIES_TEXT.finditer(value):
        parts.append(escape(value[end:match.start()]))
        parts.append('<tspan font-style="italic">'+escape(match.group())+'</tspan>')
        end=match.end()
    parts.append(escape(value[end:]))
    return ''.join(parts)
def label_width(value,mm_per_char):
    # Conservative ink box includes italic overhang in addition to advance width.
    return len(value)*mm_per_char+(.6 if SPECIES_TEXT.search(value) else 0)
def text(x,y,s,size=6,anchor='start',weight='normal',fill='#222'):
    assert size >= 6
    o.append(f'<text x="{x:.3f}" y="{y:.3f}" font-family="Helvetica, Arial, sans-serif" font-size="{size*25.4/72:.4f}" text-anchor="{anchor}" font-weight="{weight}" fill="{fill}">{styled_text(s)}</text>')
def rect(x,y,w,h,c):
    o.append(f'<rect x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" fill="{c}"/>')
def path(d,c='none',stroke='none',sw=.2,opacity=1):
    o.append(f'<path d="{d}" fill="{c}" stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"/>')
def line(x,y,xx,yy,c='#333',sw=.2,opacity=1):
    path(f'M{x:.3f} {y:.3f}L{xx:.3f} {yy:.3f}',stroke=c,sw=sw,opacity=opacity)
def ramp(t):
    t=max(0,min(1,t)); a=(209,244,186); b=(244,186,186)
    return '#%02x%02x%02x'%tuple(round(x+(y-x)*t) for x,y in zip(a,b))
CX,CY=49,56

def pol(r,a):
    a=math.radians(a); return CX+r*math.sin(a),CY-r*math.cos(a)
def arc(r,a,b):
    x,y=pol(r,b); return f'A{r:.4f} {r:.4f} 0 {int(abs(b-a)>180)} {int(b>=a)} {x:.4f} {y:.4f}'
def wedge(r,s,a,b):
    x,y=pol(s,a); xx,yy=pol(r,b)
    return f'M{x:.4f} {y:.4f}{arc(s,a,b)}L{xx:.4f} {yy:.4f}{arc(r,b,a)}Z'
def panel(letter,x,y,title):
    text(x,y,letter,10,weight='bold'); text(x+5,y,title,7,weight='bold')
def gradient(gid,stops):
    o.append(f'<defs><linearGradient id="{gid}">'+''.join(f'<stop offset="{j/(len(stops)-1)}" stop-color="{c}"/>' for j,c in enumerate(stops))+'</linearGradient></defs>')
def red(t):
    t=max(0,min(1,t)); stops=['#f4baba','#e0504d','#7f1010']
    k=min(int(t*2),1); u=t*2-k
    a=[int(stops[k][j:j+2],16) for j in (1,3,5)]
    b=[int(stops[k+1][j:j+2],16) for j in (1,3,5)]
    return '#%02x%02x%02x'%tuple(round(x+(y-x)*u) for x,y in zip(a,b))

panel('a',2,6,'MAG catalogue')
nodes={int(r['id']):r for r in read('tree_nodes')}
meta={r['genome']:r for r in read('mag_metadata')}
pcol={r['phylum']:r['colour'] for r in read('phylum_colors')}
counts={r['phylum']:int(r['n']) for r in read('phylum_counts')}
nspecies={r['species']:int(r['n']) for r in read('species_counts')}
phyla=sorted(counts,key=lambda p:(-counts[p],p))
tips={i:int(r['tip_pos']) for i,r in nodes.items() if r['tip_pos']!='NA'}
kids=defaultdict(list)
for i,r in nodes.items():
    if r['parent']!='NA': kids[int(r['parent'])].append(i)
root=next(i for i,r in nodes.items() if r['parent']=='NA')
post=[]; stack=[root]
while stack:
    i=stack.pop(); post.append(i); stack.extend(kids[i])
# Open gap lies to the upper left, allowing ring names to end at its right edge.
step=310/(len(tips)-1); ang={i:p*step for i,p in tips.items()}
sub={}
for i in reversed(post):
    if i not in tips: ang[i]=(min(ang[c] for c in kids[i])+max(ang[c] for c in kids[i]))/2
    sub[i]={meta[nodes[i]['label']]['phylum']} if i in tips else set().union(*(sub[c] for c in kids[i]))
branch={i:pcol[next(iter(v))] if len(v)==1 else '#999' for i,v in sub.items()}
maxd=max(float(r['depth']) for r in nodes.values())
rad={i:2+math.sqrt(float(r['depth'])/maxd)*24 for i,r in nodes.items()}
for i,cs in kids.items():
    if not cs: continue
    x,y=pol(rad[i],min(ang[c] for c in cs))
    path(f'M{x} {y}'+arc(rad[i],min(ang[c] for c in cs),max(ang[c] for c in cs)),stroke=branch[i],sw=.12)
    for c in cs: line(*pol(rad[i],ang[c]),*pol(rad[c],ang[c]),branch[c],.12)
maxmb=max(float(m['length'])/1e6 for m in meta.values())
gc=[float(m['gc']) for m in meta.values()]; gclo,gchi=min(gc),max(gc)
for mb in (2,4,6):
    r=36+mb*.65; x,y=pol(r,0)
    path(f'M{x} {y}'+arc(r,0,310),stroke='#ccc',sw=.1)
ordered=sorted(tips,key=tips.get)
for i in ordered:
    m=meta[nodes[i]['label']]; a,b=ang[i]-step*.51,ang[i]+step*.51
    path(wedge(27,28.5,a,b),'#1e6e55' if m['novel']=='TRUE' else '#ededed')
    # Bar length retains completeness (50–100%); colour is contamination (0–10%).
    path(wedge(29.2,29.2+(float(m['completeness'])-50)/50*2.8,a,b),red(float(m['contamination'])/10))
    path(wedge(36,36+float(m['length'])/1e6*.65,a,b),'#999')
    pv=[float(m['prev_'+sp]) for sp in COL]; total=sum(pv); r0=41.4
    if total:
        for sp,v in zip(COL,pv):
            r1=r0+v/total*2.5
            if r1>r0: path(wedge(r0,r1,a,b),COL[sp])
            r0=r1
    else: path(wedge(41.4,43.9,a,b),'#eee')
# Unsmoothed GC values: exact values, joined in tree order.
gpts=[pol(32.8+(float(meta[nodes[i]['label']]['gc'])-gclo)/(gchi-gclo)*2.5,ang[i]) for i in ordered]
linepts='L'.join(f'{x:.4f} {y:.4f}' for x,y in gpts)
x0,y0=pol(32.8,0); x1,y1=pol(32.8,310)
path(f'M{x0} {y0}L{linepts}L{x1} {y1}'+arc(32.8,310,0)+'Z','#cfe8cf')
path('M'+linepts,stroke='#1f7a3a',sw=.15)
for lab,r in [('Novel species',27.7),('Genome quality',30.5),('GC content',34),('Genome size',38.5),('Species composition',42.6)]:
    text(CX-1,CY-r+.7,lab,6,'end')
# Largest contiguous run per phylum. Long runs receive rotated labels on the
# inside of the clade; small-clade names have leaders to a collision-spaced list.
runs=[]; run=[]
for i in ordered:
    p=meta[nodes[i]['label']]['phylum']
    if run and p!=run[0]: runs.append(run); run=[]
    if not run: run=[p,[]]
    run[1].append(i)
runs.append(run)
best={}
for p,ids in runs:
    if p not in best or len(ids)>len(best[p]): best[p]=ids
small=[]
for p,ids in best.items():
    a=(ang[ids[0]]+ang[ids[-1]])/2; span=len(ids)*step
    if span>32 and p!='p__Pseudomonadota':
        # Curved phylum labels sit just inside the tips, upright on both halves.
        delta=max(span*.43,len(p.removeprefix('p__'))*1.1/23*180/math.pi/2+2)
        a0,a1=a-delta,a+delta
        if 90<a<270: a0,a1=a1,a0
        x,y=pol(23,a0); gid='clade-'+str(ids[0])
        o.append(f'<defs><path id="{gid}" d="M{x} {y}'+arc(23,a0,a1)+'"/></defs>')
        o.append(f'<text font-family="Helvetica, Arial, sans-serif" font-size="{6*25.4/72}" font-weight="bold" fill="{pcol[p]}" stroke="white" stroke-width=".5" paint-order="stroke"><textPath href="#{gid}" startOffset="50%" text-anchor="middle">'+escape(p.removeprefix('p__'))+'</textPath></text>')
    else: small.append((p,a))
# Place small-clade labels in the inner fan with short leaders. Greedy placement
# favours the nearest free box; all boxes stay inside the tree and avoid each other.
occupied=[]
leaders=[]; labels=[]
for p,a in sorted(small,key=lambda v:-len(v[0])):
    name=p.removeprefix('p__'); width=len(name)*1.04
    tx,ty=pol(18,a); options=[]
    for yy in [CY-16+j*2.8 for j in range(13)]:
        for xx in [CX-17+j for j in range(35)]:
            box=(xx-width/2-.4,yy-1.9,xx+width/2+.4,yy+.5)
            if any(math.hypot(bx-CX,by-CY)>18.5 for bx in (box[0],box[2]) for by in (box[1],box[3])): continue
            if any(box[0]<b[2] and box[2]>b[0] and box[1]<b[3] and box[3]>b[1] for b in occupied): continue
            options.append(((xx-tx)**2+(yy-ty)**2,xx,yy,box))
    if not options: raise ValueError('Cannot place clade label '+name)
    _,xx,yy,box=min(options); occupied.append(box)
    px,py=pol(25.6,a); ex=max(box[0],min(px,box[2])); ey=max(box[1],min(py,box[3]))
    leaders.append((px,py,ex,ey,pcol[p])); labels.append((xx,yy,name,pcol[p]))
for px,py,ex,ey,c in leaders: line(px,py,ex,ey,c,.16)
for xx,yy,name,c in labels:
    o.append(f'<text x="{xx}" y="{yy}" font-family="Helvetica, Arial, sans-serif" font-size="{6*25.4/72}" text-anchor="middle" fill="{c}" stroke="white" stroke-width=".65" paint-order="stroke">'+escape(name)+'</text>')
text(50,105,f'{len(tips)} MAGs; {sum(m["novel"]=="TRUE" for m in meta.values())} novel species',6,'middle')
print(f'PANEL a: {len(tips)} tips; five rings; {len(best)} clade labels; square-root branches')

panel('b',101,6,'Wild composition')
comp=read('sample_composition'); bysample=defaultdict(dict); species={}
for r in comp:
    bysample[r['sample']][r['phylum']]=float(r['rel']); species[r['sample']]=r['species']
tiles=defaultdict(list)
for r in read('sample_mag_abundance'): tiles[r['sample']].append(r)
rank={nodes[i]['label']:tips[i] for i in tips}; prank={p:j for j,p in enumerate(phyla)}
dominant=max(phyla,key=lambda p:sum(v.get(p,0) for v in bysample.values()))
# Place the ordering phylum at the bottom, like the reference stacked bars.
prank={p:j for j,p in enumerate([dominant]+[p for p in phyla if p!=dominant])}
B_X0,B_X1=109,175
HEADER_TOP=CY-43.9  # outer species-composition edge at theta=0
# Helvetica/Arial capital height at 6 pt; checked against Inkscape's ink bounds.
HEADER_BASELINE=HEADER_TOP+1.543
x=B_X0; y=HEADER_TOP+5; h=72-y; bw=(B_X1-B_X0-3)/len(bysample)
for sp in COL:
    ids=sorted((s for s in bysample if species[s]==sp),key=lambda s:(-bysample[s].get(dominant,0),s)); w=len(ids)*bw
    o.append(f'<g id="header-{sp}">')
    text(x+w/2,HEADER_BASELINE,sp.replace('Podarcis_','P. ')+f' (n={len(ids)})',6,'middle',fill=COL[sp])
    o.append('</g>'); line(x,y-1.4,x+w,y-1.4,COL[sp],.5)
    for j,s in enumerate(ids):
        yy=y+h
        assert abs(sum(float(r['rel']) for r in tiles[s])-1)<1e-10
        for r in sorted(tiles[s],key=lambda r:(prank[r['phylum']],rank.get(r['genome'],10**6))):
            hh=float(r['rel'])*h; yy-=hh
            path(f'M{x+j*bw:.4f} {yy:.4f}h{bw:.4f}v{hh:.4f}h{-bw:.4f}Z',pcol[r['phylum']],'#fff',.035)
    line(x,y+h,x+w,y+h); x+=w+3
for v in (0,50,100): text(B_X0-1.4,y+h*(1-v/100)+.7,v,6,'end')
line(B_X0,y,B_X0,y+h)
o.append('<g transform="translate(101 44.5) rotate(-90)">'); text(0,0,'Relative abundance (%)',6,'middle'); o.append('</g>')
text((B_X0+B_X1)/2,77,'Wild samples',6,'middle')
print(f'PANEL b: {sum(map(len,tiles.values()))} nonzero MAG tiles; ordering by {dominant}')

# Panel b ends at the same baseline as panel a's catalogue caption.
TOP_ROW_BOTTOM=105
text(101,82,'Phylum (number of MAGs)',6,weight='bold')
for j,p in enumerate(phyla):
    xx=101+(j//5)*25; yy=TOP_ROW_BOTTOM-4*4.5+(j%5)*4.5
    o.append(f'<g id="phylum-key-{j}">')
    rect(xx,yy-2,2,2,pcol[p]); text(xx+3,yy,p.removeprefix('p__')+f' ({counts[p]})'); o.append('</g>')
print(f'ALIGNMENT: header top = outer ring top {HEADER_TOP:.2f} mm; bar height {h:.2f} mm')
print(f'ALIGNMENT: a caption and b legend bottom = {TOP_ROW_BOTTOM} mm')

# Ring descriptions belong in the manuscript caption; names remain in the fan gap.
print('LEGEND: phylum counts below b; ring descriptions in caption')
o.append('<g>')

# Function palette is matched by the exact prefix before the first underscore.
# Display-only exclusions: retain all significant tests/means in source TSVs.
EXCLUDED_FUNCTIONS=frozenset({'Cellular structure','Appendages'})
all_significant=read('functional_capacity')
excluded=[r for r in all_significant if r['Function'].split('_',1)[0] in EXCLUDED_FUNCTIONS]
fun=[r for r in all_significant if r['Function'].split('_',1)[0] not in EXCLUDED_FUNCTIONS]
assert len(excluded)==2 and {r['trait'] for r in excluded}=={'S0103','S0202'}
print(f'DISPLAY FILTER: {len(all_significant)} significant; excluded {len(excluded)} (Teichoic acid, Pilus); showing {len(fun)} bars')
with (OUT.parents[1]/'data/gift_colors.tsv').open() as f:
    palette_rows=list(csv.DictReader(f,delimiter='\t'))
assert len({r['Function'] for r in palette_rows})==len(palette_rows)
function_col={r['Function']:r['Color'] for r in palette_rows}
# Display-only override (data/gift_colors.tsv is not modified): keeps the warm
# (biosynthesis) / cool (degradation) scheme but separates similar hues by hue and
# lightness. CIELAB dE between the closest displayed pair rises from 10.9 to 27.8
# (normal vision) and from 6.5 to 10.3 (simulated protan/deutan).
DISPLAY_COLOURS={
    'Nucleic acid biosynthesis':'#f9c145','Amino acid biosynthesis':'#f79039',
    'Amino acid derivative biosynthesis':'#b5421a','Organic anion biosynthesis':'#d01040',
    'Vitamin biosynthesis':'#7a1650','Aromatic compound biosynthesis':'#d8359a',
    'Metallophore biosynthesis':'#ff7f9e','Antibiotic biosynthesis':'#f5b3a0',
    'Polysaccharide degradation':'#2a64ad','Sugar degradation':'#62d4ea',
    'Amino acid degradation':'#c7e9a4','Nitrogen compound degradation':'#2fa56b'}
assert set(DISPLAY_COLOURS)<=set(function_col)
function_col.update(DISPLAY_COLOURS)
present_functions={r['Function'].split('_',1)[0] for r in fun}
functions=[r['Function'] for r in palette_rows if r['Function'] in present_functions]
functions+=sorted(present_functions-set(functions))
missing=[f for f in functions if f not in function_col]
# Assign an explicitly reported, unused fallback to each unexpected annotation.
for j,f in enumerate(missing):
    fallback=f'#{(0x5b21b6+j*0x17912b)%0xffffff:06x}'
    while fallback.lower() in {c.lower() for c in function_col.values()}:
        fallback=f'#{(int(fallback[1:],16)+0x10101)%0xffffff:06x}'
    function_col[f]=fallback
    print(f'WARNING: function colour missing for {f!r}; fallback {fallback}')
def plot_label(name):
    function,element=name.split('_',1)
    verb=next((v for v in ('biosynthesis','degradation') if v in function.lower()),'')
    return element+(' '+verb if verb else '')
assert plot_label('Amino acid biosynthesis_Tryptophan')=='Tryptophan biosynthesis'
assert plot_label('Polysaccharide degradation_Arabinan')=='Arabinan degradation'
assert plot_label('Appendages_Pilus')=='Pilus'
assert plot_label('Antibiotic biosynthesis_Fosfomycin')=='Fosfomycin biosynthesis'
label_rows=[dict(trait=r['trait'],Function=r['Function'],function=r['Function'].split('_',1)[0],label=plot_label(r['Function']),colour=function_col[r['Function'].split('_',1)[0]]) for r in fun]
with (OUT/'functional_plot_labels.tsv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['trait','Function','function','label','colour'],delimiter='\t',lineterminator='\n')
    writer.writeheader(); writer.writerows(label_rows)
function_groups={
    'Biosynthesis':[f for f in functions if 'biosynthesis' in f.lower()],
    'Degradation':[f for f in functions if 'degradation' in f.lower()]}
# Warm reds/oranges/pinks versus cool greens/blues, using conventional HSV hue
# sectors (warm <=70 or >=300 degrees; cool 70..260). Report, never recolour.
palette_exceptions=[]
for group in ('Biosynthesis','Degradation'):
    for f in function_groups[group]:
        colour=function_col[f]; rgb=[int(colour[k:k+2],16)/255 for k in (1,3,5)]
        hue,saturation,_=colorsys.rgb_to_hsv(*rgb); hue*=360
        warm=hue<=70 or hue>=300
        cool=70<hue<260
        if saturation<.1 or not (warm if group=='Biosynthesis' else cool):
            palette_exceptions.append((f,colour,round(hue,1)))
print('PALETTE CHECK:',palette_exceptions if palette_exceptions else 'all 8 biosynthesis warm; all 4 degradation cool; no Other functions displayed')
print(f'FUNCTION COLOURS: {len(functions)} exact-name matches; {len(missing)} fallbacks')

# c: source points and species-centroid spokes; no ellipses in the source plot.
panel('c',2,117,'Neutral beta diversity')
pts=read('nmds_points'); x0,y0,w,h=13,127,64,40
xs=[float(r['NMDS1']) for r in pts]; ys=[float(r['NMDS2']) for r in pts]
xmin,xmax=min(xs),max(xs); ymin,ymax=min(ys),max(ys)
xpad=(xmax-xmin)*.08; ypad=(ymax-ymin)*.08
xmin-=xpad; xmax+=xpad; ymin-=ypad; ymax+=ypad
X=lambda v:x0+(v-xmin)/(xmax-xmin)*w
Y=lambda v:y0+h-(v-ymin)/(ymax-ymin)*h
line(x0,y0,x0,y0+h); line(x0,y0+h,x0+w,y0+h)
for k in range(5):
    v=xmin+(xmax-xmin)*k/4; text(X(v),y0+h+4,f'{v:.2f}',6,'middle')
    v=ymin+(ymax-ymin)*k/4; text(x0-2,Y(v)+.7,f'{v:.2f}',6,'end')
def symbol(x,y,kind,c):
    if kind=='Cold_control':
        o.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r=".9" fill="{c}"/>')
    elif kind=='Warm_control': rect(x-.85,y-.85,1.7,1.7,c)
    elif kind=='Cold_intervention': path(f'M{x} {y-1.1}l1 1.9h-2Z',c)
    else: raise ValueError('Unexpected experimental group: '+kind)
for r in pts:
    line(X(float(r['x_cen'])),Y(float(r['y_cen'])),X(float(r['NMDS1'])),Y(float(r['NMDS2'])),COL[r['species']],.2,.2)
for r in pts: symbol(X(float(r['NMDS1'])),Y(float(r['NMDS2'])),r['type'],COL[r['species']])
text(x0+w/2,y0+h+9,'NMDS1',6,'middle')
o.append(f'<g id="nmds2-axis-label" transform="translate(2.5 {y0+h/2}) rotate(-90)">')
text(0,0,'NMDS2',6,'middle'); o.append('</g>')
text(x0+w/2,y0-3,'Stress = '+format(float(read('nmds_fit')[0]['stress']),'.3f'),6,'middle')
# Place names beside their data-derived centroids in free plot space. Exclude
# point symbols and all spoke segments (with margins), rather than hide them.
def segment_hits_box(a,b,box):
    t0,t1=0.,1.
    for av,bv,lo,hi in ((a[0],b[0],box[0],box[2]),(a[1],b[1],box[1],box[3])):
        delta=bv-av
        if abs(delta)<1e-12:
            if av<lo or av>hi: return False
        else:
            u,v=(lo-av)/delta,(hi-av)/delta
            t0=max(t0,min(u,v)); t1=min(t1,max(u,v))
            if t0>t1: return False
    return True
point_xy=[(X(float(r['NMDS1'])),Y(float(r['NMDS2']))) for r in pts]
spokes=[((X(float(r['x_cen'])),Y(float(r['y_cen']))),point) for r,point in zip(pts,point_xy)]
name_boxes=[]
for sp in COL:
    label=sp.replace('Podarcis_','P. '); width=label_width(label,1.08)
    member=next(r for r in pts if r['species']==sp)
    cx,cy=X(float(member['x_cen'])),Y(float(member['y_cen']))
    options=[]
    for j in range(45):
        yy=y0+3+j*.8
        for k in range(60):
            xx=x0+1+width/2+k
            box=(xx-width/2-.5,yy-1.9,xx+width/2+.5,yy+.6)
            if box[2]>x0+w-1 or box[3]>y0+h-1: continue
            if any(box[0]-1.6<px<box[2]+1.6 and box[1]-1.6<py<box[3]+1.6 for px,py in point_xy): continue
            expanded=(box[0]-.5,box[1]-.5,box[2]+.5,box[3]+.5)
            if any(segment_hits_box(a,b,expanded) for a,b in spokes): continue
            if any(box[0]<b[2] and box[2]>b[0] and box[1]<b[3] and box[3]>b[1] for b in name_boxes): continue
            options.append(((xx-cx)**2+(yy-cy)**2,xx,yy,box))
    if not options: raise ValueError('No free cluster-label position: '+sp)
    _,xx,yy,box=min(options); name_boxes.append(box)
    o.append(f'<g id="cluster-name-{sp}">')
    text(xx,yy,label,6,'middle',fill=COL[sp]); o.append('</g>')
    print(f'CLUSTER LABEL: {label} at ({xx:.2f}, {yy:.2f}) mm; no point/spoke overlaps')
for j,t in enumerate(['Cold_control','Warm_control','Cold_intervention']):
    symbol(x0+1+j*23,y0+h+14.5,t,'#555'); text(x0+3.5+j*23,y0+h+15.2,t.replace('_',' '))
# Heading-free function key spans the full c column, before the gutter at x=88.
legend_functions=function_groups['Biosynthesis']+function_groups['Degradation']
assert len(legend_functions)==12 and len(function_groups['Biosynthesis'])==8
legend_columns=[legend_functions[:6],legend_functions[6:]]
LEGEND_X,LEGEND_W=2,84
LEGEND_TOP,LEGEND_PAD=186,2
LEGEND_ROW=4
legend_height=6*LEGEND_ROW+2*LEGEND_PAD
o.append(f'<rect id="function-legend-box" x="{LEGEND_X}" y="{LEGEND_TOP}" width="{LEGEND_W}" height="{legend_height}" rx="2" fill="#f2f2f2"/>')
for col,entries in enumerate(legend_columns):
    xx=LEGEND_X+LEGEND_PAD+col*(LEGEND_W/2)
    for j,f in enumerate(entries):
        centre=LEGEND_TOP+LEGEND_PAD+(j+.5)*LEGEND_ROW
        o.append(f'<g id="function-key-{functions.index(f)}">')
        rect(xx,centre-1,2,2,function_col[f]); text(xx+2.6,centre+.7,f,6)
        o.append('</g>')
C_BOTTOM=LEGEND_TOP+legend_height+2
print(f'FUNCTION BOX: x={LEGEND_X}..{LEGEND_X+LEGEND_W} mm; 12 single-line entries, 6+6; zero wraps')
print(f'PANEL c: {len(pts)} points; {len(functions)} functions below a 40 mm NMDS plot')

# d: exact plotted differences, sorted using fct_reorder(Function, Difference).
panel('d',88,117,'Functional capacity')
# Direction communicates the higher species; fill communicates GIFT function.
groups=defaultdict(list)
for r in fun: groups[r['Function']].append(float(r['Difference']))
# R fct_reorder uses median by default; equal labels are stacked by geom_col.
def median(v):
    v=sorted(v); n=len(v); return (v[(n-1)//2]+v[n//2])/2
order=sorted(groups,key=lambda s:(median(groups[s]),s))
rows=[]
for name in reversed(order):
    rows.append((name,[plot_label(name)],2.4))
total=sum(r[2] for r in rows)
assert total<122, f'Functional labels need {total} mm; revise layout'
base_left,right=139,175; top=122
original_lo=min(0,min(sum(v for v in groups[s] if v<0) for s in groups)); original_hi=max(0,max(sum(v for v in groups[s] if v>0) for s in groups))
pad=(original_hi-original_lo)*.06
original_lo-=pad; hi=original_hi+pad
# Preserve the previous zero position and positive units/mm exactly. Expand the
# chart leftwards (and move y labels), instead of squeezing the positive bars.
MM_PER_UNIT=(right-base_left)/(hi-original_lo)
ZERO_X=base_left-original_lo*MM_PER_UNIT
F=lambda v:ZERO_X+v*MM_PER_UNIT
lo=-.12; left=F(lo)
assert left>120 and abs(F(hi)-right)<1e-10
for v in (-.12,-.08,-.04,0,.04,.08,.12,.16,.20):
    if v>hi: continue
    line(F(v),top,F(v),top+total,'#ddd',.15)
    text(F(v),top+total+4,'0' if v==0 else f'{v:.2f}',6,'middle')
bar_boxes=[]
line(F(0),top,F(0),top+total,'#555',.2)
y=top
for name,labels,rh in rows:
    for j,label in enumerate(labels): text(left-2,y+rh/2+(j-(len(labels)-1)/2)*2.3+.7,label,6,'end')
    neg=pos=0
    for r in (r for r in fun if r['Function']==name):
        v=float(r['Difference']); start=neg if v<0 else pos
        if v<0: neg+=v
        else: pos+=v
        bx,by,bw,bh=min(F(start),F(start+v)),y+rh*.2,abs(F(start+v)-F(start)),rh*.6
        bar_boxes.append((bx,by,bx+bw,by+bh))
        rect(bx,by,bw,bh,function_col[name.split('_',1)[0]])
    y+=rh
# Direction notes are placed only in the empty opposite-sign corners.
direction_boxes=[]
for gid,xx,yy,anchor,sp,lines in [
    ('liolepis',left+1,top+4.5,'start','Podarcis_liolepis',['← higher in','P. liolepis']),
    ('muralis',right-1,top+total-7,'end','Podarcis_muralis',['higher in','P. muralis →'])]:
    width=max(label_width(t,1.12) for t in lines)
    box=(xx-.4 if anchor=='start' else xx-width-.4,yy-2.1,
         xx+width+.4 if anchor=='start' else xx+.4,yy+3.2)
    assert box[0]>=left and box[2]<=right
    assert all(not(box[0]<b[2] and box[2]>b[0] and box[1]<b[3] and box[3]>b[1]) for b in bar_boxes),gid+' direction note intersects bar'
    direction_boxes.append(box)
    o.append(f'<g id="direction-{gid}" stroke="white" stroke-width=".55" paint-order="stroke" stroke-linejoin="round">')
    for j,label in enumerate(lines): text(xx,yy+j*2.5,label,6,anchor,fill=COL[sp])
    o.append('</g>')
print(f'AXIS: {lo:.2f}..{hi:.5f}; zero x={ZERO_X:.5f} mm; {MM_PER_UNIT:.5f} mm/unit unchanged; both direction boxes clear of all {len(bar_boxes)} bars')
text((left+right)/2,top+total+9,'Mean difference (P. muralis − P. liolepis)',6,'middle')
o.append('</g>')
H=max(top+total+12,C_BOTTOM)
print(f'PANEL d: {len(fun)} significant pathways, {len(rows)} labels, height {total:.1f} mm')
o[0]=f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}mm" height="{H:.1f}mm" viewBox="0 0 {W} {H:.1f}">'
o[1]=f'<rect width="{W}" height="{H:.1f}" fill="white"/>'
o.append('</svg>')
svg=OUT/'figure2.svg'; svg.write_text('\n'.join(o))
# Avoid Quartz GUI registration in headless macOS execution.
env = os.environ.copy(); env['GDK_BACKEND'] = 'x11'
subprocess.run([INKSCAPE,str(svg),'--export-type=png','--export-dpi=160','--export-filename='+str(OUT/'figure2_preview.png')],check=True,env=env)
FINAL=OUT.parent/'final'; FINAL.mkdir(exist_ok=True)
subprocess.run([INKSCAPE,str(svg),'--export-type=pdf','--export-text-to-path','--export-filename='+str(FINAL/'figure2.pdf')],check=True,env=env)
subprocess.run([INKSCAPE,str(svg),'--export-type=png','--export-dpi=300','--export-filename='+str(FINAL/'figure2.png')],check=True,env=env)
subprocess.run(['magick',str(FINAL/'figure2.png'),'-background','white','-alpha','remove','-alpha','off','-compress','LZW',str(FINAL/'figure2.tiff')],check=True)
# Render the exported PDF independently through Ghostscript (bundled with Inkscape).
gs='/Applications/Inkscape.app/Contents/Resources/bin/gs'
subprocess.run([gs,'-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=png16m','-dTextAlphaBits=4','-dGraphicsAlphaBits=4','-r160','-sOutputFile='+str(OUT/'figure2_pdf_preview.png'),str(FINAL/'figure2.pdf')],check=True,env=env)
print(f'WROTE SVG, previews, outlined PDF, 300 dpi PNG and LZW TIFF; {W} × {H:.1f} mm')
