# One-time source rewrite helper; retained as a job intermediate.
from pathlib import Path
p=Path(__file__).with_name('build_figure2.py')
s=p.read_text(); start=s.index('TCX,TCY,TS='); end=s.index('# c: source points')
top='''CX,CY=49,56

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
rad={i:2+math.sqrt(float(r['depth'])/maxd)*21 for i,r in nodes.items()}
for i,cs in kids.items():
    if not cs: continue
    x,y=pol(rad[i],min(ang[c] for c in cs))
    path(f'M{x} {y}'+arc(rad[i],min(ang[c] for c in cs),max(ang[c] for c in cs)),stroke=branch[i],sw=.12)
    for c in cs: line(*pol(rad[i],ang[c]),*pol(rad[c],ang[c]),branch[c],.12)
maxmb=max(float(m['length'])/1e6 for m in meta.values())
gc=[float(m['gc']) for m in meta.values()]; gclo,gchi=min(gc),max(gc)
for mb in (2,4,6):
    r=33+mb*.65; x,y=pol(r,0)
    path(f'M{x} {y}'+arc(r,0,310),stroke='#ccc',sw=.1)
    xx,yy=pol(r,313); text(xx-1,yy+.6,f'{mb}',6,'end',fill='#777')
ordered=sorted(tips,key=tips.get)
for i in ordered:
    m=meta[nodes[i]['label']]; a,b=ang[i]-step*.51,ang[i]+step*.51
    path(wedge(24,25.5,a,b),'#1e6e55' if m['novel']=='TRUE' else '#ededed')
    # Bar length retains completeness (50–100%); colour is contamination (0–10%).
    path(wedge(26.2,26.2+(float(m['completeness'])-50)/50*2.8,a,b),red(float(m['contamination'])/10))
    path(wedge(33,33+float(m['length'])/1e6*.65,a,b),'#999')
    pv=[float(m['prev_'+sp]) for sp in COL]; total=sum(pv); r0=38.4
    if total:
        for sp,v in zip(COL,pv):
            r1=r0+v/total*2.5
            if r1>r0: path(wedge(r0,r1,a,b),COL[sp])
            r0=r1
    else: path(wedge(38.4,40.9,a,b),'#eee')
    if float(m['prev_all'])>0: path(wedge(41.7,41.7+float(m['prev_all'])/100*2.8,a,b),'#777')
# Unsmoothed GC values: exact values, joined in tree order.
gpts=[pol(29.8+(float(meta[nodes[i]['label']]['gc'])-gclo)/(gchi-gclo)*2.5,ang[i]) for i in ordered]
linepts='L'.join(f'{x:.4f} {y:.4f}' for x,y in gpts)
x0,y0=pol(29.8,0); x1,y1=pol(29.8,310)
path(f'M{x0} {y0}L{linepts}L{x1} {y1}'+arc(29.8,310,0)+'Z','#cfe8cf')
path('M'+linepts,stroke='#1f7a3a',sw=.15)
for lab,r in [('Novel species',24.7),('Genome quality',27.5),('GC content',31),('Genome size',35.5),('Species composition',39.6),('Overall prevalence',43.5)]:
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
    if span>32:
        x,y=pol(17.5,a); rot=a-90
        if 90<a<270: rot+=180
        o.append(f'<g transform="rotate({rot:.2f} {x:.3f} {y:.3f})">')
        text(x,y,p.removeprefix('p__'),6,'middle',weight='bold',fill=pcol[p]); o.append('</g>')
    else: small.append((p,a))
# Leader labels sit below the fan: no labels over rings or branches.
for j,(p,a) in enumerate(sorted(small,key=lambda v:pol(23,v[1])[0])):
    col=j%3; row=j//3; tx=3+col*32; ty=103+row*3.5
    xx,yy=pol(23,a)
    line(xx,yy,tx+13,ty-2,pcol[p],.12,.6)
    text(tx,ty,p.removeprefix('p__'),6,fill=pcol[p])
text(50,98,f'{len(tips)} MAGs',6,'middle')
print(f'PANEL a: {len(tips)} tips; six rings; {len(best)} clade labels; square-root branches')

panel('b',101,6,'Wild composition')
comp=read('sample_composition'); bysample=defaultdict(dict); species={}
for r in comp:
    bysample[r['sample']][r['phylum']]=float(r['rel']); species[r['sample']]=r['species']
tiles=defaultdict(list)
for r in read('sample_mag_abundance'): tiles[r['sample']].append(r)
rank={nodes[i]['label']:tips[i] for i in tips}; prank={p:j for j,p in enumerate(phyla)}
dominant=max(phyla,key=lambda p:sum(v.get(p,0) for v in bysample.values()))
x=115; y=27; h=63; bw=(60-3)/len(bysample)
for sp in COL:
    ids=sorted((s for s in bysample if species[s]==sp),key=lambda s:(-bysample[s].get(dominant,0),s)); w=len(ids)*bw
    text(x+w/2,18,sp.replace('Podarcis_','P. '),6,'middle',fill=COL[sp])
    text(x+w/2,22,f'n = {len(ids)}',6,'middle',fill=COL[sp]); line(x,24,x+w,24,COL[sp],.5)
    for j,s in enumerate(ids):
        yy=y+h
        assert abs(sum(float(r['rel']) for r in tiles[s])-1)<1e-10
        for r in sorted(tiles[s],key=lambda r:(prank[r['phylum']],rank.get(r['genome'],10**6))):
            hh=float(r['rel'])*h; yy-=hh
            path(f'M{x+j*bw:.4f} {yy:.4f}h{bw:.4f}v{hh:.4f}h{-bw:.4f}Z',pcol[r['phylum']],'#fff',.035)
    line(x,y+h,x+w,y+h); x+=w+3
for v in (0,50,100): text(113,y+h*(1-v/100)+.7,v,6,'end')
line(115,y,115,y+h)
o.append('<g transform="translate(106 59) rotate(-90)">'); text(0,0,'Relative abundance (%)',6,'middle'); o.append('</g>')
text(145,95,'Wild samples',6,'middle')
print(f'PANEL b: {sum(map(len,tiles.values()))} nonzero MAG tiles; ordering by {dominant}')

# One unified legend block, spanning the page below the top row.
LEGEND_Y=121
text(3,LEGEND_Y,'Phylum (number of MAGs)',7,weight='bold')
for j,p in enumerate(phyla):
    xx=3+(j%4)*44; yy=LEGEND_Y+5+(j//4)*4
    rect(xx,yy-2,2.5,2.5,pcol[p]); text(xx+3.5,yy,p.removeprefix('p__')+f' ({counts[p]})')
text(3,145,'Rings (inside to outside)',7,weight='bold')
novel=sum(m['novel']=='TRUE' for m in meta.values())
text(3,151,'Species novelty',6,weight='bold')
rect(3,154,2.5,2.5,'#1e6e55'); text(7,156.2,f'Novel ({novel})')
rect(25,154,2.5,2.5,'#ededed'); text(29,156.2,f'Known ({len(meta)-novel})')
text(3,163,'Genome quality',6,weight='bold')
text(3,167,'Completeness: bar height, 50–100%')
gradient('quality',['#f4baba','#e0504d','#7f1010']); rect(3,170,30,2.5,'url(#quality)')
text(3,176,'Contamination: 0–10%')
text(62,151,'GC content',6,weight='bold')
path('M62 157l4 -2 4 1 4 -2 4 2',stroke='#1f7a3a',sw=.25)
text(82,157,f'{gclo:.2f}–{gchi:.2f}')
text(62,163,'Genome size',6,weight='bold'); rect(62,166,20,2.5,'#999')
text(62,173,'Gridlines every 2 Mb')
text(120,151,'Species composition',6,weight='bold')
for j,sp in enumerate(COL):
    rect(120,154+j*4,2.5,2.5,COL[sp]); text(124,156.2+j*4,sp.replace('Podarcis_','P. ')+f' (n = {nspecies[sp]})')
text(120,167,'Overall prevalence',6,weight='bold'); rect(120,170,20,2.5,'#777')
text(120,176,'0–100% of wild samples')
print('LEGEND: unified four-column phylum key and three-column rings key')
o.append('<g transform="translate(0 68)">')

'''
s=s[:start]+top+s[end:]
s=s.replace('H=top+total+12',"o.append('</g>')\nH=top+total+12+68")
s=s.replace("print('WROTE figure2.svg and figure2_preview.png (178 mm wide; minimum text 6 pt)')",'''FINAL=OUT.parent/'final'; FINAL.mkdir(exist_ok=True)
subprocess.run([INKSCAPE,str(svg),'--export-type=pdf','--export-text-to-path','--export-filename='+str(FINAL/'figure2.pdf')],check=True,env=env)
subprocess.run([INKSCAPE,str(svg),'--export-type=png','--export-dpi=300','--export-filename='+str(FINAL/'figure2.png')],check=True,env=env)
subprocess.run(['magick',str(FINAL/'figure2.png'),'-background','white','-alpha','remove','-alpha','off','-compress','LZW',str(FINAL/'figure2.tiff')],check=True)
# Render the exported PDF independently through Ghostscript (bundled with Inkscape).
gs='/Applications/Inkscape.app/Contents/Resources/bin/gs'
subprocess.run([gs,'-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=png16m','-r160','-sOutputFile='+str(OUT/'figure2_pdf_preview.png'),str(FINAL/'figure2.pdf')],check=True,env=env)
print(f'WROTE SVG, previews, outlined PDF, 300 dpi PNG and LZW TIFF; {W} × {H:.1f} mm')''')
p.write_text(s)
