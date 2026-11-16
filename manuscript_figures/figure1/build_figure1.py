#!/usr/bin/env python3
"""Manuscript Figure 1: experimental design (a) and sampling timeline (b).

Run from the repository root:
    python3 -I manuscript_figures/figure1/build_figure1.py

Drawn at final size (180 mm wide, units in mm). Icons are the supplied PNG assets
(assets/), everything else is vector. Exports SVG, outlined PDF, 300 dpi PNG and TIFF
to manuscript_figures/final/ via Inkscape and ImageMagick.
"""
import base64, os, struct, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ASSETS = HERE / 'assets'
FINAL = HERE.parent / 'final'
INK = '/Applications/Inkscape.app/Contents/MacOS/inkscape'
FONT = "Helvetica, Arial, sans-serif"

W = 180.0
HL = 1.64   # arrowhead length (mm): lines stop here so the head's flat base sits on the line end
# Species colours (as in Fig. 2) and group colours (new: three distinct hues)
ORANGE, TEAL = '#d57d2c', '#008080'
WC, CI, CC = '#d57d2c', '#7b3294', '#008080'          # warm-control, cold-intervention, cold-control
WC_T, CI_T, CC_T = '#fbeedd', '#eadff1', '#d9eeee'    # light tints
GREY, DARK, MID = '#f1f1f1', '#222222', '#666666'

o = []


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def rich(s):
    """Text with *italic* segments (species names)."""
    parts = s.split('*')
    out = ''
    for i, p in enumerate(parts):
        if not p:
            continue
        out += f'<tspan font-style="italic">{esc(p)}</tspan>' if i % 2 else esc(p)
    return out


def text(x, y, s, size=2.2, anchor='start', weight='normal', fill=DARK, rot=None):
    tr = f' transform="rotate({rot} {x} {y})"' if rot is not None else ''
    o.append(f'<text x="{x:.2f}" y="{y:.2f}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
             f'fill="{fill}" text-anchor="{anchor}"{tr}>{rich(s)}</text>')


def rect(x, y, w, h, fill='none', stroke='none', sw=.3, rx=0, opacity=1):
    o.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" fill="{fill}" '
             f'stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"/>')


def line(x1, y1, x2, y2, c='#999', sw=.25, dash=None, marker=None):
    d = f' stroke-dasharray="{dash}"' if dash else ''
    m = f' marker-end="url(#{marker})"' if marker else ''
    o.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{c}" stroke-width="{sw}"{d}{m}/>')


def path(d, stroke, sw=.3, marker=None, opacity=1, fill='none'):
    m = f' marker-end="url(#{marker})"' if marker else ''
    o.append(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"{m}/>')


_cache = {}
used = set()


def asset(name):
    if name not in _cache:
        data = (ASSETS / f'{name}.png').read_bytes()
        w, h = struct.unpack('>II', data[16:24])
        _cache[name] = (base64.b64encode(data).decode(), h / w)
    return _cache[name]


def image(name, x, y, w, rot=None, around=None):
    b64, ar = asset(name)
    used.add(name)
    tr = f'rotate({rot} {around[0]:.2f} {around[1]:.2f}) ' if rot is not None else ''
    o.append(f'<use href="#img_{name}" transform="{tr}translate({x:.2f} {y:.2f}) scale({w:.3f})"/>')
    return w * ar


def lizard(species, cx, cy, length, vertical=True):
    name = 'podarcis_orange' if species == 'liolepis' else 'podarcis_teal'
    ar = asset(name)[1]
    if vertical:   # head up
        image(name, cx - length / 2, cy - length * ar / 2, length, rot=-90, around=(cx, cy))
    else:
        image(name, cx - length / 2, cy - length * ar / 2, length)


def defs():
    o.append('<defs>')
    for gid, top, bot in (('grad_wc', ORANGE, WC), ('grad_ci', TEAL, CI), ('grad_cc', TEAL, CC)):
        o.append(f'<linearGradient id="{gid}" x1="0" y1="0" x2="0" y2="1">'
                 f'<stop offset="0" stop-color="{top}" stop-opacity=".25"/><stop offset="1" stop-color="{bot}" stop-opacity=".95"/>'
                 f'</linearGradient>')
    for name in sorted(used):
        b64, ar = asset(name)
        o.append(f'<image id="img_{name}" width="1" height="{ar:.5f}" preserveAspectRatio="none" href="data:image/png;base64,{b64}"/>')
    for mid, col in (('arr_o', ORANGE), ('arr_t', TEAL), ('arr_g', '#777777'), ('arr_p', CI)):
        o.append(f'<marker id="{mid}" viewBox="0 0 6 6" refX="0" refY="3" markerWidth="1.7" markerHeight="1.7" '
                 f'markerUnits="userSpaceOnUse" orient="auto"><path d="M0,0.3 L5.8,3 L0,5.7 z" fill="{col}"/></marker>')
    o.append('</defs>')


def panel_letter(letter, y):
    text(0.4, y, letter, size=4.2, weight='bold')


# ---------------------------------------------------------------------------
# geometry of panel a
# ---------------------------------------------------------------------------
PITCH, BW, BH = 5.8, 5.2, 14.0
X_WC = 10.0
X_CI = X_WC + 9 * PITCH + 6.0 - 0.6          # species gap
X_CC = X_CI + 9 * PITCH + 3.0 - 0.6          # group gap inside P. muralis
GROUPS = [('WC', X_WC, WC, WC_T, 'liolepis'), ('CI', X_CI, CI, CI_T, 'muralis'), ('CC', X_CC, CC, CC_T, 'muralis')]


def boxes(y, grouped):
    for key, x0, col, tint, sp in GROUPS:
        for i in range(9):
            x = x0 + i * PITCH
            if grouped:
                rect(x, y, BW, BH, fill=tint, stroke=col, sw=.5, rx=1.2)
            else:
                rect(x, y, BW, BH, fill='#ffffff', stroke=ORANGE if sp == 'liolepis' else TEAL, sw=.3, rx=1.2)
            lizard(sp, x + BW / 2, y + BH / 2, 11.2)


def box_cx(x0, i):
    return x0 + i * PITCH + BW / 2


def band(y, h, fill, label):
    rect(8.0, y, W - 8.0, h, fill=fill)
    text(4.6, y + h / 2, label, size=2.6, weight='bold', fill=MID, anchor='middle', rot=-90)


def build_a():
    panel_letter('a', 4.0)
    # --- Wild
    band(6, 21, GREY, 'Wild')
    lizard('liolepis', 24, 15.5, 22, vertical=False)
    text(38, 12.5, '*P. liolepis*', size=2.6, weight='bold', fill=ORANGE)
    text(38, 16.2, '9 individuals', size=2.2)
    text(38, 19.9, 'Warm and dry zone', size=2.2, fill=MID)
    xmid = (X_CI + X_CC + 8 * PITCH + BW) / 2            # middle of the CI + CC columns
    lizard('muralis', xmid - 12.0, 15.5, 22, vertical=False)
    text(xmid + 1.8, 12.5, '*P. muralis*', size=2.6, weight='bold', fill=TEAL)
    text(xmid + 1.8, 16.2, '18 individuals', size=2.2)
    text(xmid + 1.8, 19.9, 'Cold and wet zone', size=2.2, fill=MID)
    text(W - 1, 9.5, 'S1', size=2.2, anchor='end', fill=MID)

    # wide downward arrows: each wild group -> the experimental group(s) it supplied.
    # Gradient runs from the source species colour (top) to the destination group colour (bottom).
    def wide_arrow(x0, x1, y0, gid):
        xm = (x0 + x1) / 2
        sh, hd = 3.6, 5.2
        pts = [(x0 + 3.2, y0), (x1 - 3.2, y0), (x1 - 3.2, y0 + sh), (x1, y0 + sh), (xm, y0 + sh + hd),
               (x0, y0 + sh), (x0 + 3.2, y0 + sh)]
        d = 'M' + ' L'.join(f'{px:.2f},{py:.2f}' for px, py in pts) + ' Z'
        o.append(f'<path d="{d}" fill="url(#{gid})" stroke="none"/>')
    wide_arrow(X_WC, X_WC + 8 * PITCH + BW, 22.0, 'grad_wc')
    wide_arrow(X_CI, X_CI + 8 * PITCH + BW, 22.0, 'grad_ci')
    wide_arrow(X_CC, X_CC + 8 * PITCH + BW, 22.0, 'grad_cc')
    o.append('<g transform="translate(0 4)">')

    # --- Acclimation
    band(27, 25, '#ffffff', 'Acclimation')
    text(X_WC, 31.2, '28 °C / 40% RH', size=2.3, weight='bold', fill=ORANGE)
    text(X_CI, 31.2, '16 °C / 80% RH', size=2.3, weight='bold', fill=TEAL)
    text(X_CC, 31.2, '16 °C / 80% RH', size=2.3, weight='bold', fill=TEAL)
    text(W - 1, 31.2, 'S2 + respirometry', size=2.2, anchor='end', fill=MID)
    boxes(34, grouped=False)

    # --- Antibiotics + FMT
    band(52, 41, GREY, 'Antibiotics + FMT')
    text(X_WC, 56.2, '28 °C / 40% RH', size=2.3, weight='bold', fill=ORANGE)
    text(X_CI, 56.2, '16 °C / 80% RH', size=2.3, weight='bold', fill=TEAL)
    text(X_CC, 56.2, '16 °C / 80% RH', size=2.3, weight='bold', fill=TEAL)
    text(W - 1, 56.2, 'S3–S5', size=2.2, anchor='end', fill=MID)
    yb = 70.0
    boxes(yb, grouped=True)
    # interspecific donor arrows: P. liolepis donors -> cold-intervention recipients
    for i in range(9):
        xa, xb = box_cx(X_WC, i), box_cx(X_CI, i)
        path(f'M{xa:.2f},{yb - .6:.2f} C{xa:.2f},{yb - 11:.2f} {xb:.2f},{yb - 11:.2f} {xb:.2f},{yb - .6 - HL:.2f}',
             ORANGE, sw=.3, marker='arr_o', opacity=.85)
    # intraspecific rings below the boxes
    for key, x0, col, tint, sp in GROUPS:
        if key == 'CI':
            continue
        mk = 'arr_o' if sp == 'liolepis' else 'arr_t'
        for i in range(8):
            xa, xb = box_cx(x0, i), box_cx(x0, i + 1)
            path(f'M{xa:.2f},{yb + BH + .6:.2f} C{xa:.2f},{yb + BH + 3.6:.2f} {xb:.2f},{yb + BH + 3.6:.2f} {xb:.2f},{yb + BH + .6 + HL:.2f}',
                 col, sw=.3, marker=mk)
        xa, xb = box_cx(x0, 8), box_cx(x0, 0)
        path(f'M{xa:.2f},{yb + BH + .6:.2f} C{xa:.2f},{yb + BH + 8.4:.2f} {xb:.2f},{yb + BH + 8.4:.2f} {xb:.2f},{yb + BH + .6 + HL:.2f}',
             col, sw=.3, marker=mk)
    # arrow key (donor species)
    text(X_CI + 9 * PITCH / 2 - 2.2, yb - 11.7, '', size=2.2)

    # --- Post-FMT (band first, so the annotations below stay visible)
    band(93, 41, '#ffffff', 'Post-FMT')
    # gradual environmental change for the two cold-adapted groups
    xs = [x0 + 9 * PITCH / 2 - .3 for key, x0, col, tint, sp in GROUPS[1:]]
    for xm in xs:
        line(xm, 94.0, xm, 101.8 - HL, c='#777', sw=.35, marker='arr_g')
    text(sum(xs) / 2, 96.0, 'Gradual change (0.5 \u00b0C h\u207b\u00b9, 1% h\u207b\u00b9)', size=2.2,
         anchor='middle', weight='bold', fill=DARK)
    text(X_WC, 100.7, '28 °C / 40% RH', size=2.3, weight='bold', fill=ORANGE)
    text(X_CI, 100.7, '28 °C / 40% RH', size=2.3, weight='bold', fill=ORANGE)
    text(X_CC, 100.7, '28 °C / 40% RH', size=2.3, weight='bold', fill=ORANGE)
    text(W - 1, 100.7, 'S6–S7 + respirometry', size=2.2, anchor='end', fill=MID)
    boxes(103.2, grouped=True)
    # group brackets and labels
    yl = 103.2 + BH + 2.2
    names = {'WC': ('Warm-control (WC)', 'FMT from *P. liolepis* donors'),
             'CI': ('Cold-intervention (CI)', 'FMT from *P. liolepis* donors'),
             'CC': ('Cold-control (CC)', 'FMT from *P. muralis* donors')}
    for key, x0, col, tint, sp in GROUPS:
        xe = x0 + 8 * PITCH + BW
        rect(x0, yl, xe - x0, 1.1, fill=col, rx=.5)
        text((x0 + xe) / 2, yl + 4.9, f'{names[key][0]}', size=2.4, weight='bold', fill=col, anchor='middle')
        text((x0 + xe) / 2, yl + 8.3, f'n = 9 · {names[key][1]}', size=2.1, anchor='middle', fill=DARK)
    # key for arrows
    yk = 131.4
    path(f'M{X_WC:.2f},{yk:.2f} L{X_WC + 6 - HL:.2f},{yk:.2f}', ORANGE, sw=.3, marker='arr_o')
    text(X_WC + 7.4, yk + .8, 'FMT donor → recipient, *P. liolepis* donors', size=2.1)
    path(f'M{X_WC + 62:.2f},{yk:.2f} L{X_WC + 68 - HL:.2f},{yk:.2f}', TEAL, sw=.3, marker='arr_t')
    text(X_WC + 69.4, yk + .8, 'FMT donor → recipient, *P. muralis* donors', size=2.1)
    text(X_WC + 125, yk + .8, 'Frame colour = experimental group', size=2.1, fill=MID)
    o.append('</g>')


# ---------------------------------------------------------------------------
# panel b
# ---------------------------------------------------------------------------
def build_b(b0):
    panel_letter('b', b0 - 0.5)
    X0, WK = 12.0, 40.0
    ax_y, ax_h = b0 + 25.0, 3.6

    def tx(t):
        return X0 + t * WK

    # axis with weeks
    for k in range(4):
        rect(tx(k), ax_y, WK, ax_h, fill='#dcdcdc' if k % 2 == 0 else '#ececec')
        text(tx(k) + WK / 2, ax_y + 2.55, f'Week {k + 1}', size=2.2, anchor='middle', fill=MID)
    # phase bars
    bars = [(0.15, 1.0, 'Acclimation', '#8a8a8a'), (1.05, 1.3, '', '#5f5f5f'), (1.35, 1.6, '', '#a05a1c')]
    for t0, t1, lab, col in bars:
        rect(tx(t0), ax_y - 4.3, tx(t1) - tx(t0), 3.4, fill=col, rx=.6)
    text((tx(.15) + tx(1.0)) / 2, ax_y - 1.9, 'Acclimation', size=2.2, anchor='middle', fill='#ffffff', weight='bold')
    text(tx(1.05) + .3, ax_y - 6.0, 'Antibiotics', size=2.2, anchor='start', weight='bold', fill='#5f5f5f')
    text(tx(1.475), ax_y - 6.3, 'FMT', size=2.2, anchor='middle', weight='bold', fill='#a05a1c')

    # event markers (x position, sample index or None)
    S = {1: .0, 2: 1.0, 3: 1.3, 4: 1.4, 5: 1.6, 6: 2.7, 7: 3.7}
    for k, t in S.items():
        line(tx(t) + (.4 if k == 1 else 0), ax_y - 1.2 if k in (3, 4, 5) else ax_y - 8.0, tx(t) + (.4 if k == 1 else 0),
             ax_y, c='#b0b0b0', sw=.3)
    # labels above the axis
    text(X0 - .2, b0 + 4.0, 'Capture', size=2.2, weight='bold')
    text(X0 - .2, b0 + 7.2, 'Sample S1', size=2.2)
    line(X0 + .4, b0 + 8.6, X0 + .4, ax_y, c='#b0b0b0', sw=.3)
    text(tx(.12), b0 + 14.2, 'Deworming bath', size=2.2)
    line(tx(.12), b0 + 15.0, tx(.12), ax_y - 4.3, c='#b0b0b0', sw=.3)
    text(tx(1.0), b0 + 4.0, 'Respirometry', size=2.2, anchor='middle', weight='bold')
    text(tx(1.0), b0 + 7.2, 'Sample S2', size=2.2, anchor='middle')
    line(tx(1.0), b0 + 8.6, tx(1.0), ax_y, c='#b0b0b0', sw=.3)
    text(tx(1.7) - .2, b0 + 4.0, 'Change of conditions', size=2.2, weight='bold')
    text(tx(1.7) - .2, b0 + 7.2, 'to 28 °C / 40% RH', size=2.2)
    line(tx(1.7), b0 + 8.6, tx(1.7), ax_y, c='#b0b0b0', sw=.3)
    for k, nm in ((6, 'FMT1'), (7, 'FMT2')):
        text(tx(S[k]), b0 + 4.0, 'Respirometry', size=2.2, anchor='middle', weight='bold')
        text(tx(S[k]), b0 + 7.2, f'Sample S{k} ({nm})', size=2.2, anchor='middle')
        line(tx(S[k]), b0 + 8.6, tx(S[k]), ax_y, c='#b0b0b0', sw=.3)
    text(tx(3.95) + 1.5, b0 + 14.2, 'Euthanasia', size=2.2, anchor='middle')
    line(tx(3.95), b0 + 15.0, tx(3.95), ax_y, c='#b0b0b0', sw=.3)

    # sample row: seven equally spaced tubes with the names used in the Results
    names = {1: 'Wild', 2: 'Acclimation', 3: 'Antibiotics', 4: 'Transplant 1', 5: 'Transplant 2', 6: 'FMT1', 7: 'FMT2'}
    XS = {k: 21.0 + (k - 1) * 23.0 for k in range(1, 8)}
    ytube = ax_y + ax_h + 7.0
    for k in range(1, 8):
        # leader from the timeline to the sample icon
        t = S[k] + (0.01 if k == 1 else 0)
        line(tx(t) + (.4 if k == 1 else 0), ax_y + ax_h, XS[k], ytube - 1.0, c='#b0b0b0', sw=.25)
        th = image('eppendorf_tube', XS[k] - 6.0, ytube, 5.6)
        image('faecal_pellet', XS[k] - .2, ytube + th - 3.3, 5.8)
        text(XS[k], ytube + th + 3.2, f'S{k}', size=2.4, weight='bold', anchor='middle')
        text(XS[k], ytube + th + 6.4, names[k], size=2.2, anchor='middle')
    yresp = ytube + 7.8 + 11.5
    for k in (2, 6, 7):
        ch = image('respirometry_chamber', XS[k] - 7.2, yresp, 14.4)
        text(XS[k], yresp + ch + 3.4, 'Oxygen consumption', size=2.2, anchor='middle')
        line(XS[k], ytube + 7.8 + 8.0, XS[k], yresp - .8, c='#b0b0b0', sw=.3, dash='.8,.6')
    text(3.4, ytube + 6.5, 'Faecal samples', size=2.4, weight='bold', fill=MID, anchor='middle', rot=-90)
    text(3.4, yresp + 6.5, 'Respirometry', size=2.4, weight='bold', fill=MID, anchor='middle', rot=-90)
    return yresp + 14.4 * 0.399 + 7.0


def main():
    # panel a ends at ~134 mm; panel b starts below
    build_a()
    b0 = 144.0
    h = build_b(b0)
    H = h + 3.0
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}mm" height="{H:.1f}mm" viewBox="0 0 {W} {H:.2f}">',
           f'<rect width="{W}" height="{H:.2f}" fill="#ffffff"/>']
    defs_o, body = [], o[:]
    o.clear(); defs(); defs_o = o[:]
    out = '\n'.join(svg + defs_o + body + ['</svg>'])
    work = HERE
    (work / 'figure1.svg').write_text(out)
    FINAL.mkdir(exist_ok=True)
    env = dict(os.environ, GDK_BACKEND='x11')
    ink = lambda *a: subprocess.run([INK, *a], env=env, capture_output=True, text=True)
    ink(str(work / 'figure1.svg'), '--export-type=pdf', '--export-text-to-path', f'--export-filename={FINAL / "figure1.pdf"}')
    ink(str(work / 'figure1.svg'), '--export-type=png', '--export-dpi=300', '--export-background=white', f'--export-filename={FINAL / "figure1.png"}')
    subprocess.run(['magick', str(FINAL / 'figure1.png'), '-background', 'white', '-alpha', 'remove', '-alpha', 'off',
                    '-compress', 'lzw', str(FINAL / 'figure1.tiff')], check=True)
    ink(str(work / 'figure1.svg'), '--export-type=png', '--export-dpi=110', '--export-background=white', f'--export-filename={work / "figure1_preview.png"}')
    print(f'WROTE figure1.svg ({W} x {H:.1f} mm), PDF, PNG 300 dpi, TIFF, preview')


if __name__ == '__main__':
    main()
