#!/usr/bin/env python3
"""Manuscript Figure 4 (oxygen consumption), same style as Figure 3a.

Run from the repository root:
    python3 -I manuscript_figures/figure4/build_figure4.py

Draws the panel in R (vector SVG at final size, 180 mm wide) and exports an outlined PDF,
a 300 dpi PNG and a TIFF to manuscript_figures/final/.
"""
import os, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FINAL = HERE.parent / 'final'
INK = '/Applications/Inkscape.app/Contents/MacOS/inkscape'

def run(*a, **k):
    r = subprocess.run(a, capture_output=True, text=True, cwd=k.get('cwd'))
    if r.returncode:
        sys.exit(f'FAILED: {a[:3]}\n{r.stdout}\n{r.stderr}')
    return r

print(run('Rscript', str(HERE / 'make_figure4.R'), cwd=ROOT).stdout.strip())
FINAL.mkdir(exist_ok=True)
env = dict(os.environ, GDK_BACKEND='x11')
ink = lambda *a: subprocess.run([INK, *a], env=env, capture_output=True, text=True)
src = str(HERE / 'figure4.svg')
ink(src, '--export-type=pdf', '--export-text-to-path', f'--export-filename={FINAL / "figure4.pdf"}')
ink(src, '--export-type=png', '--export-dpi=300', '--export-background=white', f'--export-filename={FINAL / "figure4.png"}')
run('magick', str(FINAL / 'figure4.png'), '-background', 'white', '-alpha', 'remove', '-alpha', 'off', '-compress', 'lzw', str(FINAL / 'figure4.tiff'))
ink(src, '--export-type=png', '--export-dpi=110', '--export-background=white', f'--export-filename={HERE / "figure4_preview.png"}')
print('WROTE figure4.svg, PDF, PNG 300 dpi, TIFF, preview')
