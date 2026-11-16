#!/usr/bin/env python3
"""Manuscript Figure 3 (alpha diversity + dbRDA centroid trajectories).

Run from the repository root:
    python3 -I manuscript_figures/figure3/build_figure3.py

Draws both panels in R (make_panels.R, vector SVG at final size), composes them with the
figures skill's compose.py (180 mm wide, panel letters a/b) and exports an outlined PDF,
a 300 dpi PNG and a TIFF to manuscript_figures/final/.
"""
import os, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FINAL = HERE.parent / 'final'
INK = '/Applications/Inkscape.app/Contents/MacOS/inkscape'
COMPOSE = Path('/Users/anttonalberdi/My Drive/Supervision/Garazi PhD/Manuscripts/Lizards_adaptation_2024/.claude/skills/figures/compose.py')

def run(*a, **k):
    r = subprocess.run(a, capture_output=True, text=True, cwd=k.get('cwd'), env=k.get('env'))
    if r.returncode:
        sys.exit(f'FAILED: {a[:3]}\n{r.stdout}\n{r.stderr}')
    return r

run('Rscript', str(HERE / 'make_panels.R'), cwd=ROOT)
r = run(sys.executable, '-I', str(COMPOSE), str(HERE / 'figure3.svg'), '--width', '180', '--cols', '1',
        '--gap', '2', '--label-size', '9', str(HERE / 'fig3a.svg'), str(HERE / 'fig3b.svg'))
print(r.stdout.strip())
FINAL.mkdir(exist_ok=True)
env = dict(os.environ, GDK_BACKEND='x11')
ink = lambda *a: subprocess.run([INK, *a], env=env, capture_output=True, text=True)
src = str(HERE / 'figure3.svg')
ink(src, '--export-type=pdf', '--export-text-to-path', f'--export-filename={FINAL / "figure3.pdf"}')
ink(src, '--export-type=png', '--export-dpi=300', '--export-background=white', f'--export-filename={FINAL / "figure3.png"}')
run('magick', str(FINAL / 'figure3.png'), '-background', 'white', '-alpha', 'remove', '-alpha', 'off', '-compress', 'lzw', str(FINAL / 'figure3.tiff'))
ink(src, '--export-type=png', '--export-dpi=110', '--export-background=white', f'--export-filename={HERE / "figure3_preview.png"}')
print('WROTE figure3.svg, PDF, PNG 300 dpi, TIFF, preview')
