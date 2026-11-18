#!/usr/bin/env python3
"""Regenerate the curated workflow schematics as SVG, vector PDF and PNG.

The drawing modules use Python's standard library. Optional Pillow/DejaVu Sans
checks text widths. PDF/PNG export requires rsvg-convert from librsvg.
"""
import argparse
from pathlib import Path
import shutil
import subprocess
from assembly_diagram import draw as draw_assembly
from qc_diagram import draw as draw_qc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--svg-only', action='store_true')
    parser.add_argument('--png-width', type=int, default=3200)
    args = parser.parse_args()
    if args.png_width <= 0:
        parser.error('--png-width must be positive')
    renderer = shutil.which('rsvg-convert')
    if not args.svg_only and not renderer:
        parser.error('PDF/PNG export needs rsvg-convert; use --svg-only for SVG output')
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name, draw in [('workflow_overview', draw_assembly), ('quality_control', draw_qc)]:
        svg = args.outdir / (name + '.svg')
        draw(svg)
        if not args.svg_only:
            for format in ['pdf', 'png']:
                command = [renderer, '--format', format, '--output', str(svg.with_suffix('.' + format))]
                if format == 'png':
                    command.extend(['--width', str(args.png_width)])
                subprocess.run(command + [str(svg)], check=True)
        print('Created ' + name + ': SVG' + ('' if args.svg_only else ', PDF, PNG'))


if __name__ == '__main__':
    main()
