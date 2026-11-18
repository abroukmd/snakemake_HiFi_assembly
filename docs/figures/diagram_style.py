"""Shared vector drawing primitives for the release workflow schematics."""

from html import escape
import math
from pathlib import Path


C = {
    "paper": "#FFFFFF", "ink": "#192A35", "muted": "#5C6870",
    "line": "#A7ADAE", "input": "#2B9FF3", "hif": "#67BA48",
    "lja": "#F1A43D", "mga": "#A365D3", "mito": "#17ABA5",
    "qc": "#9DD6AB", "bypass": "#EC6869",
}


def tint(color, fraction=0.12):
    value = C.get(color, color).lstrip("#")
    channels = [int(value[i:i + 2], 16) for i in [0, 2, 4]]
    return "#" + "".join(f"{round(c * fraction + 255 * (1 - fraction)):02x}" for c in channels)


class Diagram:
    def __init__(self, width, height, title, description):
        self.width, self.height = width, height
        self.parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img" aria-labelledby="figure-title figure-description">',
            f'<title id="figure-title">{escape(title)}</title>',
            f'<desc id="figure-description">{escape(description)}</desc>',
            '<defs>',
        ]
        for key, value in C.items():
            self.parts.append(
                f'<marker id="arrow-{key}" viewBox="0 0 14 14" refX="12" refY="7" '
                'markerWidth="17" markerHeight="17" orient="auto-start-reverse" '
                f'markerUnits="userSpaceOnUse"><path d="M 1 1 L 12 7 L 1 13 Z" fill="{value}"/></marker>'
            )
        self.parts.extend(['</defs>', '<g font-family="DejaVu Sans, Arial, sans-serif">'])
        self.rect(0, 0, width, height, "paper", radius=0)
        self._font = None
        try:
            from PIL import ImageFont
            font_root = Path("/usr/share/fonts/truetype/dejavu")
            if all((font_root / name).is_file() for name in ["DejaVuSans.ttf", "DejaVuSans-Bold.ttf"]):
                ImageFont.truetype(str(font_root / "DejaVuSans.ttf"), 17)
                self._font = ImageFont
        except ImportError:
            pass

    def measure(self, text, size=20, bold=False):
        filename = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
        path = Path("/usr/share/fonts/truetype/dejavu") / filename
        if self._font and path.is_file():
            return self._font.truetype(str(path), size).getlength(text)
        return len(text) * size * (0.64 if bold else 0.55)

    def rect(self, x, y, w, h, fill="paper", stroke=None, radius=16, dash=False, width=2.5):
        attributes = f'fill="{C.get(fill, fill)}"'
        if stroke:
            attributes += f' stroke="{C.get(stroke, stroke)}" stroke-width="{width}"'
        if dash:
            attributes += ' stroke-dasharray="8 6"'
        self.parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" {attributes}/>')

    def text(self, x, y, text, size=20, color="ink", bold=False, anchor="start", max_width=None):
        if max_width and self._font and self.measure(text, size, bold) > max_width:
            raise ValueError(f"Text exceeds available width {max_width}: {text}")
        weight = ' font-weight="700"' if bold else ""
        self.parts.append(f'<text x="{x}" y="{y}" font-size="{size}" fill="{C.get(color, color)}" '
                          f'text-anchor="{anchor}"{weight}>{escape(text)}</text>')

    def label(self, x, y, text, color="muted", size=17, anchor="middle"):
        # Keep the SVG geometry identical with or without optional Pillow.
        # Glyph measurements are a validation aid, not a layout dependency.
        width = round(len(text) * size * 0.60, 3)
        left = x - width / 2 if anchor == "middle" else x - width if anchor == "end" else x
        self.rect(left - 6, y - size, width + 12, size + 8, "paper", radius=3)
        self.text(x, y, text, size, color, anchor=anchor)

    def node(self, x, y, w, h, title, subtitle=None, color=None, optional=False, output=False, size=22):
        stroke = color if output or optional else "line"
        fill = tint(color or "qc", 0.07) if output else "paper"
        self.rect(x, y, w, h, fill, stroke or "line", radius=17, dash=optional, width=3)
        lines = str(title).split("\n")
        details = [subtitle] if isinstance(subtitle, str) else list(subtitle or [])
        offsets = [i * (size + 4) for i in range(len(lines))]
        for i in range(len(details)):
            offsets.append((len(lines) - 1) * (size + 4) + (size + 17) / 2 + 5 + i * 22)
        last_size = 17 if details else size
        # Centre the visible ink block, rather than adding a whole font-size to
        # a line-height block. This keeps subtitles away from the lower border.
        baseline = y + h / 2 - (offsets[-1] - size * 0.78 + last_size * 0.23) / 2
        for offset, line in zip(offsets[:len(lines)], lines):
            self.text(x + w / 2, baseline + offset, line, size, "ink", anchor="middle", max_width=w - 24)
        for offset, line in zip(offsets[len(lines):], details):
            self.text(x + w / 2, baseline + offset, line, 17, "muted", anchor="middle", max_width=w - 20)

    def input(self, cx, cy, rx, ry, lines, optional=False, color="input", size=19):
        dashed = ' stroke-dasharray="5 5"' if optional else ""
        self.parts.append(f'<ellipse cx="{cx}" cy="{cy}" rx="{rx}" ry="{ry}" fill="white" '
                          f'stroke="{C.get(color, color)}" stroke-width="3.5"{dashed}/>')
        lines = [lines] if isinstance(lines, str) else list(lines)
        y = cy - (len(lines) - 1) * (size + 4) / 2 + size * 0.34
        for i, line in enumerate(lines):
            self.text(cx, y + i * (size + 4), line, size, anchor="middle", max_width=rx * 1.80)

    def route(self, points, color="hif", dash=False, arrow=True, width=5):
        points = [(float(x), float(y)) for x, y in points]
        points = [p for i, p in enumerate(points) if i == 0 or p != points[i - 1]]
        if len(points) < 2:
            raise ValueError("A route needs two distinct points")
        data = f"M {points[0][0]:g} {points[0][1]:g}"
        for i in range(1, len(points) - 1):
            a, b, c = points[i - 1], points[i], points[i + 1]
            before, after = math.dist(a, b), math.dist(b, c)
            radius = min(13, before / 2, after / 2)
            p = (b[0] - (b[0] - a[0]) * radius / before,
                 b[1] - (b[1] - a[1]) * radius / before)
            q = (b[0] + (c[0] - b[0]) * radius / after,
                 b[1] + (c[1] - b[1]) * radius / after)
            data += f" L {p[0]:g} {p[1]:g} Q {b[0]:g} {b[1]:g} {q[0]:g} {q[1]:g}"
        data += f" L {points[-1][0]:g} {points[-1][1]:g}"
        self.parts.append(f'<path d="{data}" fill="none" stroke="white" stroke-width="{width + 5}" '
                          'stroke-linecap="round" stroke-linejoin="round"/>')
        attributes = f'fill="none" stroke="{C.get(color, color)}" stroke-width="{width}" stroke-linecap="round" stroke-linejoin="round"'
        if dash:
            attributes += ' stroke-dasharray="9 7"'
        if arrow:
            if color not in C:
                raise ValueError("Arrow colors must be palette names")
            attributes += f' marker-end="url(#arrow-{color})"'
        self.parts.append(f'<path d="{data}" {attributes}/>')

    def group(self, x, y, w, h, title, color="qc"):
        self.rect(x, y, w, h, tint(color, 0.045), tint(color, 0.45), radius=22, dash=True, width=1.6)
        self.text(x + 20, y + 29, title, 19, "muted", bold=True, max_width=w - 40)

    def checkpoint(self, cx, cy, label="QC"):
        self.parts.append(f'<circle cx="{cx}" cy="{cy}" r="16" fill="{tint("qc", 0.25)}" '
                          f'stroke="{C["qc"]}" stroke-width="2.5"/>')
        self.text(cx, cy + 5, label, 12, "ink", True, "middle")

    def legend(self, x, y, items, title="Paths and settings", gap=60):
        self.text(x, y, title, 21, bold=True)
        for i, item in enumerate(items):
            color, label, detail, *dashed = item
            yy = y + 38 + gap * i
            self.route([(x, yy - 5), (x + 55, yy - 5)], color, dash=bool(dashed and dashed[0]), arrow=False, width=5)
            self.text(x + 72, yy, label, 19, bold=True)
            if detail:
                self.text(x + 72, yy + 23, detail, 16, "muted")

    def motif(self, cx, cy, scale=1):
        """An original multicolour DNA ribbon, drawn entirely as vector paths."""
        self.parts.append(f'<g transform="translate({cx} {cy}) scale({scale})" aria-label="Original HiFi DNA ribbon motif">')
        colors = ["hif", "lja", "mga", "input", "mito"]
        for i in range(18):
            x = -165 + i * 330 / 17
            wave = 50 * math.sin((x + 170) * math.pi / 170)
            if abs(wave) < 10:
                continue
            self.parts.append(f'<path d="M {x:g} {-wave:g} L {x:g} {wave:g}" stroke="{C["ink"]}" '
                              'stroke-width="9" stroke-linecap="round"/>')
            self.parts.append(f'<path d="M {x:g} {-wave:g} L {x:g} {wave:g}" stroke="{C[colors[i % len(colors)]]}" '
                              'stroke-width="5" stroke-linecap="round"/>')
        for sign, base in [(1, "input"), (-1, "mga")]:
            points = [(x, sign * 50 * math.sin((x + 170) * math.pi / 170)) for x in range(-170, 171, 3)]
            data = "M " + " L ".join(f"{x:g} {y:g}" for x, y in points)
            self.parts.append(f'<path d="{data}" fill="none" stroke="{C["ink"]}" stroke-width="18" '
                              'stroke-linejoin="round" stroke-linecap="round"/>')
            self.parts.append(f'<path d="{data}" fill="none" stroke="{C[base]}" stroke-width="11" '
                              'stroke-linejoin="round" stroke-linecap="round"/>')
            for j, color in enumerate(["hif", "lja", "mito"]):
                self.parts.append(f'<path d="{data}" fill="none" stroke="{C[color]}" stroke-width="11" '
                                  f'stroke-dasharray="25 105" stroke-dashoffset="{-j * 42}" stroke-linecap="butt"/>')
        self.text(0, 124, "HiFi", 69, "ink", True, "middle")
        self.text(0, 156, "GENOME ASSEMBLY", 18, "muted", anchor="middle")
        self.parts.append('</g>')

    def save(self, path):
        self.parts.extend(['</g>', '</svg>'])
        Path(path).write_text("\n".join(self.parts) + "\n", encoding="utf-8")
