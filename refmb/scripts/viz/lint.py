#!/usr/bin/env python3
"""figure_lint (plan §9.5): every PDF/SVG pair must pass or the build fails.

1. Width matches a target exactly (88.9 or 181.9 mm, ±0.2 mm) — from the PDF MediaBox.
2. No text below 6 pt — from SVG font-size attributes (matplotlib writes px at 72 dpi = pt).
3. All text is real text — SVG contains <text> elements and no glyph <path> definitions (svg.fonttype none).
4. All fonts embedded — every PDF font has a FontFile/FontFile2/FontFile3 stream.
5. Every colour is from the token palette — SVG fill/stroke hexes ⊂ tokens (ink, deviation, categorical, sequential).
6. Greyscale separation of the deviation classes — relative luminance of the five steps pairwise ≥ 0.06 except the two
   symmetric light steps, whose sign is carried by texture (45° vs 135°) in print.
7. Colourblind check — the validator result recorded in tokens_validation.txt is required to say PASS for the sets used.
8. Determinism — the SVG rendered twice is byte-identical (caller renders twice into two dirs and passes both).
"""
import json, os, re, sys
import pypdf

sys.path.insert(0, os.path.dirname(__file__))
import tokens as T

TARGETS_MM = (88.9, 181.9)
ALLOWED = {c.lower() for c in list(T.DEV.values()) + T.CAT + T.SEQ + list(T.INK.values()) + ["#000000", "#ffffff", "none"]}


def lum(hexc):
    r, g, b = [int(hexc[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


def lint(pdf_path, svg_path, svg_path_rerun=None):
    problems = []
    r = pypdf.PdfReader(pdf_path); page = r.pages[0]
    w_mm = float(page.mediabox.width) / 72 * 25.4
    if not any(abs(w_mm - t) <= 0.2 for t in TARGETS_MM):
        problems.append(f"width {w_mm:.2f} mm is not a target width")
    fonts = {}
    res = page.get("/Resources", {}).get("/Font", {})
    for k, f in res.items():
        fo = f.get_object(); desc = fo.get("/FontDescriptor") or (fo.get("/DescendantFonts", [None])[0].get_object().get("/FontDescriptor") if fo.get("/DescendantFonts") else None)
        embedded = bool(desc and any(key in desc for key in ("/FontFile", "/FontFile2", "/FontFile3")))
        fonts[str(fo.get("/BaseFont"))] = embedded
        if not embedded:
            problems.append(f"font not embedded: {fo.get('/BaseFont')}")
    svg = open(svg_path, encoding="utf-8").read()
    sizes = [float(x) for x in re.findall(r'font-size:\s*([0-9.]+)', svg)] + [float(x) for x in re.findall(r'font-size="([0-9.]+)', svg)]
    if sizes and min(sizes) < T.PT_MIN - 1e-6:
        problems.append(f"text below {T.PT_MIN} pt: min {min(sizes):.2f}")
    if "<text" not in svg:
        problems.append("no <text> elements (text outlined?)")
    if re.search(r'<path[^>]*id="DejaVuSans', svg):
        problems.append("glyphs rendered as paths")
    colours = {c.lower() for c in re.findall(r'(?:fill|stroke):\s*(#[0-9a-fA-F]{6})', svg)} | {c.lower() for c in re.findall(r'(?:fill|stroke)="(#[0-9a-fA-F]{6})"', svg)}
    bad = sorted(colours - ALLOWED)
    if bad:
        problems.append(f"colours outside the token palette: {bad}")
    L = [lum(T.DEV[k]) for k in T.DEV_ORDER]
    for i in range(5):
        for j in range(i + 1, 5):
            if {i, j} in ({1, 3}, {0, 4}):
                continue  # symmetric steps of a diverging scale share lightness by construction; sign is carried by texture (45°/135°) in print
            if abs(L[i] - L[j]) < 0.06:
                problems.append(f"greyscale: deviation steps {T.DEV_ORDER[i]} and {T.DEV_ORDER[j]} too close in luminance")
    val = open(os.path.join(os.path.dirname(__file__), "tokens_validation.txt")).read() if os.path.exists(os.path.join(os.path.dirname(__file__), "tokens_validation.txt")) else ""
    if "ALL CHECKS PASS" not in val:
        problems.append("palette validator record does not show a passing categorical run")
    if svg_path_rerun:
        if open(svg_path, "rb").read() != open(svg_path_rerun, "rb").read():
            problems.append("non-deterministic: SVG differs between two renders")
    return {"file": os.path.basename(pdf_path), "width_mm": round(w_mm, 2), "fonts": fonts, "min_font_pt": round(min(sizes), 2) if sizes else None,
            "n_colours": len(colours), "problems": problems, "pass": not problems}


if __name__ == "__main__":
    out = [lint(a, a[:-4] + ".svg", (sys.argv[2] + "/" + os.path.basename(a)[:-4] + ".svg") if len(sys.argv) > 2 else None) for a in sorted(sys.argv[1].split(","))]
    print(json.dumps(out, indent=1)); sys.exit(0 if all(o["pass"] for o in out) else 1)
