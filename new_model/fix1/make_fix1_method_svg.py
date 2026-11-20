# Drawing of fix 1 (informative k-mers), same style as final_script/figures/make_kmer_method_svg.py.
# Numbers: pilot (new_model/pilot_fix1/README.md) and the full build (new_model/fix1/README.md).
# Usage (any python3): python new_model/fix1/make_fix1_method_svg.py  ->  new_model/fix1/fix1_method.svg
import os
W, Hh = 1100, 1100
B73 = "#3A6FB5"; SEG = "#B8860B"; DON = "#BB5566"; OK = "#2E8B57"; GRAY = "#9AA3AD"; INK = "#1F2933"; SOFT = "#5B6670"
LIGHT = "#E6EAF0"
o = []
def t(x, y, s, size=15, w="normal", fill=INK, anchor="start"):
    o.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{w}" fill="{fill}" text-anchor="{anchor}">{s}</text>')
def rect(x, y, w, h, fill, rx=0, stroke="none", sw=0, extra=""):
    o.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" {extra}/>')
def line(x1, y1, x2, y2, c, w=2, extra=""):
    o.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{c}" stroke-width="{w}" {extra}/>')
def circ(x, y, r, fill): o.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}"/>')
def panel(y, num, title, sub):
    circ(44, y + 14, 14, INK); t(44, y + 19, str(num), 14, "bold", "#fff", "middle")
    t(68, y + 19, title, 18, "bold"); t(68, y + 40, sub, 13.5, fill=SOFT)
def stick(x1, x2, y, c, w=4):
    line(x1, y, x2, y, c, w, 'stroke-linecap="round"'); circ(x1, y, 4.5, c); circ(x2, y, 4.5, c)
def check(x, y): o.append(f'<path d="M{x-6},{y} l4,5 l9,-10" fill="none" stroke="{OK}" stroke-width="3" stroke-linecap="round"/>')
def cross(x, y): o.append(f'<path d="M{x-6},{y-6} l12,12 M{x+6},{y-6} l-12,12" stroke="{DON}" stroke-width="3" stroke-linecap="round"/>')

def bars(x0, base, h, vals, errs, title, contrast):
    """three bars B73 / SEG / DONOR with a +-noise whisker, and the contrast"""
    t(x0 + 165, base - h - 34, title, 14, "bold", anchor="middle")
    line(x0 - 10, base, x0 + 340, base, SOFT, 1.5)
    line(x0 - 10, base - h, x0 + 340, base - h, SOFT, 1, 'stroke-dasharray="5,5"')
    t(x0 - 16, base - h + 5, "1.0", 12, fill=SOFT, anchor="end"); t(x0 - 16, base + 5, "0", 12, fill=SOFT, anchor="end")
    for i, (lab, v, e, c) in enumerate(zip(["B73", "SEG", "DONOR"], vals, errs, [B73, SEG, DON])):
        x = x0 + 15 + i * 110
        rect(x, base - h * v, 70, h * v, c, 4)
        cx = x + 35; line(cx, base - h * (v + e), cx, base - h * (v - e), INK, 1.6)
        line(cx - 8, base - h * (v + e), cx + 8, base - h * (v + e), INK, 1.6); line(cx - 8, base - h * (v - e), cx + 8, base - h * (v - e), INK, 1.6)
        t(x + 35, base + 20, lab, 13.5, "bold", fill=c, anchor="middle")
        t(x + 35, base - h * (v + e) - 7, f"{v:.2f}", 12.5, "bold", anchor="middle")
    t(x0 + 165, base + 46, contrast, 13.5, "bold", fill=OK if "3.2" in contrast else DON, anchor="middle")

o.append(f'<rect x="0" y="0" width="{W}" height="{Hh}" fill="#FFFFFF"/>')
t(W / 2, 38, "Fix 1: count only the informative k-mers, where they help", 22, "bold", anchor="middle")
t(W / 2, 62, "folder new_model/fix1 · same pools, same HMM · only the k-mers counted in some windows change", 14, fill=SOFT, anchor="middle")

# ---------- panel 1: the problem ----------
y = 90
panel(y, 1, "The problem: windows where the teosinte keeps most B73 k-mers",
      "e.g. parviglumis chr7 54-65 Mb: the teosinte has 90% of the B73 single-copy k-mers. Those shared k-mers are found in every pool,")
t(68, y + 58, "whatever its state, so they add counting noise and no signal. DONOR and SEG end up less than one noise unit apart.", 13.5, fill=SOFT)
gy = y + 125
rect(80, gy - 6, 470, 12, LIGHT, 6); t(70, gy + 5, "B73", 14, "bold", anchor="end")
xs = [90 + i * 46 for i in range(10)]
for i, x in enumerate(xs):
    inf = i in (6,)
    stick(x, x + 36, gy - 22, DON if inf else GRAY)
stick(600, 640, gy - 30, GRAY); t(650, gy - 25, "shared with the teosinte: found in B73, SEG and DONOR pools", 13)
stick(600, 640, gy - 4, DON); t(650, gy + 1, "informative: missing in the teosinte, found only from B73 DNA", 13)

# ---------- panel 2: informative k-mers ----------
y = 250
panel(y, 2, "Keep only the informative k-mers (fx1_1, fx1_2)",
      "Informative k-mers of a taxon = B73 single-copy k-mers found in NONE of that taxon's assemblies (KMC: B73 list minus each assembly).")
t(68, y + 58, "Every pool is recounted against its taxon's list (as step 17); the 10 B73 controls are counted against all five lists.", 13.5, fill=SOFT)
gy = y + 100
labels = [("B73", None), ("teosinte assembly 1", {0, 1, 2, 3, 4, 5, 7, 8, 9}), ("teosinte assembly 2", {0, 1, 2, 3, 4, 5, 7, 8, 9})]
for r_, (lab, has) in enumerate(labels):
    yy = gy + r_ * 34
    t(200, yy + 5, lab, 13, "bold" if has is None else "normal", fill=INK if has is None else SOFT, anchor="end")
    for i in range(10):
        x = 215 + i * 46
        if has is None: stick(x, x + 36, yy, DON if i == 6 else GRAY)
        else:
            (check if i in has else cross)(x + 18, yy)
o.append(f'<path d="M700,{gy+30} l40,0" stroke="{INK}" stroke-width="2.5" marker-end="url(#arr)"/>')
t(760, gy + 22, "keep the k-mers that every", 13.5); t(760, gy + 40, "assembly lacks (here 1 of 10)", 13.5)
stick(770, 810, gy + 66, DON); t(822, gy + 71, "informative only", 13.5, "bold", fill=DON)
t(80, gy + 108, "Fewer k-mers per window (more counting noise), but every one of them carries signal: DONOR pools lose most of them, SEG pools about half.", 13.5)

# ---------- panel 3: score again ----------
y = 470
panel(y, 3, "Score the window again with the informative k-mers",
      "Same normalisation as step 20 (pool ÷ B73 controls). Our donors still carry some of these k-mers, so their DONOR level is learned (panel 4).")
base = y + 205
bars(110, base, 110, [1.0, 0.95, 0.90], [0.10, 0.10, 0.10], "ALL B73 k-mers (current)", "contrast 0.5")
bars(620, base, 110, [1.0, 0.70, 0.35], [0.11, 0.11, 0.11], "INFORMATIVE k-mers only (fix 1)", "contrast 3.2")
t(W / 2, base + 76, "contrast = gap between SEG and DONOR ÷ noise of the window (whiskers = ±1 noise unit). Contrasts: parviglumis chr7 54-65 Mb, pilot; bar heights illustrative", 12.5, fill=SOFT, anchor="middle")

# ---------- panel 4: choose window by window, then HMM ----------
y = 790
panel(y, 4, "Choose window by window, then run the same HMM (fx1_3 to fx1_6)",
      "A window switches if B73 has ≥ 200 informative k-mers there AND the new contrast is ≥ 1.2 × the old one (measured from all pools of the taxon).")
t(68, y + 58, "Switched windows take their DONOR level and its spread from our pools called DONOR there (the assemblies have none of these k-mers).", 13.5, fill=SOFT)
sy = y + 92
import random
random.seed(3)
x = 80; t(70, sy + 13, "windows", 12.5, fill=SOFT, anchor="end")
for i in range(94):
    hot = 52 <= i <= 66
    sw = hot or random.random() < 0.38
    rect(x, sy, 8, 18, DON if sw else LIGHT); x += 10
t(80 + 59 * 10, sy - 6, "hotspot: ~all switched", 12, "bold", fill=DON, anchor="middle")
rect(80, sy + 30, 14, 12, DON); t(100, sy + 41, "informative k-mers (33-53% of the measured windows, by taxon)", 12.5)
rect(520, sy + 30, 14, 12, LIGHT); t(540, sy + 41, "all B73 k-mers, assembly level (as before)", 12.5)
# example calls
cy = sy + 72
def track(y_, lab, segs):
    t(250, y_ + 13, lab, 12.5, "bold", anchor="end")
    for a, b, c in segs: rect(260 + a, y_, b - a, 18, c)
t(80, cy - 10, "one pool, chr7 (schematic)", 12.5, fill=SOFT)
track(cy, "current calls", [(0, 120, B73), (120, 330, DON), (330, 470, SEG), (470, 690, DON), (690, 760, B73)])
track(cy + 26, "with fix 1", [(0, 120, B73), (120, 690, DON), (690, 760, B73)])
track(cy + 52, "targeted sequencing", [(0, 120, B73), (120, 690, DON), (690, 760, B73)])
t(260 + 400, cy + 92, "the SEG island inside DONOR disappears (chr7 54-65 Mb: 28 → 7 pools)", 13, "bold", fill=OK, anchor="middle")
t(80, cy + 120, "Targeted sequencing (85 plots): SEG tracts that are AB 74% → 86% (AA, i.e. no teosinte, 18% → 6%); DONOR tracts that are BB 91% → 96%;", 13)
t(80, cy + 140, "B73 tracts that are AA 97% → 97%. Not fixed: chr6 21-23 and chr1 297 Mb (a wrong DONOR level, not low contrast) and parviglumis chr3 75-78 Mb.", 13)
o.insert(0, f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{Hh}" viewBox="0 0 {W} {Hh}" font-family="Helvetica, Arial, sans-serif">'
            f'<defs><marker id="arr" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="{INK}"/></marker></defs>')
o.append("</svg>")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fix1_method.svg")
open(out, "w").write("".join(o)); print("wrote", out)
