"""Slide version of the CU1 level diagram: add the mirror image of the laser arrow, mF=0 -> mF'=-1.

    python scripts/add_polswitch_arrow_cu1_level_diagram.py

With polarization switching (X/Z) the laser drives both mF' = +1 and mF' = -1; the diagram from
scripts/cleanup_level_diagrams.py only draws the arrow to +|mF'|. This copies that arrow's pixels, mirrored about the
mF = 0 column, so the new arrow is identical in style. Nothing else is changed (no text added).
Input: the top 1290 px of docs/state_cleanup/figures/slides/CU1_level_diagram_R1_F2_150Vcm_cleanup_19MHz.png
(the part shown on the slide). Writes ..._19MHz_polswitch.png next to it.
"""
from pathlib import Path

import numpy as np
from PIL import Image

FIG = Path(__file__).resolve().parents[1] / "docs" / "state_cleanup" / "figures" / "slides"
SRC = FIG / "CU1_level_diagram_R1_F2_150Vcm_cleanup_19MHz.png"
OUT = FIG / "CU1_level_diagram_R1_F2_150Vcm_cleanup_19MHz_polswitch.png"
CROP_H = 1290
FOOT, TIP = (842.5, 821.0), (1001.5, 484.0)   # existing arrow, pixels
AXIS2 = 1686                                  # 2 x mirror axis (mF = 0 column): box edges 640/730 <-> 956/1046
BAND = 12                                     # half-width of the region around the arrow that is copied

im = Image.open(SRC).convert("RGB")
a = np.asarray(im)[:CROP_H].copy()
yy, xx = np.mgrid[:a.shape[0], :a.shape[1]]
(x0, y0), (x1, y1) = FOOT, TIP
t = np.clip(((xx - x0) * (x1 - x0) + (yy - y0) * (y1 - y0)) / ((x1 - x0) ** 2 + (y1 - y0) ** 2), 0, 1)
near = np.hypot(xx - (x0 + t * (x1 - x0)), yy - (y0 + t * (y1 - y0))) <= BAND
px = a.astype(int)
bluish = (px[..., 2] - px[..., 0] > 30) & (px[..., 2] - px[..., 1] > 20)
ys, xs = np.nonzero(near & bluish)
a[ys, AXIS2 - xs] = a[ys, xs]
Image.fromarray(a).save(OUT)
print("wrote", OUT, a.shape, "| pixels copied:", len(xs))
