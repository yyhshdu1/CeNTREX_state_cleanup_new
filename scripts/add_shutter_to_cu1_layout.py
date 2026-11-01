"""Crop the CU1 part of the user's CU1 + CU2 optical-layout drawing, add an SRS shutter after the SHG cavity and a
photodiode at the end of the multipass (it monitors the shutter state).

    python scripts/add_shutter_to_cu1_layout.py <layout.png>

Nothing in the drawing is changed except the added shutter (with label) and photodiode (symbol only). The crop is upscaled 2x (Lanczos) so
that the new label is sharp. Writes docs/state_cleanup/figures/slides/cu1_optical_layout_shutter.png.
"""
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "state_cleanup" / "figures" / "slides" / "cu1_optical_layout_shutter.png"
CROP = (0, 40, 1078, 225)          # CU1 row of the drawing (pixels of the original)
BEAM_Y, SHUTTER_X = 116, 185       # beam row after the SHG, shutter position (original pixels)
PD_X, PD_Y = 1024, 58              # end of the multipass output beam: photodiode (original pixels)
S = 2                               # upscale factor

src = Image.open(sys.argv[1]).convert("RGB")
im = src.crop(CROP)
im = im.resize((im.width * S, im.height * S), Image.LANCZOS)
d = ImageDraw.Draw(im)
x, y = SHUTTER_X * S, (BEAM_Y - CROP[1]) * S

# shutter head: dark housing across the beam with a light blade
w, h = 26 * S, 34 * S
d.rounded_rectangle((x - w // 2, y - h // 2, x + w // 2, y + h // 2), radius=4 * S, fill="#3b3f46", outline="#111111",
                    width=S)
d.rectangle((x - 2 * S, y - 12 * S, x + 2 * S, y + 12 * S), fill="#d9dde3")
try:
    font = ImageFont.truetype("arial.ttf", 15 * S)
except OSError:
    font = ImageFont.load_default()
label = "SRS shutter"
tw = d.textlength(label, font=font)
d.text((x - tw / 2, y - h // 2 - 22 * S), label, fill="#000000", font=font)

# photodiode at the end of the multipass output beam (monitors shutter on/off); symbol only, no label
px, py = PD_X * S, (PD_Y - CROP[1]) * S
pw, ph = 18 * S, 28 * S
d.rectangle((px + pw // 3, py - ph // 2, px + pw, py + ph // 2), fill="#3b3f46", outline="#111111", width=S)
d.pieslice((px, py - ph // 3, px + 2 * pw // 3, py + ph // 3), 90, 270, fill="#d9dde3", outline="#111111", width=S)

OUT.parent.mkdir(parents=True, exist_ok=True)
im.save(OUT)
print("wrote", OUT, im.size)
