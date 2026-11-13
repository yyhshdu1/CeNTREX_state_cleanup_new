"""Saturation of the DET probe: photons per molecule of each ground group vs probe power.

    python -I scripts/det_power_saturation.py                       # the six CU2-test / planned lines
    python -I scripts/det_power_saturation.py only=P3o2F1,Q1o2F0    # just these lines, own figure
    python -I scripts/det_power_saturation.py only=P3o2F1,Q1o2F0 E=150   # DET plates on: 150 V/cm along k
    -> tmp/detection/power_scan/*.json (probe OBE per line, all powers), and
       docs/state_cleanup/figures/slides/DET_power_saturation[_<lines>].{png,pdf,json}

Same probe OBE as cu2_detection_test.yields (DET plates off: E = 0, earth field 0.5 G; sigma+ about k, run with the
axis along k and each vertical-axis level rotated onto k, POL "k"; flat top 1.5 x 1.5 cm, v_z 156 m/s, 96 us), only
over more powers. J=2 lines: group yield = eta_off of cu2_test_line_combinations, i.e. the mean over the lens
sublevels (F1=5/2 F=3: m_F = 0, +-1; F1=5/2 F=2: m_F = 0; F1=3/2: all m_F); at 24 mW this reproduces eta_off of
numbers.json. J=1 lines (Q(1)): mean over all m_F of each J=1 group (no lens start is defined for J=1).
Both weighted with the lens v_x distribution (Gaussian 0.65 m/s, |v_x| <= 1 m/s, CU2 grid +-0.2, +-0.6, +-1.0 m/s).
Saturation measure: local slope s = d ln(eta) / d ln(P): s = 1 linear in power, s = 0 fully saturated.
E=<V/cm>: DET plates on, field along k. sigma+ about k is then sigma+ about the field axis and the field turns
adiabatically from CU2 (vertical) to DET without passing zero, so every level keeps its m about the field: the
yields are read without the |d(pi/2)|^2 rotation. Opposite-parity B levels are kept (Stark mixing). The E = 0
curves are drawn faintly for comparison.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import cu2_detection_test as D  # noqa: E402
import cu2_test_line_combinations as L  # noqa: E402

POWERS = [0.3e-3, 0.5e-3, 1e-3, 1.5e-3, 2e-3, 3e-3, 5e-3, 7.5e-3, 10e-3, 15e-3, 20e-3, 24e-3, 30e-3, 40e-3, 60e-3,
          100e-3]
TITLE = {"F4": "A  R(2) F1'=7/2 F'=4", "R7o2F3": "B  R(2) F1'=7/2 F'=3", "F1p5F3": "C  R(2) F1'=5/2 F'=3",
         "F1p5F2": "D  R(2) F1'=5/2 F'=2", "P1o2F1": "P(2) F1'=1/2 F'=1 (planned)", "P1o2F0": "P(2) F1'=1/2 F'=0",
         "P3o2F1": "P(2) F1'=3/2 F'=1", "Q1o2F0": "Q(1) F1'=1/2 F'=0 (J=1 line)"}
DEFAULT = ["F4", "R7o2F3", "F1p5F3", "F1p5F2", "P1o2F1", "P1o2F0"]
MAIN = {"F4": (2, 2.5, 3), "R7o2F3": (2, 2.5, 2), "F1p5F3": (2, 1.5, 2), "F1p5F2": (2, 1.5, 1),
        "P1o2F1": (2, 1.5, 2), "P1o2F0": (2, 1.5, 1), "P3o2F1": (2, 2.5, 2), "Q1o2F0": (1, 0.5, 1)}
GROUPS = {2: [(2.5, 3), (2.5, 2), (1.5, 2), (1.5, 1)], 1: [(0.5, 0), (0.5, 1), (1.5, 1), (1.5, 2)]}
GCOL = {(2, 2.5, 3): "#2e8b57", (2, 2.5, 2): "#c0392b", (2, 1.5, 2): "#7d4fa8", (2, 1.5, 1): "#8c564b",
        (1, 0.5, 0): "#1f77b4", (1, 0.5, 1): "#ff7f0e", (1, 1.5, 1): "#17becf", (1, 1.5, 2): "#bcbd22"}
SCAN_DIR = ROOT / "tmp" / "detection" / "power_scan"
E_FIELD = 0.0                                  # V/cm along k
OUT = ROOT / "docs" / "state_cleanup" / "figures" / "slides"


def frac(x):
    return f"{int(2 * x)}/2" if x != int(x) else f"{int(x)}"


def gname(key):
    J, F1, F = key
    return (f"J=1 F1={frac(F1)} F={F}" if J == 1 else f"F1={frac(F1)} F={F}")


def groups_of(line):
    """[(key, levels)] of the line's ground J; J=2: lens m_F for F1=5/2, all m_F for F1=3/2; J=1: all m_F."""
    J = D.ground_J(line)
    out = []
    for F1, F in GROUPS[J]:
        ms = L.SF.GROUPS[(F1, F)] if (J == 2 and F1 == 2.5) else tuple(range(-F, F + 1))
        out.append(((J, F1, F), [(F1, F, m) for m in ms]))
    return out


def scan_cache(line, pol=None, E=None):
    E = E_FIELD if E is None else E
    return SCAN_DIR / (f"cu2_r2_{line}_yields_polspk_pscan" + (f"_E{E:g}" if E else "") + ".json")


def run_missing(lines, E=None):
    """Probe OBE for the powers not yet in the line's cache; new powers are merged into it."""
    E = E_FIELD if E is None else E
    tmp = SCAN_DIR / "_partial.json"
    D.CACHE_DIR, D.E_DET, D.POL, D.THREADS = SCAN_DIR, E, "spk", 16   # the run: sigma+ about its z axis = k
    for line in lines:
        f = scan_cache(line, E=E)
        old = json.loads(f.read_text()) if f.exists() else None
        have = set(old["powers_W"]) if old else set()
        todo = sorted(p for p in POWERS if p not in have)
        if not todo:
            continue
        D.cache, D.POWERS = (lambda line, pol=None: tmp), todo
        try:
            D.SOLVER = "dopri5"
            D.yields(line)
        except ValueError as err:                                      # dopri5 stiffness check at high power
            if "Stiff" not in str(err):
                raise
            print(f"  {line}: dopri5 reports stiffness ({err}); trying tsit5", flush=True)
            try:
                D.SOLVER = "tsit5"
                D.yields(line)
            except ValueError as err2:
                if "Stiff" not in str(err2):
                    raise
                print(f"  {line}: tsit5 too ({err2}); fixed-step dopri5, dt = 1 ns", flush=True)
                D.SOLVER, D.DT = "fixed_dopri5", 1e-9
                D.yields(line)
            D.SOLVER, D.DT = "dopri5", 2e-9
        new = json.loads(tmp.read_text())
        if old:                                                        # merge, sorted by power
            assert old["labels"] == new["labels"] and old["vx_m_s"] == new["vx_m_s"]
            ph = dict(zip(old["powers_W"], old["photons"])) | dict(zip(new["powers_W"], new["photons"]))
            new["powers_W"] = sorted(ph)
            new["photons"] = [ph[p] for p in new["powers_W"]]
            new["wall_s"] = old.get("wall_s", 0) + new["wall_s"]
        f.write_text(json.dumps(new, indent=1))
        tmp.unlink()
        print("merged", f.name, "powers (mW):", [round(p * 1e3, 1) for p in new["powers_W"]], flush=True)


def group_yields(lines, E=None):
    E = E_FIELD if E is None else E
    z = np.load(L.TRANSFER, allow_pickle=True)                        # the CU2 v_x grid and its lens weights
    vx, y0 = np.asarray(z["vx"], float), np.asarray(z["y0"], float)
    wv = L.P.acceptance_weights(vx, y0, svx=0.65).sum(1)
    D.cache, D.CACHE_DIR = (lambda line, pol=None: scan_cache(line, E=E)), SCAN_DIR
    D.POL = "k" if E == 0 else "spk"          # E = 0: rotate vertical levels onto k; plates on: m kept, no rotation
    out = {}
    for line in lines:
        gs = groups_of(line)
        eta, eta0 = np.zeros((len(POWERS), len(gs))), np.zeros((len(POWERS), len(gs)))
        for k, P in enumerate(POWERS):
            for j, (_, lev) in enumerate(gs):
                eta[k, j] = wv @ D.level_photons(line, lev, vx, power=P).mean(1)
                eta0[k, j] = D.level_photons(line, lev, [0.0], power=P).mean()
        out[line] = dict(groups=[gname(g) for g, _ in gs], eta=eta.tolist(), eta_vx0=eta0.tolist())
    return out


def slope(P, y):
    lp, ly = np.log(P), np.log(np.maximum(y, 1e-12))
    return np.gradient(ly, lp)


def plot(res, lines, name, ref=None):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": "white"})
    P = np.array(POWERS) * 1e3
    k24 = POWERS.index(24e-3)
    fig, axs = plt.subplots(2, len(lines), figsize=(3.1 * len(lines) + 0.4, 6.2), sharex=True, squeeze=False,
                            constrained_layout=True, gridspec_kw={"height_ratios": [1.5, 1]})
    summary = {}
    for c, line in enumerate(lines):
        a, b = axs[0, c], axs[1, c]
        eta = np.array(res[line]["eta"])
        summary[line] = {}
        for j, (g, _) in enumerate(groups_of(line)):
            if eta[:, j].max() < 0.02:
                continue
            main = g == MAIN[line]
            lw, ls = (2.2, "-") if main else (1.2, "--")
            if ref is not None:
                a.plot(P, np.array(ref[line]["eta"])[:, j], ls, color=GCOL[g], lw=lw * 0.8, alpha=0.3)
            a.plot(P, eta[:, j], ls, color=GCOL[g], lw=lw, marker="o", ms=3 if main else 2,
                   label=gname(g) + (" (main)" if main else ""))
            s = slope(P, eta[:, j])
            b.plot(P, s, ls, color=GCOL[g], lw=lw)
            p_half = float(np.interp(0.5, s[::-1], P[::-1])) if s.min() < 0.5 < s.max() else None
            summary[line][gname(g)] = dict(eta_24mW=float(eta[k24, j]), slope_24mW=float(s[k24]),
                                           frac_of_Pmax_at_24mW=float(eta[k24, j] / eta[-1, j]),
                                           eta_max=float(eta[:, j].max()),
                                           P_at_max_mW=float(P[int(np.argmax(eta[:, j]))]), P_slope_half_mW=p_half)
            if main and p_half:
                b.axvline(p_half, color=GCOL[g], lw=0.8, ls=":")
                b.annotate(f"s = 0.5 at\n{p_half:.1f} mW", (p_half, 0.5), xytext=(4, 6), textcoords="offset points",
                           fontsize=7.5, color=GCOL[g])
        for ax in (a, b):
            ax.axvline(24, color="#999", lw=0.9)
            ax.set_xscale("log")
        a.set_title(TITLE[line], fontsize=9.5)
        a.set_ylim(bottom=0)
        a.legend(frameon=False, fontsize=7, loc="upper left")
        b.set_ylim(-0.05, 1.05)
        b.axhline(0.5, color="#ccc", lw=0.6)
        b.set_xlabel("DET power (mW)")
        tk = [0.3, 1, 3, 10, 24, 100] + ([500] if P[-1] > 400 else []) + ([1000] if P[-1] > 900 else [])
        b.set_xticks(tk, [f"{x:g}" for x in tk])
    axs[0, 0].set_ylabel("photons per molecule η\n(J=2: lens m$_F$; J=1: all m$_F$; lens v$_x$)")
    axs[1, 0].set_ylabel("local slope d ln η / d ln P\n1 = linear, 0 = saturated")
    sep = "\n" if len(lines) < 4 else " "
    field = (f"plates on, E = {E_FIELD:g} V/cm along k (faint: E = 0)" if E_FIELD else "plates off (E = 0)")
    fig.suptitle(f"DET probe saturation, {field}, σ+ about k," + sep + "flat top 1.5 × 1.5 cm, v$_z$ 156 m/s "
                 "(96 µs); grey line: 24 mW", fontsize=10)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=200)
    (OUT / f"{name}.json").write_text(json.dumps(dict(powers_W=POWERS, groups=res, summary=summary), indent=1))
    print("wrote", OUT / f"{name}.png")
    return summary


if __name__ == "__main__":
    only = next((a[5:].split(",") for a in sys.argv[1:] if a.startswith("only=")), None)
    E_FIELD = next((float(a[2:]) for a in sys.argv[1:] if a.startswith("E=")), 0.0)
    pmax = next((float(a[5:]) for a in sys.argv[1:] if a.startswith("pmax=")), 100.0)   # mW
    POWERS = POWERS + [p for p in (150e-3, 200e-3, 300e-3, 400e-3, 500e-3, 700e-3, 1000e-3) if p <= pmax * 1e-3 + 1e-12]
    lines = only or DEFAULT
    name = ("DET_power_saturation" + ("_" + "_".join(lines) if only else "") + (f"_E{E_FIELD:g}" if E_FIELD else "")
            + (f"_to{pmax:g}mW" if pmax > 100 else ""))
    run_missing(lines)
    res = group_yields(lines)
    ref = None
    if E_FIELD:
        run_missing(lines, E=0.0)                                      # E = 0 reference, same powers
    if E_FIELD and all(scan_cache(x, E=0.0).exists() for x in lines):
        ref = group_yields(lines, E=0.0)
    if not only and not E_FIELD:                                       # check against the solve's eta_off
        eoff = json.loads((L.OUT / "numbers.json").read_text())["eta_off"]
        k24 = POWERS.index(24e-3)
        for i, line in enumerate(lines[:4]):
            print(line, "eta(24 mW)", np.round(res[line]["eta"][k24], 3), " numbers.json eta_off",
                  np.round(eoff[i], 3))
    summ = plot(res, lines, name, ref)
    for line, d in summ.items():
        for g, v in d.items():
            print(f"{line:7s} {g:16s} eta24 {v['eta_24mW']:.3f}  slope24 {v['slope_24mW']:.2f}  "
                  f"eta24/eta(Pmax) {v['frac_of_Pmax_at_24mW']:.2f}  max {v['eta_max']:.3f} at {v['P_at_max_mW']:g} mW"
                  f"  s=0.5 at {v['P_slope_half_mW']}")
