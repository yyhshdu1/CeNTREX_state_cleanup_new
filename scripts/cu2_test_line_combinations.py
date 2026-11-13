"""CU2 test 2026-11-08: J=2 left after CU2 per hyperfine group, from combinations of the four detection lines --
every number and figure of docs/state_cleanup/cu2_residual_from_line_combinations_2026-11-09_v2.{html,pdf}.

    python -I scripts/cu2_test_line_combinations.py
    -> docs/state_cleanup/figures/cu2_residual_lines/{fig1_eta,fig2_combinations,fig3_onoff}.png and numbers.json

Method (cu2_residual_from_line_combinations_2026-11-09, sections 2 and 7.4):
  sum_g eta_on[l, g] n_g = R_l sum_g eta_off[l, g] f_g,  f = CU2-off populations (lens: f_t on F1=5/2 F=3 mF = 0, +-1,
  1 - f_t on F1=5/2 F=2 mF = 0; F1=3/2 empty unless stated), f_t from S_B/S_A.
  eta_off: probe photons per molecule (cu2_detection_test.level_photons, plates off, sigma+ about k, v_z 156 m/s),
           lens m_F, vx weights G(vx) (sigma 0.65 m/s, |vx| <= 1 m/s) x y0 acceptance (sigma 10.35 mm, |y0| <= 8.5 mm).
  eta_on:  the same probe yields weighted with the CU2 OBE residual of every start level, per vx, y0 and final level
           (transfer stage, NLTL2 2068 = delta0 13 MHz, f_AOM 209 MHz, 20 + 30 mW, 217 V/cm, 5 % per pass,
           v_z 156 m/s): the residual's own m_F and vx distribution. Molecules CU2 pumps into F1=3/2 keep the lens vx
           (photon recoil ~ 7 mm/s), so no separate sigma_vx is needed for them.
  Errors: Monte Carlo over the four R_l (Gaussian) and S_B/S_A = 1.56 +- 0.17, 20000 draws.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import cu2_test_start_fit as SF  # noqa: E402
import plot_cu2_test_vs_sim as P  # noqa: E402

T = SF.T
G = list(SF.GROUPS)
GN = {(2.5, 3): "F1=5/2 F=3", (2.5, 2): "F1=5/2 F=2", (1.5, 2): "F1=3/2 F=2", (1.5, 1): "F1=3/2 F=1"}
LINES = ["F4", "R7o2F3", "F1p5F3", "F1p5F2"]
LNAME = ["A R(2) F1'=7/2 F'=4", "B R(2) F1'=7/2 F'=3", "C R(2) F1'=5/2 F'=3", "D R(2) F1'=5/2 F'=2"]
R = np.array([0.0019, 0.0057, 0.026, 0.013])
RERR = np.array([0.0059, 0.0037, 0.014, 0.019])
RBA, RBA_ERR = 1.56, 0.17
NMC = 20000
TRANSFER = ROOT / "tmp" / "cu2_power" / "transfer_loss5_d13_f209_20_30_e217_acc_result.npz"
OUT = ROOT / "docs" / "state_cleanup" / "figures" / "cu2_residual_lines"


class Model:
    def __init__(self, power_mW=24.0):
        z = np.load(TRANSFER, allow_pickle=True)
        self.lab = [(float(a), int(b), int(c)) for a, b, c in zip(z["F1"], z["F"], z["mF"])]
        self.st = [(float(a), int(b), int(c)) for a, b, c in zip(z["start_F1"], z["start_F"], z["start_mF"])]
        self.vx, self.y0 = np.asarray(z["vx"], float), np.asarray(z["y0"], float)
        self.W = P.acceptance_weights(self.vx, self.y0, svx=0.65)
        self.wv = self.W.sum(1)
        self.fin = np.asarray(z["final"])
        self.mirror = [self.lab.index((l[0], l[1], -l[2])) for l in self.lab]
        T.POL = "k"
        self.Y = np.array([T.level_photons(k, self.lab, self.vx, power=1e-3 * power_mW) for k in LINES])  # [l, vx, lev]
        self.gi = {g: [i for i, l in enumerate(self.lab) if (l[0], l[1]) == g] for g in G}
        self.power = power_mW

    def start_ms(self, g):
        return SF.GROUPS[g] if g[0] == 2.5 else tuple(range(-g[1], g[1] + 1))       # F1=3/2: all m_F

    def final(self, s):
        if s[2] >= 0:
            return self.fin[self.st.index(s)]
        return self.fin[self.st.index((s[0], s[1], -s[2]))][..., self.mirror]

    def pv(self, g0, mode):
        """[vx, level] residual of group g0 (1 molecule, uniform over its start m_F)."""
        ms = self.start_ms(g0)
        out = np.zeros((len(self.vx), len(self.lab)))
        for m in ms:
            Pf = self.final((g0[0], g0[1], m)) / len(ms)
            out += np.einsum("vy,vyl->vl", self.W, Pf)
        if mode == "lensvx":                                       # same totals, lens vx
            out = np.outer(self.wv, out.sum(0))
        return out

    def parts(self, mode="res", mf="obe"):
        """per start group: (num[l, g], den[g]) of eta_on; eta_on = sum_s w_s num_s / sum_s w_s den_s."""
        res = {}
        for g0 in G:
            pv = self.pv(g0, mode)
            num = np.zeros((len(LINES), len(G))); den = np.zeros(len(G))
            for j, g in enumerate(G):
                idx = self.gi[g]
                if mf == "same":                                   # m_F as before CU2, residual vx of the group
                    tot_v = pv[:, idx].sum(1)
                    ms = self.start_ms(g)
                    w_m = np.array([1.0 / len(ms) if self.lab[i][2] in ms else 0.0 for i in idx])
                    p = np.outer(tot_v, w_m)
                else:
                    p = pv[:, idx]
                den[j] = p.sum()
                num[:, j] = np.einsum("vi,lvi->l", p, self.Y[:, :, idx])
            res[g0] = (num, den)
        return res

    def eta_off(self):
        return np.array([[np.mean([self.wv @ self.Y[l, :, self.lab.index((g[0], g[1], m))] for m in self.start_ms(g)])
                          for g in G] for l in range(len(LINES))])


def eta_on_from(parts, f):
    num = sum(w * parts[g][0] for g, w in zip(G, f) if w > 0)
    den = sum(w * parts[g][1] for g, w in zip(G, f) if w > 0)
    return num / np.where(den > 0, den, 1.0)


def ft_of(eoff, r):
    a0, a1, b0, b1 = eoff[0, 0], eoff[0, 1], eoff[1, 0], eoff[1, 1]
    return (b1 - r * a1) / (r * (a0 - a1) - (b0 - b1))


def solve(eon, eoff, f, Rv, lines=(0, 1, 2, 3)):
    """groups solved = the groups of the lines used (A->F=3, B->F=2, C->3/2 F=2, D->3/2 F=1); the others set to 0."""
    li = list(lines)
    n = np.zeros(4)
    n[li] = np.linalg.solve(eon[np.ix_(li, li)], (Rv * (eoff @ f))[li])
    return n


def run_mc(M, parts, eoff, f_of_r, lines=(0, 1, 2, 3), nmc=NMC, seed=1):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(nmc):
        Rv = R + RERR * rng.standard_normal(4)
        f = f_of_r(RBA + RBA_ERR * rng.standard_normal())
        out.append(solve(eta_on_from(parts, f), eoff, f, Rv, lines))
    return np.array(out)


def lens_f(eoff):
    return lambda r: np.array([ft_of(eoff, r), 1 - ft_of(eoff, r), 0.0, 0.0])


def summarize(n0, S, f):
    e = S.std(0); tot = S.sum(1)
    d = {GN[g]: [float(n0[j]), float(e[j])] for j, g in enumerate(G)}
    d["total"] = [float(n0.sum()), float(tot.std())]
    d["total_upper90"] = float(np.quantile(tot, 0.9))
    for j, g in enumerate(G):
        if f[j] > 0:
            d[f"eff {GN[g]}"] = [float(1 - n0[j] / f[j]), float(e[j] / f[j])]
    d["eff J=2"] = [float(1 - n0.sum()), float(tot.std())]
    return d


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    num = {}
    M = Model(24.0)
    eoff = M.eta_off()
    ft0 = ft_of(eoff, RBA)
    f0 = np.array([ft0, 1 - ft0, 0, 0])
    parts = M.parts("res", "obe")
    eon = eta_on_from(parts, f0)
    eon_lensvx = eta_on_from(M.parts("lensvx", "obe"), f0)
    num["ft"] = float(ft0)
    num["eta_off"] = eoff.tolist(); num["eta_on"] = eon.tolist(); num["eta_on_lensvx"] = eon_lensvx.tolist()
    # model residual (OBE) for the lens start
    resid = sum(w * M.pv(g, "res").sum(0) for g, w in zip(G, f0) if w > 0)
    n_mod = np.array([resid[M.gi[g]].sum() for g in G])
    num["model"] = {GN[g]: float(v) for g, v in zip(G, n_mod)}; num["model"]["total"] = float(n_mod.sum())
    # transfer matrix: from each lens group into each group (per molecule of that start)
    num["transfer"] = {GN[g0]: {GN[g]: float(M.pv(g0, "res").sum(0)[M.gi[g]].sum()) for g in G} for g0 in G[:2]}
    # main result
    n0 = solve(eon, eoff, f0, R)
    S = run_mc(M, parts, eoff, lens_f(eoff))
    num["main"] = summarize(n0, S, f0)
    num["corr"] = np.corrcoef(S.T).tolist()
    # without the residual vx (document's former main result, for comparison)
    pl = M.parts("lensvx", "obe")
    S2 = run_mc(M, pl, eoff, lens_f(eoff))
    num["lensvx"] = summarize(solve(eta_on_from(pl, f0), eoff, f0, R), S2, f0)
    # line combinations (figure 2)
    combos = {"A": (0,), "A+B": (0, 1), "A+C": (0, 2), "A+B+C": (0, 1, 2), "A+B+C+D": (0, 1, 2, 3)}
    num["combos"] = {}
    for name, ls in combos.items():
        Sc = run_mc(M, parts, eoff, lens_f(eoff), ls, nmc=5000)
        nc = solve(eon, eoff, f0, R, ls)
        num["combos"][name] = {"n": nc.tolist(), "err": Sc.std(0).tolist(), "lines": list(ls),
                               "total": float(nc.sum()), "total_err": float(Sc.sum(1).std())}
    # model on/off ratio of each line and the signal shares (figure 3)
    on = np.einsum("lvi,vi->l", M.Y, sum(w * M.pv(g, "res") for g, w in zip(G, f0) if w > 0))
    off = eoff @ f0
    num["model_R"] = (on / off).tolist()
    sh_off = eoff * f0 / off[:, None]
    sh_on = eon * n_mod / (eon @ n_mod)[:, None]
    num["share_off"] = sh_off.tolist(); num["share_on"] = sh_on.tolist()
    # section 6: f_t and F1=3/2 before CU2
    num["ft_scan"] = []
    for ft in (0.30, 0.35, 0.41, ft0, 0.50, 0.60):
        f = np.array([ft, 1 - ft, 0, 0])
        Sf = run_mc(M, parts, eoff, lambda r, f=f: f, nmc=5000)
        d = summarize(solve(eta_on_from(parts, f), eoff, f, R), Sf, f); d["ft"] = float(ft)
        num["ft_scan"].append(d)
    num["f32_before"] = []
    for g, x in (((1.5, 2), 0.005), ((1.5, 2), 0.01), ((1.5, 2), 0.02), ((1.5, 2), 0.05), ((1.5, 1), 0.01),
                 ((1.5, 1), 0.05)):
        f = np.array([0.41, 0.59, 0, 0]) * (1 - x)
        f[G.index(g)] = x
        Sf = run_mc(M, parts, eoff, lambda r, f=f: f, nmc=5000)
        d = summarize(solve(eta_on_from(parts, f), eoff, f, R), Sf, f); d["case"] = f"{GN[g]}: {100 * x:g} %"
        num["f32_before"].append(d)
    # section 5: probe power x vx of the residual x m_F after CU2
    num["syst"] = []
    for pw in (15.0, 24.0):
        Mp = M if pw == 24.0 else Model(pw)
        eo = Mp.eta_off()
        for mode, mf in (("res", "obe"), ("res", "same"), ("lensvx", "obe"), ("lensvx", "same")):
            pp = Mp.parts(mode, mf)
            ftp = ft_of(eo, RBA)
            f = np.array([ftp, 1 - ftp, 0, 0])
            Sp = run_mc(Mp, pp, eo, lens_f(eo), nmc=5000)
            d = summarize(solve(eta_on_from(pp, f), eo, f, R), Sp, f)
            d.update(power=pw, vx="CU2 OBE" if mode == "res" else "lens", mF="CU2 OBE" if mf == "obe" else "as before",
                     ft=float(ftp))
            num["syst"].append(d)
    (OUT / "numbers.json").write_text(json.dumps(num, indent=1))
    figures(num)
    print(json.dumps({k: num[k] for k in ("ft", "model", "main", "lensvx", "model_R")}, indent=1))


def figures(num):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 9,
                         "savefig.facecolor": "white"})
    cols = ["F1=5/2 F=3\n(triplet start)", "F1=5/2 F=2\n(singlet start)", "F1=3/2 F=2", "F1=3/2 F=1"]
    # fig 1: eta matrices
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.9), constrained_layout=True)
    for ax, key, title in zip(axs, ("eta_off", "eta_on"), ("CU2 off: lens m_F, lens vx",
                                                             "CU2 on: m_F and vx of the CU2-OBE residual")):
        a = np.array(num[key])
        im = ax.imshow(a, cmap="Blues", vmin=0, vmax=2.6)
        for i in range(4):
            for j in range(4):
                ax.text(j, i, f"{a[i, j]:.2f}", ha="center", va="center", color="white" if a[i, j] > 1.3 else "#111")
        ax.set_xticks(range(4), cols, fontsize=7.5); ax.set_yticks(range(4), LNAME if key == "eta_off" else "ABCD",
                                                                     fontsize=8)
        ax.set_title(title, fontsize=9)
    fig.colorbar(im, ax=axs, shrink=0.8, label="photons per molecule")
    fig.suptitle("Probe yields η[line, group]: 24 mW, σ+ about k, E = 0", fontsize=10)
    fig.savefig(OUT / "fig1_eta.png", dpi=180)
    plt.close(fig)
    # fig 2: combinations
    fig, axs = plt.subplots(1, 5, figsize=(12.5, 3.6), sharey=True, constrained_layout=True)
    order = ["A", "A+B", "A+C", "A+B+C", "A+B+C+D"]
    colors = {"A": "#9ecae1", "A+B": "#6baed6", "A+C": "#f4a261", "A+B+C": "#3182bd", "A+B+C+D": "#08306b"}
    titles = ["F1=5/2 F=3 (triplet start)", "F1=5/2 F=2 (singlet start)", "F1=3/2 F=2", "F1=3/2 F=1",
              "sum of the solved groups"]
    model = [num["model"][GN[g]] for g in G] + [num["model"]["total"]]
    for j, ax in enumerate(axs):
        xs = []
        for c in order:
            d = num["combos"][c]
            if j < 4 and j not in d["lines"]:
                continue
            v, e = (d["total"], d["total_err"]) if j == 4 else (d["n"][j], d["err"][j])
            k = len(xs)
            ax.errorbar(k, 100 * v, 100 * e, fmt="s" if j == 4 else "o", color=colors[c], ms=5, capsize=0)
            xs.append(c)
        ax.axhline(100 * model[j], color="#d62728", ls="--", lw=1, label="CU2 OBE, vz 156 m/s" if j == 0 else None)
        ax.axhline(0, color="#999", lw=0.6)
        ax.set_xticks(range(len(xs)), xs, rotation=90, fontsize=7.5)
        ax.set_xlim(-0.7, len(xs) - 0.3)
        ax.set_title(titles[j], fontsize=8.5)
    axs[0].set_ylabel("left after CU2 (% of the J=2 start)")
    axs[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    axs[0].set_ylim(-1.6, 2.3)
    fig.suptitle("J=2 population left after CU2 (NLTL2 2068), solved from line combinations (labels: lines used)",
                 fontsize=9.5)
    fig.savefig(OUT / "fig2_combinations.png", dpi=180)
    plt.close(fig)
    # fig 3: on/off and signal composition
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.0), gridspec_kw={"width_ratios": [1, 1.25]},
                                 constrained_layout=True)
    x = np.arange(4)
    a1.errorbar(x - 0.08, 100 * R, 100 * RERR, fmt="o", color="#1f77b4", capsize=0, label="measured (template fit)")
    a1.plot(x + 0.08, 100 * np.array(num["model_R"]), "D", color="#d62728", label="CU2 OBE residual → probe OBE")
    a1.axhline(0, color="#999", lw=0.6)
    a1.set_xticks(x, [n.replace(" R(2)", "\nR(2)") for n in LNAME], fontsize=7.5)
    a1.set_ylabel("S_on / S_off (%)"); a1.set_title("on/off ratio of each line", fontsize=9)
    a1.legend(frameon=False, fontsize=8)
    gcol = ["#2e8b57", "#c0392b", "#7d4fa8", "#8c564b"]
    for i in range(4):
        for row, sh, light in ((i * 2, num["share_off"][i], True), (i * 2 + 1, num["share_on"][i], False)):
            left = 0
            for j in range(4):
                a2.barh(-row, 100 * sh[j], left=left, color=gcol[j], alpha=0.4 if light else 1.0, height=0.85)
                left += 100 * sh[j]
    a2.set_yticks([-(2 * i + 0.5) for i in range(4)], "ABCD")
    a2.set_xlabel("share of the line's signal (%)"); a2.set_xlim(0, 100)
    a2.set_title("signal by ground group: pale = CU2 off, dark = CU2 on (model residual)", fontsize=9)
    from matplotlib.patches import Patch
    a2.legend([Patch(color=c) for c in gcol], [GN[g] for g in G], frameon=False, fontsize=7.5, ncol=4,
              loc="upper center", bbox_to_anchor=(0.5, -0.13))
    fig.savefig(OUT / "fig3_onoff.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
