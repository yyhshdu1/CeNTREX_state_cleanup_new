"""Phase EOM + polarization EOM at the same frequency: why the sidebands of each polarization move.

Jones model. A resonant Pockels cell (polarization EOM) with eigenaxes x', y' and linear input at 45 deg applies a
retardance delta(t) = delta0 sin(w t + theta): x' gets +delta/2, y' gets -delta/2. The phase EOM before it applies
phi(t) = beta sin(w t) to both. So in the x', y' basis the two phases simply add:

    E_x' ~ exp(i[beta sin wt + (delta0/2) sin(wt + theta)]),   E_y' ~ exp(i[beta sin wt - (delta0/2) sin(wt + theta)])

and each is a pure phase modulation with depth |beta +- (delta0/2) e^{i theta}| (same frequency -> phasor sum).
The two switched polarizations are A, B = (E_x' +- i E_y')/sqrt2 (circular after the cell; orthogonal linear after a
quarter-wave plate). Their intensities I_A,B = (1 +- sin delta)/2 are gated at w, in a fixed phase relation to the
instantaneous frequency beta w cos(wt) of the phase EOM. A therefore sees mostly one half of the frequency sweep and B
the other half: mean frequency  <w>_A,B = +- beta w J1(delta0) sin(theta)   (exact; checked numerically to 1e-15).

    python scripts/eom_pol_phase_sidebands.py        # paper values: beta 3.8, w = Gamma, delta0 = pi/2
    -> docs/state_cleanup/figures/EOM_pol_phase_sidebands.{png,pdf}
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import j1

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "docs" / "state_cleanup" / "figures"

import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK, MUTED = "#262624", "#6b6b68"
CA, CB = "#2a78d6", "#eb6834"
mpl.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
                     "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "axes.labelcolor": INK, "text.color": INK})

GAMMA_MHZ = 1.608          # Gamma/2pi
BETA, DELTA0 = 3.8, np.pi / 2
NT = 4096                  # samples per modulation period (one period gives the exact line spectrum)


def fields(theta, beta=BETA, delta0=DELTA0, common=0.0):
    """A, B over one period; common = fraction of delta/2 that the crystal adds to both axes (common-mode PM)."""
    x = 2 * np.pi * np.arange(NT) / NT
    d = delta0 * np.sin(x + theta)
    phi = beta * np.sin(x) + common * d / 2
    ex, ey = np.exp(1j * (phi + d / 2)) / np.sqrt(2), np.exp(1j * (phi - d / 2)) / np.sqrt(2)
    return (ex + 1j * ey) / np.sqrt(2), (ex - 1j * ey) / np.sqrt(2)


def lines(e, nmax=14):
    c = np.fft.fft(e) / NT
    n = np.arange(-nmax, nmax + 1)
    return n, np.abs(c[n]) ** 2


def centroid(e):
    n, p = lines(e, 40)
    return (n * p).sum() / p.sum()


def main():
    thetas = np.linspace(0, 2 * np.pi, 73)
    num = np.array([[centroid(e) for e in fields(t)] for t in thetas])
    ana = BETA * j1(DELTA0) * np.sin(thetas)
    err = np.max(np.abs(num[:, 0] - ana))
    print(f"centroid of A, numeric vs beta J1(delta0) sin(theta): max difference {err:.2e} (units of w)")
    print(f"max shift {BETA * j1(DELTA0) * GAMMA_MHZ:.2f} MHz per polarization (A-B {2 * BETA * j1(DELTA0) * GAMMA_MHZ:.2f} MHz)")
    for t in (0, np.pi / 2):
        print(f"theta {np.degrees(t):3.0f}: eigen-axis depths {abs(BETA + DELTA0 / 2 * np.exp(1j * t)):.2f}, "
              f"{abs(BETA - DELTA0 / 2 * np.exp(1j * t)):.2f}")

    fig = plt.figure(figsize=(11, 6.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.05], hspace=0.5, wspace=0.32)
    cases = [(None, r"phase EOM only ($\beta=3.8$)"), (0.0, r"+ pol. EOM, $\theta=0$"),
             (np.pi / 2, r"+ pol. EOM, $\theta=90^\circ$")]
    for k, (theta, title) in enumerate(cases):
        ax = fig.add_subplot(gs[k // 2, k % 2]) if k < 2 else fig.add_subplot(gs[1, 0])
        if theta is None:
            x = 2 * np.pi * np.arange(NT) / NT
            e = np.exp(1j * BETA * np.sin(x))
            n, p = lines(e)
            ax.vlines(n * GAMMA_MHZ, 0, p, color=MUTED, lw=2)
        else:
            A, B = fields(theta)
            for e, col, lab, off in ((A, CA, "A", -0.17), (B, CB, "B", 0.17)):
                n, p = lines(e)
                ax.vlines(n * GAMMA_MHZ + off, 0, p, color=col, lw=2, label=f"pol. {lab}: mean {centroid(e) * GAMMA_MHZ + 0.0:+.2f} MHz".replace("-0.00", "0.00"))
            ax.legend(fontsize=7.5, frameon=False, loc="upper right")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("frequency − carrier (MHz)")
        ax.set_ylabel("line power")
        ax.set_xlim(-14 * GAMMA_MHZ, 14 * GAMMA_MHZ)
        ax.set_ylim(0, 0.32)
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(np.degrees(thetas), num[:, 0] * GAMMA_MHZ, color=CA, lw=2, label="pol. A")
    ax.plot(np.degrees(thetas), num[:, 1] * GAMMA_MHZ, color=CB, lw=2, label="pol. B")
    ax.plot(np.degrees(thetas), ana * GAMMA_MHZ, color=INK, lw=0.8, ls="--", label=r"$\pm\beta\,\omega\,J_1(\delta_0)\sin\theta$")
    ax.axhline(0, color=MUTED, lw=0.5)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel(r"RF phase $\theta$ of pol. EOM vs phase EOM (deg)")
    ax.set_ylabel("mean frequency (MHz)")
    ax.set_title("mean frequency of each polarization", fontsize=9)
    ax.set_ylim(-5.6, 5.6)
    ax.legend(fontsize=7, frameon=False, loc="upper center", ncol=3, handlelength=1.5, columnspacing=1)
    ax = fig.add_subplot(gs[:, 2]); ax.axis("off")
    txt = (r"$\bf{Model}$ (Jones, one period, exact line spectrum)" "\n"
           r"phase EOM: $\phi=\beta\sin\omega t$, $\beta=3.8$, $\omega=\Gamma$" "\n"
           r"pol. EOM: retardance $\delta=\delta_0\sin(\omega t+\theta)$, $\delta_0=\pi/2$" "\n"
           r"   input at 45$^\circ$ to its axes $x', y'$" "\n\n"
           r"$\bf{1.}$ Along $x'$, $y'$ the phases add:" "\n"
           r"   $\phi_\pm=\beta\sin\omega t\pm\frac{\delta_0}{2}\sin(\omega t+\theta)$" "\n"
           r"   = PM with depth $|\beta\pm\frac{\delta_0}{2}e^{i\theta}|$" "\n"
           r"   ($\theta=0$: 4.59 and 3.01; $\theta=90^\circ$: 3.88 both)" "\n\n"
           r"$\bf{2.}$ The switched pair A, B $=(E_{x'}\pm iE_{y'})/\sqrt{2}$" "\n"
           r"   has $I_{A,B}=(1\pm\sin\delta)/2$, gated at $\omega$ in" "\n"
           r"   step with the PM frequency sweep $\beta\omega\cos\omega t$." "\n"
           r"   A gets one half of the sweep, B the other:" "\n"
           r"   $\langle\nu\rangle_{A,B}=\pm\beta\,\nu_{mod}\,J_1(\delta_0)\sin\theta$" "\n"
           r"   $=\pm3.46$ MHz $\cdot\sin\theta$  (A$-$B up to 6.9 MHz $\approx4\Gamma$)" "\n\n"
           r"$\bf{3.}$ Not locked ($\omega_{pol}=\omega+\Delta$): $\theta=\Delta t$," "\n"
           r"   the shift oscillates at $\Delta$ and the sidebands" "\n"
           r"   of each polarization walk back and forth." "\n\n"
           r"$\bf{4.}$ Total spectrum (no polarizer) does not move:" "\n"
           r"   $I_A+I_B$ is constant. Only a polarization-" "\n"
           r"   selective view (polarizer, or the molecule's" "\n"
           r"   transitions) sees the shift.")
    ax.text(0.0, 1.0, txt, va="top", ha="left", fontsize=8.6, transform=ax.transAxes, linespacing=1.35)
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"EOM_pol_phase_sidebands.{ext}", dpi=200, bbox_inches="tight")
    print("wrote", FIG / "EOM_pol_phase_sidebands.png")


if __name__ == "__main__":
    main()
