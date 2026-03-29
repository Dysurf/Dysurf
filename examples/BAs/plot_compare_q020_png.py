#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from matplotlib.colors import LogNorm


ROOT = Path(__file__).resolve().parent
Q_MAX = 0.5
E_MAX = 120.0
Q_CUT = 0.20
DELTA_E = 0.05


def load_matrix(path):
    return np.loadtxt(path)


def load_qvals(path):
    omega = np.loadtxt(path)
    return omega[:, 0]


def crop_q_window(sqe, qvals, qmax):
    mask = qvals <= qmax + 1.0e-12
    return sqe[mask], qvals[mask]


def crop_e_window(sqe, emax, delta_e):
    max_idx = min(int(round(emax / delta_e)), sqe.shape[1] - 1)
    energies = np.arange(max_idx + 1, dtype=float) * delta_e
    return sqe[:, : max_idx + 1], energies


def pick_q_cut_index(qvals, q_target):
    return int(np.argmin(np.abs(qvals - q_target)))

def main():
    sqe_born = load_matrix(ROOT / "SQE_300K_born.dat")
    sqe_noborn = load_matrix(ROOT / "SQE_300K_no_born.dat")
    qvals = load_qvals(ROOT / "omega_born.dat")

    sqe_born, qvals = crop_q_window(sqe_born, qvals, Q_MAX)
    sqe_noborn, _ = crop_q_window(sqe_noborn, load_qvals(ROOT / "omega_no_born.dat"), Q_MAX)

    sqe_born, energies = crop_e_window(sqe_born, E_MAX, DELTA_E)
    sqe_noborn, _ = crop_e_window(sqe_noborn, E_MAX, DELTA_E)

    q_idx = pick_q_cut_index(qvals, Q_CUT)
    q_actual = float(qvals[q_idx])

    positive = np.concatenate([sqe_born[sqe_born > 0.0], sqe_noborn[sqe_noborn > 0.0]])
    vmin = max(float(np.min(positive)), 1.0e-12)
    vmax = float(np.max(positive))

    fig = plt.figure(figsize=(12, 9))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 0.85], hspace=0.32, wspace=0.28)

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    extent = [0.0, Q_MAX, 0.0, E_MAX]
    common_imshow = dict(
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap="magma",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        interpolation="nearest",
    )

    im0 = ax0.imshow(sqe_born.T, **common_imshow)
    ax0.set_title("With BORN")
    ax0.set_xlabel("q")
    ax0.set_ylabel("Energy (meV)")
    ax0.set_xlim(0.0, Q_MAX)
    ax0.set_ylim(0.0, E_MAX)
    ax0.axvline(q_actual, color="cyan", lw=1.2, ls="--", alpha=0.9)

    im1 = ax1.imshow(sqe_noborn.T, **common_imshow)
    ax1.set_title("Without BORN")
    ax1.set_xlabel("q")
    ax1.set_ylabel("Energy (meV)")
    ax1.set_xlim(0.0, Q_MAX)
    ax1.set_ylim(0.0, E_MAX)
    ax1.axvline(q_actual, color="cyan", lw=1.2, ls="--", alpha=0.9)

    cbar = fig.colorbar(im1, ax=[ax0, ax1], shrink=0.92, pad=0.03)
    cbar.set_label(r"$S(\mathbf{Q}, E)$")
    
    born_cut = np.clip(sqe_born[q_idx], vmin, None)
    noborn_cut = np.clip(sqe_noborn[q_idx], vmin, None)
    diff_cut = sqe_born[q_idx] - sqe_noborn[q_idx]

    ax2.plot(energies, born_cut, lw=2.0, color="#d62728", label="With BORN")
    ax2.plot(energies, noborn_cut, lw=2.0, color="#1f77b4", label="Without BORN")
    ax2.set_title(f"Energy Cut at q = {Q_CUT:.2f}")
    ax2.set_xlabel("Energy (meV)")
    ax2.set_ylabel(r"$S(\mathbf{Q}, E)$")
    ax2.set_xlim(0.0, E_MAX)
    ax2.set_yscale("log")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(frameon=False)

    ax3.plot(energies, diff_cut, lw=2.0, color="#2ca02c", label="With BORN - Without BORN")
    ax3.axhline(0.0, color="0.4", lw=1.0, ls="--")
    ax3.set_title(f"Difference at q = {Q_CUT:.2f}")
    ax3.set_xlabel("Energy (meV)")
    ax3.set_ylabel(r"$\Delta S(\mathbf{Q}, E)$")
    ax3.set_xlim(0.0, E_MAX)
    linthresh = max(np.max(np.abs(diff_cut)) * 1.0e-3, 1.0e-8)
    ax3.set_yscale("symlog", linthreshy=linthresh)
    ax3.grid(True, which="both", alpha=0.25)
    ax3.legend(frameon=False)

    # fig.suptitle("BAs Comparison Using Current Input", fontsize=14)

    out_png = ROOT / "BAs_compare_q0p20.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    with open(ROOT / "BAs_compare_q0p20_summary.txt", "w", encoding="utf-8") as handle:
        handle.write("BAs compare summary\n")
        handle.write(f"Q range: 0.0 to {Q_MAX:.1f}\n")
        handle.write(f"Energy range: 0.0 to {E_MAX:.1f} meV\n")
        handle.write(f"Requested q cut: {Q_CUT:.4f}\n")
        handle.write(f"Actual q cut: {q_actual:.6f}\n")
        handle.write(f"Number of q points in window: {len(qvals)}\n")
        handle.write("Line-cut y-scale: log\n")
        handle.write(f"Difference y-scale: symlog (linthresh={linthresh:.6e})\n")
        handle.write(f"Max abs difference at q-cut: {np.max(np.abs(diff_cut)):.6e}\n")


if __name__ == "__main__":
    main()
