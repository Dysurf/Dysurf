from pathlib import Path
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
PAT = re.compile(r"(?<=\d)([+-]\d{2,3})(?=\s|$)")

# ---------- real grid parameters ----------
ne = 320
deltaE = 0.05

nqh = 100
deltaH = 0.005


def load_matrix(path: Path) -> np.ndarray:
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            fixed = PAT.sub(r"E\1", line)
            vals = np.fromstring(fixed, sep=" ")
            if vals.size:
                rows.append(vals)
    return np.vstack(rows)


def safe_log10(arr: np.ndarray, floor: float = 1e-12) -> np.ndarray:
    return np.log10(np.maximum(arr, floor))


def summarize(ref: np.ndarray, arr: np.ndarray):
    diff = arr - ref
    rel_l2 = np.linalg.norm(diff) / max(np.linalg.norm(ref), 1e-30)
    return {
        "shape": arr.shape,
        "peak": float(np.max(arr)),
        "sum": float(np.sum(arr)),
        "rel_l2_vs_cncs12": float(rel_l2),
    }


def build_axes(arr: np.ndarray):
    nq, ne_arr = arr.shape

    if ne_arr not in (320, 321):
        print(f"[WARN] unexpected energy dimension: {ne_arr}")
    if nq != 100:
        print(f"[WARN] unexpected Q dimension: {nq}")

    e_axis = np.arange(ne_arr) * deltaE
    q_axis = np.arange(nq, dtype=float) * deltaH
    return q_axis, e_axis


def main():
    cn = load_matrix(ROOT / "SQE_300K_CNCS12.dat")
    psf20 = load_matrix(ROOT / "SQE_300K_psf_hcol_20.dat")
    psf80 = load_matrix(ROOT / "SQE_300K_psf_hcol_80_120.dat")

    q_axis, e_axis = build_axes(cn)
    nq, ne_arr = cn.shape  # ✅ FIXED: define nq and ne_arr from actual data shape
    qmid = nq // 3
    qmid_value = q_axis[qmid]

    with open(ROOT / "compare_cncs12_psf_hcol.txt", "w", encoding="utf-8") as f:
        f.write("CNCS12 vs PSF(hcol-only) comparison\n")
        f.write("=================================\n")
        f.write(f"Q range = 0 to {q_axis[-1]:.6f}\n")
        f.write(f"E range = 0 to {e_axis[-1]:.6f} meV\n")
        f.write(f"qmid index = {qmid}\n")
        f.write(f"qmid value = {qmid_value:.6f}\n\n")

        for name, arr in [
            ("CNCS12", cn),
            ("psf-[20,20,20,20]", psf20),
            ("psf-[80,120,80,120]", psf80),
        ]:
            s = summarize(cn, arr) if name != "CNCS12" else {
                "shape": arr.shape,
                "peak": float(np.max(arr)),
                "sum": float(np.sum(arr)),
                "rel_l2_vs_cncs12": 0.0
            }
            f.write(f"[{name}]\n")
            f.write(f"shape = {s['shape']}\n")
            f.write(f"peak = {s['peak']:.12e}\n")
            f.write(f"sum = {s['sum']:.12e}\n")
            f.write(f"rel_l2_vs_cncs12 = {s['rel_l2_vs_cncs12']:.12e}\n\n")

    datasets = [
        ("CNCS12", cn),
        ("psf-[20,20,20,20]", psf20),
        ("psf-[80,120,80,120]", psf80),
    ]

    vmin = 0.0
    vmax = 7.5

    # ✅ FIXED: use ne_arr and nq from actual shape
    extent = [0.0, nq * deltaH,0.0, ne_arr * deltaE]

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for j, (title, arr) in enumerate(datasets):
        arr_log = safe_log10(arr)
        im = axes[0, j].imshow(
            arr_log.T,
            origin="lower",
            aspect="auto",
            cmap="magma",
            vmin=vmin,
            vmax=vmax,
            extent=extent,
        )
        axes[0, j].set_title(f"{title} (log10)")
        axes[0, j].set_ylabel("E (meV)")
        axes[0, j].set_xlabel("q")
        fig.colorbar(im, ax=axes[0, j], fraction=0.046, pad=0.04, label="log10(Intensity)")

        axes[1, j].plot(e_axis, arr[qmid, :], lw=2.0)
        axes[1, j].set_title(f"{title}: q cut at q = {qmid_value:.3f}")
        axes[1, j].set_xlabel("E (meV)")
        axes[1, j].set_ylabel("Intensity")

    fig.savefig(ROOT / "compare_cncs12_psf_hcol.png", dpi=180)

    fig2, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    ax.plot(e_axis, cn[qmid, :], lw=2.0, label="CNCS12")
    ax.plot(e_axis, psf20[qmid, :], lw=2.0, label="psf-[20,20,20,20]")
    ax.plot(e_axis, psf80[qmid, :], lw=2.0, label="psf-[80,120,80,120]")
    ax.set_xlabel("E (meV)")
    ax.set_ylabel("Intensity")
    ax.set_title(f"q cut overlay at q = {qmid_value:.3f}")
    ax.legend(frameon=False)
    fig2.savefig(ROOT / "compare_cncs12_psf_hcol_linecut.png", dpi=180)

    fig3, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    for ax, (title, arr) in zip(axes, datasets):
        arr_log = safe_log10(arr)
        im = ax.imshow(
            arr_log.T,
            origin="lower",
            aspect="auto",
            cmap="magma",
            vmin=vmin,
            vmax=vmax,
            extent=extent,
        )
        ax.set_title(f"{title} (log10)")
        ax.set_ylabel(r"$E$ (meV)")
        ax.set_xlabel("q")
        fig3.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="log10(Intensity)")
    fig3.savefig(ROOT / "compare_cncs12_psf_hcol_maps_only.png", dpi=180)


if __name__ == "__main__":
    main()