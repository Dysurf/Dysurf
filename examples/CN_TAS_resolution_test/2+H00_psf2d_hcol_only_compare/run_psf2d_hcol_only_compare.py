from pathlib import Path
import re
import subprocess
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
ACC = ROOT.parent / "2+H00_psf2d_minimal_acceptance"
EXE = ROOT.parents[2] / "src" / "src_gfortran" / "validate_psf2d_parametrization"
PAT = re.compile(r"(?<=\d)([+-]\d{2,3})(?=\s|$)")

PRESETS = [
    ("hcol_narrow", "[20,20,20,20]"),
    ("hcol_base",   "[50,80,50,120]"),
    ("hcol_wide",   "[80,120,80,120]"),
]


def load_matrix(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            fixed = PAT.sub(r"E\1", line)
            vals = np.fromstring(fixed, sep=" ")
            if vals.size:
                rows.append(vals)
    if not rows:
        raise ValueError(f"No numeric rows loaded from {path}")
    return np.vstack(rows)


def parse_vector_from_file(path: Path, key: str) -> list:
    text = path.read_text(encoding="utf-8")
    match = re.search(rf"^{re.escape(key)}\s*=\s*\[(.*?)\]\s*$", text, re.MULTILINE)
    if not match:
        raise ValueError(f"Missing {key} in {path}")
    vals = [float(x.strip()) for x in match.group(1).split(",")]
    if len(vals) != 3:
        raise ValueError(f"{key} in {path} must have 3 components")
    return vals


def load_triplets(path: Path):
    arr = load_matrix(path)
    if arr.shape[1] < 3:
        raise ValueError(f"Expected at least 3 columns in {path}, got shape {arr.shape}")

    dq = np.unique(arr[:, 0])
    de = np.unique(arr[:, 1])

    expected = dq.size * de.size
    if arr.shape[0] != expected:
        raise ValueError(
            f"Triplet table in {path} cannot be reshaped cleanly: "
            f"rows={arr.shape[0]}, unique dq={dq.size}, unique dE={de.size}, expected={expected}"
        )

    z = arr[:, 2].reshape(de.size, dq.size)
    return dq, de, z


def choose_representative_point():
    shear_file = ACC / "PSF_shear_base_300K.dat"
    if not shear_file.exists():
        return 6, 142
    shear = load_matrix(shear_file)

    margin_q = 5
    margin_e = 10
    if shear.shape[0] <= 2 * margin_q or shear.shape[1] <= 2 * margin_e:
        raise ValueError(
            f"Shear field too small for requested margins: shape={shear.shape}, "
            f"margin_q={margin_q}, margin_e={margin_e}"
        )

    core = np.abs(shear[margin_q:-margin_q, margin_e:-margin_e])
    idx = np.unravel_index(np.argmax(core), core.shape)

    # Fortran-like 1-based index for downstream executable
    iq = idx[0] + margin_q + 1
    ie = idx[1] + margin_e + 1
    return iq, ie


def write_selection_summary(iq: int, ie: int) -> None:
    h0 = 2.0 + 0.005 * (iq - 1)
    e0 = 0.05 * (ie - 1)
    orient_file = ROOT / "instrument_hcol_base.txt"
    u = parse_vector_from_file(orient_file, "u")
    v = parse_vector_from_file(orient_file, "v")
    text = f"""hcol-only local PSF comparison selection
=================================
selection_mode = automatic
selection_rule = max_abs_shear_inside_margin
iq = {iq}
iE = {ie}
q0_hkl = [{h0:.6f}, 0.000000, 0.000000]
E0_meV = {e0:.6f}
u = {u}
v = {v}

Parameter files:
  instrument_hcol_narrow.txt
  instrument_hcol_base.txt
  instrument_hcol_wide.txt
  example_nominal_point.txt
"""
    (ROOT / "psf2d_hcol_only_selection.txt").write_text(text, encoding="utf-8")


def run_exports(iq: int, ie: int) -> None:
    if not EXE.exists():
        raise SystemExit(f"validate_psf2d_parametrization executable not found: {EXE}")

    for preset, _ in PRESETS:
        prefix = f"psf2d_{preset}"
        orient_file = ROOT / f"instrument_{preset}.txt"
        u = parse_vector_from_file(orient_file, "u")
        v = parse_vector_from_file(orient_file, "v")
        cmd = [str(EXE), str(iq), str(ie), preset, prefix] + [str(x) for x in (u + v)]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, cwd=ROOT, check=True)


def parse_params_file(path: Path) -> dict:
    data = {}
    if not path.exists():
        return data

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()

            try:
                data[key] = float(value)
            except ValueError:
                data[key] = value
    return data


def write_summary(iq: int, ie: int) -> None:
    out = ROOT / "psf2d_hcol_only_summary.txt"
    with open(out, "w", encoding="utf-8") as f:
        f.write("Local PSF comparison with hcol-only variation\n")
        f.write("===========================================\n")
        f.write(f"iq = {iq}\n")
        f.write(f"iE = {ie}\n\n")
        f.write("Only hcol is changed across the three presets.\n")
        f.write("All mosaic values are fixed at mono=35', ana=35', sample=60'.\n\n")

        for preset, label in PRESETS:
            params_path = ROOT / f"psf2d_{preset}_params.txt"
            params = parse_params_file(params_path)

            f.write(f"[{preset}] hcol = {label}\n")
            if not params:
                f.write("  params file missing or empty\n\n")
                continue

            keys_of_interest = [
                "q0",
                "E0",
                "sigma_q",
                "sigma_e_left",
                "sigma_e_right",
                "shear",
                "relative_l2",
                "rmse",
            ]
            for key in keys_of_interest:
                if key in params:
                    f.write(f"  {key} = {params[key]}\n")

            f.write("\nFull params dump:\n")
            f.write(params_path.read_text(encoding="utf-8"))
            f.write("\n")


def plot_compare() -> None:
    fig, axes = plt.subplots(5, 3, figsize=(12, 15), constrained_layout=True)

    for j, (preset, label) in enumerate(PRESETS):
        dq_raw, de_raw, raw = load_triplets(ROOT / f"psf2d_{preset}_raw.dat")
        dq_fit, de_fit, fit = load_triplets(ROOT / f"psf2d_{preset}_fit.dat")
        dq_cut = load_matrix(ROOT / f"psf2d_{preset}_dq_cut.dat")
        de_cut = load_matrix(ROOT / f"psf2d_{preset}_dE_cut.dat")
        params = parse_params_file(ROOT / f"psf2d_{preset}_params.txt")

        if raw.shape != fit.shape:
            raise ValueError(f"raw/fit shape mismatch for {preset}: {raw.shape} vs {fit.shape}")

        title_prefix = f"hcol={label}"
        vmin = min(np.min(raw), np.min(fit))
        vmax = max(np.max(raw), np.max(fit))
        residual = fit - raw
        rabs = np.max(np.abs(residual))

        im0 = axes[0, j].imshow(
            raw,
            origin="lower",
            aspect="auto",
            cmap="magma",
            extent=[dq_raw[0], dq_raw[-1], de_raw[0], de_raw[-1]],
            vmin=vmin,
            vmax=vmax,
        )
        axes[0, j].set_title(f"{title_prefix}\nraw")
        axes[0, j].set_xlabel("dQ")
        axes[0, j].set_ylabel("dE")
        fig.colorbar(im0, ax=axes[0, j], fraction=0.046, pad=0.04)

        im1 = axes[1, j].imshow(
            fit,
            origin="lower",
            aspect="auto",
            cmap="magma",
            extent=[dq_fit[0], dq_fit[-1], de_fit[0], de_fit[-1]],
            vmin=vmin,
            vmax=vmax,
        )
        axes[1, j].set_title(f"{title_prefix}\nfit")
        axes[1, j].set_xlabel("dQ")
        axes[1, j].set_ylabel("dE")
        fig.colorbar(im1, ax=axes[1, j], fraction=0.046, pad=0.04)

        im2 = axes[2, j].imshow(
            residual,
            origin="lower",
            aspect="auto",
            cmap="coolwarm",
            extent=[dq_fit[0], dq_fit[-1], de_fit[0], de_fit[-1]],
            vmin=-rabs,
            vmax=rabs,
        )
        axes[2, j].set_title(f"{title_prefix}\nresidual (fit - raw)")
        axes[2, j].set_xlabel(r"$dQ$")
        axes[2, j].set_ylabel(r"$dE$")
        fig.colorbar(im2, ax=axes[2, j], fraction=0.046, pad=0.04)

        axes[3, j].plot(dq_cut[:, 0], dq_cut[:, 1], lw=2.0, label="raw")
        axes[3, j].plot(dq_cut[:, 0], dq_cut[:, 2], lw=1.6, ls="--", label="fit")
        axes[3, j].set_title(f"{title_prefix}\ndQ cut")
        axes[3, j].set_xlabel(r"$dQ$")
        axes[3, j].set_ylabel("Intensity")
        axes[3, j].legend(frameon=False, fontsize=8)

        axes[4, j].plot(de_cut[:, 0], de_cut[:, 1], lw=2.0, label="raw")
        axes[4, j].plot(de_cut[:, 0], de_cut[:, 2], lw=1.6, ls="--", label="fit")
        axes[4, j].set_title(f"{title_prefix}\ndE cut")
        axes[4, j].set_xlabel(r"$dE$")
        axes[4, j].set_ylabel("Intensity")
        axes[4, j].legend(frameon=False, fontsize=8)

        if params:
            info_lines = []
            for key in ["sigma_q", "sigma_e_left", "sigma_e_right", "shear", "relative_l2", "rmse"]:
                if key in params:
                    val = params[key]
                    if isinstance(val, float):
                        info_lines.append(f"{key}={val:.4g}")
                    else:
                        info_lines.append(f"{key}={val}")
            if info_lines:
                axes[1, j].text(
                    0.02,
                    0.98,
                    "\n".join(info_lines),
                    transform=axes[1, j].transAxes,
                    va="top",
                    ha="left",
                    fontsize=8,
                    color="white",
                    bbox=dict(facecolor="black", alpha=0.35, edgecolor="none"),
                )

    fig.savefig(ROOT / "psf2d_hcol_only_compare.png", dpi=180)


def plot_overlay() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

    colors = {
        "hcol_narrow": "#d62728",
        "hcol_base":   "#1f77b4",
        "hcol_wide":   "#2ca02c",
    }

    for preset, label in PRESETS:
        dq_cut = load_matrix(ROOT / f"psf2d_{preset}_dq_cut.dat")
        de_cut = load_matrix(ROOT / f"psf2d_{preset}_dE_cut.dat")
        color = colors[preset]

        axes[0].plot(dq_cut[:, 0], dq_cut[:, 1], marker="o", ls=" ", color=color, label=f"{label} raw")
        axes[0].plot(dq_cut[:, 0], dq_cut[:, 2], lw=1.2, ls="--", color=color, label=f"{label} fit")

        axes[1].plot(de_cut[:, 0], de_cut[:, 1], marker="o", ls=" ", color=color, label=f"{label} raw")
        axes[1].plot(de_cut[:, 0], de_cut[:, 2], lw=1.2, ls="--", color=color, label=f"{label} fit")

    axes[0].set_title("Overlay dQ cuts")
    axes[0].set_xlabel(r"$dQ$")
    axes[0].set_ylabel("Intensity")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].set_title("Overlay dE cuts")
    axes[1].set_xlabel(r"$dE$")
    axes[1].set_ylabel("Intensity")
    axes[1].legend(frameon=False, fontsize=8)

    fig.savefig(ROOT / "psf2d_hcol_only_overlay.png", dpi=180)


def main() -> None:
    iq, ie = choose_representative_point()
    print(f"Representative point: iq={iq}, iE={ie}")

    write_selection_summary(iq, ie)
    run_exports(iq, ie)
    write_summary(iq, ie)
    plot_compare()
    plot_overlay()

    print("Done.")
    print("Generated:")
    print("  - psf2d_hcol_only_selection.txt")
    print("  - psf2d_hcol_only_summary.txt")
    print("  - psf2d_hcol_only_compare.png")
    print("  - psf2d_hcol_only_overlay.png")


if __name__ == "__main__":
    main()
