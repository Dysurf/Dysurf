from pathlib import Path
import re
import subprocess
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

#  q0 = (2.025, 0, 0)
#  E0 = 7.05 meV

ROOT = Path(__file__).resolve().parent
SRC = ROOT.parents[2] / "src" / "src_gfortran"
ACC = ROOT.parent / "2+H00_psf2d_minimal_acceptance"
PAT = re.compile(r"(?<=\d)([+-]\d{2,3})(?=\s|$)")


def load_matrix(path: Path) -> np.ndarray:
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            fixed = PAT.sub(r"E\1", line)
            vals = np.fromstring(fixed, sep=" ")
            if vals.size:
                rows.append(vals)
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

def choose_representative_point():
    shear_file = ACC / "PSF_shear_base_300K.dat"
    if not shear_file.exists():
        return 6, 142
    shear = load_matrix(shear_file)
    margin_q = 5
    margin_e = 10
    core = np.abs(shear[margin_q:-margin_q, margin_e:-margin_e])
    idx = np.unravel_index(np.argmax(core), core.shape)
    iq = idx[0] + margin_q + 1
    ie = idx[1] + margin_e + 1
    return iq, ie

def write_selection_summary(iq: int, ie: int) -> None:
    h0 = 2.0 + 0.005 * (iq - 1)
    e0 = 0.05 * (ie - 1)
    orient_file = ROOT / "instrument_base.txt"
    u = parse_vector_from_file(orient_file, "u")
    v = parse_vector_from_file(orient_file, "v")
    text = f"""Local PSF example selection
===========================
selection_mode = automatic
selection_rule = max_abs_shear_inside_margin
iq = {iq}
iE = {ie}
q0_hkl = [{h0:.6f}, 0.000000, 0.000000]
E0_meV = {e0:.6f}
preset = base
u = {u}
v = {v}

Instrument parameter file:
  instrument_base.txt

Nominal point definition file:
  example_nominal_point.txt
"""
    (ROOT / "psf2d_example_selection.txt").write_text(text, encoding="utf-8")

def load_triplets(path: Path):
    arr = load_matrix(path)
    dq = np.unique(arr[:, 0])
    de = np.unique(arr[:, 1])
    z = arr[:, 2].reshape(de.size, dq.size)
    return dq, de, z

def plot_example() -> None:
    dq_raw, de_raw, raw = load_triplets(ROOT / "psf2d_example_raw.dat")
    dq_fit, de_fit, fit = load_triplets(ROOT / "psf2d_example_fit.dat")
    dq_cut = load_matrix(ROOT / "psf2d_example_dq_cut.dat")
    de_cut = load_matrix(ROOT / "psf2d_example_dE_cut.dat")

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    im0 = axes[0, 0].imshow(raw, origin="lower", aspect="auto", cmap="magma",
                            extent=[dq_raw[0], dq_raw[-1], de_raw[0], de_raw[-1]])
    axes[0, 0].set_title("CN method")
    axes[0, 0].set_xlabel(r"$dQ$")
    axes[0, 0].set_ylabel(r"$dE$")
    plt.colorbar(im0, ax=axes[0, 0], fraction=0.046, pad=0.04)

    im1 = axes[0, 1].imshow(fit, origin="lower", aspect="auto", cmap="magma",
                            extent=[dq_fit[0], dq_fit[-1], de_fit[0], de_fit[-1]])
    axes[0, 1].set_title("Fitted sheared 2D PSF")
    axes[0, 1].set_xlabel(r"$dQ$")
    axes[0, 1].set_ylabel(r"$dE$")
    plt.colorbar(im1, ax=axes[0, 1], fraction=0.046, pad=0.04)

    axes[1, 0].plot(dq_cut[:, 0], dq_cut[:, 1], marker="o",ls=" ", label="CN method")
    axes[1, 0].plot(dq_cut[:, 0], dq_cut[:, 2], lw=1.8, label="fit")
    axes[1, 0].set_title(r"Momentum resolution")
    axes[1, 0].set_xlabel(r"$dQ$")
    axes[1, 0].set_ylabel("Intensity")
    axes[1, 0].legend(frameon=False)

    axes[1, 1].plot(de_cut[:, 0], de_cut[:, 1], ls=" ", marker="o", label="CN method")
    axes[1, 1].plot(de_cut[:, 0], de_cut[:, 2], lw=1.8, label="fit")
    axes[1, 1].set_title(r"Energy resolution")
    axes[1, 1].set_xlabel(r"$dE$")
    axes[1, 1].set_ylabel("Intensity")
    axes[1, 1].legend(frameon=False)

    fig.savefig(ROOT / "psf2d_example.png", dpi=180)

def main() -> None:
    iq, ie = choose_representative_point()
    write_selection_summary(iq, ie)
    exe = SRC / "validate_psf2d_parametrization"
    if not exe.exists():
        raise SystemExit("validate_psf2d_parametrization executable not found. Build it first.")
    orient_file = ROOT / "instrument_base.txt"
    u = parse_vector_from_file(orient_file, "u")
    v = parse_vector_from_file(orient_file, "v")
    cmd = [str(exe), str(iq), str(ie), "base", "psf2d_example"] + [str(x) for x in (u + v)]
    subprocess.run(cmd, cwd=ROOT, check=True)
    plot_example()

if __name__ == "__main__":
    main()
