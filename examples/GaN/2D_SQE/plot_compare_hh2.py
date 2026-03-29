import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams["font.family"] = "Times New Roman"

SCALE = 15.957691216057308
EXTENT = [-3.0, -1.0, 0, 50]
CLIM = (3, 9)


def load_logged_data(path):
    data = np.loadtxt(path) * SCALE
    return np.log(data)


def decorate_axis(ax, label):
    ax.text(-3.1, 51, r"$\mathregular{\Gamma}$", fontsize=13)
    ax.text(-2.6, 51, "M", fontsize=13)
    ax.text(-2.1, 51, r"$\mathregular{\Gamma}$", fontsize=13)
    ax.text(-1.6, 51, "M", fontsize=13)
    ax.text(-1.1, 51, r"$\mathregular{\Gamma}$", fontsize=13)
    ax.text(-2.95, 46.5, label, fontsize=13, color="white")
    ax.set_yticks([])


data_300k = load_logged_data("h-h2/SQE_300K.dat")
data_14k = load_logged_data("h-h2_14K/SQE_14K.dat")

fig, axes = plt.subplots(
    2, 1, figsize=(5.0, 5.6), sharex=True, constrained_layout=True
)

images = []
for ax, data, label in zip(axes, [data_300k, data_14k], ["300K", "14K"]):
    image = ax.imshow(
        np.rot90(data),
        interpolation="bicubic",
        extent=EXTENT,
        aspect=0.05,
        cmap="jet",
        vmin=CLIM[0],
        vmax=CLIM[1],
    )
    decorate_axis(ax, label)
    images.append(image)

axes[1].set_xlabel("[H-H2] (r.l.u.)", fontsize=13)
for ax in axes:
    ax.set_ylabel("E (meV)", fontsize=13)

fig.colorbar(images[0], ax=axes, shrink=0.95, pad=0.02, label="log S(Q,E)")
fig.savefig("h-h2_compare_14K_300K.png", dpi=300)
