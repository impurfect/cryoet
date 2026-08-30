"""Shared plot styling. Import for side effects, use COLOR for consistency."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLOR = {
    "etomo": "#1f77b4", "aretomo": "#ff7f0e",
    "warp": "#1f77b4", "pytom": "#2ca02c",
    "both": "#9467bd", "grey": "#8c94a0",
}

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.5,
    "legend.frameon": False,
})


def save(fig, name):
    from config import OUT
    path = OUT / "plots" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    print(f"  {path.name}")
