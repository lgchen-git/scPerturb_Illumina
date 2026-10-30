"""Shared Nature/NPG-inspired figure palette and restrained Matplotlib defaults."""

NPG_CATEGORICAL = (
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
    "#7E6148", "#B09C85", "#6F99AD", "#8A7AA9",
    "#C9A227", "#6E7F80", "#D4537E", "#378ADD",
)
NPG_EVIDENCE = {
    "LOW_CELL": "#8491B4", "LOW_GUIDE_SUPPORT": "#F39B7F",
    "GUIDE_DOMINATED": "#7E6148", "DISCORDANT": "#3C5488",
    "PARTIAL_SUPPORT": "#00A087", "SUPPORTED": "#E64B35",
    "INCONCLUSIVE": "#B09C85",
}
NPG_SIGNED = ("#3C5488", "#F7F7F7", "#E64B35")
NPG_ENRICHMENT = ("#E5E7EB", "#91D1C2", "#00A087", "#DC0000")

def apply_matplotlib_style() -> None:
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titleweight": "regular",
        "axes.labelcolor": "#333333", "text.color": "#222222",
        "axes.edgecolor": "#555555", "axes.facecolor": "white",
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "axes.grid": False, "grid.color": "#E5E7EB",
    })
