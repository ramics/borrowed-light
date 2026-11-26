"""Forest plot for the README Results section.

Inputs (aggregate numbers only, copied from the spike outputs):
  real_splits.csv   from spikes/power_rethink.py   (participant-level splits)
  loo_transfer.csv  from spikes/summary_pool.py    (leave-one-trial-out, published HRs)
Run: python docs/figures/make_results_figure.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).parent
INK, MUTED, ACCENT = "#191919", "#666666", "#257de0"

SPLIT_LABELS = {
    "png_village_group": "PNG P. vivax, Ilaita → other villages",
    "png_age_yr": "PNG P. vivax, age < 2.5 → ≥ 2.5 y",
    "hk_strain": "HK H3N2, Perth → Victoria epidemics",
    "hk_age.y": "HK H3N2, age < 40 → ≥ 40",
    "h505_age": "HVTN 505, age < 30 → ≥ 30",
    "h505_BMI": "HVTN 505, BMI < 25 → ≥ 25",
}
FAMILY_LABELS = {"covid_id50": "nAb ID50", "covid_spike": "spike IgG"}


def panel(ax, labels, est, se, highlight, xlabel, title):
    y = list(range(len(labels)))[::-1]
    for yi, e, s, h in zip(y, est, se, highlight):
        c = ACCENT if h else INK
        ax.plot([e - 1.96 * s, e + 1.96 * s], [yi, yi], color=c, lw=2, solid_capstyle="round")
        ax.plot(e, yi, "o", color=c, ms=8, mec="white", mew=2)
    ax.axvline(0, color=MUTED, lw=1, ls="--", zorder=0)
    ax.set_yticks(y, labels, color=INK)
    ax.set_xlabel(xlabel, color=INK)
    ax.set_title(title, loc="left", color=INK, fontsize=11, fontweight="bold")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", colors=MUTED)
    ax.grid(axis="x", color="#e5e5e5", lw=0.8, zorder=0)
    ax.set_axisbelow(True)


def main():
    rs = pd.read_csv(HERE / "real_splits.csv")
    loo = pd.read_csv(HERE / "loo_transfer.csv")
    loo_labels = [f"{FAMILY_LABELS[f]}, others → {h}" for f, h in zip(loo.family, loo.held_out)]

    fig, axes = plt.subplots(
        2, 1, figsize=(8, 7), gridspec_kw={"height_ratios": [len(rs), len(loo)], "hspace": 0.55}
    )
    panel(
        axes[0],
        [SPLIT_LABELS[s] for s in rs.split],
        rs.attenuation, rs.se_attenuation, rs.p_perm_attenuation < 0.05,
        "Attenuation of the marker's log-OR per SD (positive = weaker in target)",
        "A. Within-study splits, participant data",
    )
    panel(
        axes[1], loo_labels, loo.attenuation, loo.se_attenuation, [False] * len(loo),
        "Attenuation of the log-HR per 10-fold (positive = weaker in held-out trial)",
        "B. Between COVID-19 vaccine trials, published HRs",
    )
    fig.text(0.01, 0.005, "Points: estimate; bars: 95% CI (±1.96 SE). Blue: permutation p < 0.05.",
             color=MUTED, fontsize=8)
    fig.savefig(HERE / "results_forest.png", dpi=200, bbox_inches="tight", facecolor="white")


if __name__ == "__main__":
    main()
