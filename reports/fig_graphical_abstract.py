#!/usr/bin/env python3
"""Graphical abstract for the anchored multi-omics framework — a self-contained schematic (no data needed).

Four stages, left to right: (1) anchor on established knowledge (a zero-parameter textbook prior), (2) admit
genome-wide data only as a non-negative gated residual (never below the anchor), (3) mine the residual to
discover an anchor-orthogonal axis, (4) read transportability — the same hypothesis can be NOVEL or REDUNDANT
depending on the measurement platform. Writes reports/figs/graphical_abstract.png.

Run:  python reports/fig_graphical_abstract.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__)); FIGS = os.path.join(HERE, "figs"); os.makedirs(FIGS, exist_ok=True)
BLUE, RED, GREY, INK = "#2c7fb8", "#d95f5f", "#bdbdbd", "#1a2b3c"


def box(ax, x, w, title, color, sub):
    ax.add_patch(FancyBboxPatch((x, 1.05), w, 2.5, boxstyle="round,pad=0.02,rounding_size=0.12",
                                linewidth=1.5, edgecolor=color, facecolor=color + "22"))
    ax.text(x + w / 2, 3.28, title, ha="center", va="center", fontsize=11.5, fontweight="bold", color=INK)
    for i, s in enumerate(sub):
        ax.text(x + w / 2, 2.95 - i * 0.32, s, ha="center", va="center", fontsize=8.6, color="#33485c")


def arrow(ax, x0, x1):
    ax.add_patch(FancyArrowPatch((x0, 2.3), (x1, 2.3), arrowstyle="-|>", mutation_scale=18,
                                 linewidth=2, color="#8aa0b3"))


def minibars(ax, cx, cy, heights, colors, w=0.12, gap=0.05, scale=0.9):
    n = len(heights); total = n * w + (n - 1) * gap; x = cx - total / 2
    for h, c in zip(heights, colors):
        ax.add_patch(plt.Rectangle((x, cy), w, h * scale, facecolor=c, edgecolor="none"))
        x += w + gap


def main():
    fig, ax = plt.subplots(figsize=(12.6, 4.4))
    ax.set_xlim(0, 20); ax.set_ylim(0, 4.7); ax.axis("off")

    ax.text(10, 4.45, "Anchored multi-omics integration & knowledge-anchored residual discovery",
            ha="center", va="center", fontsize=13.5, fontweight="bold", color=INK)

    box(ax, 0.3, 4.2, "1 · Anchor on known biology", BLUE,
        ["zero-parameter textbook prior", "(e.g. proliferation, ERBB2", "amplicon, TMB)"])
    minibars(ax, 2.4, 1.25, [0.55], [BLUE], w=1.6)
    box(ax, 5.2, 4.2, "2 · Gate the residual", "#555555",
        ["logit(anchor) + β·data,  β ≥ 0", "never below the best view;", "gate shut if redundant"])
    box(ax, 10.1, 4.2, "3 · Discover orthogonal axis", BLUE,
        ["mine what the anchor misses", "basal/keratinization; PD-L1;", "verified + reproduces"])
    box(ax, 15.0, 4.7, "4 · Read transportability", RED,
        ["NOVEL vs REDUNDANT vs INERT", "verdict set by anchor–hyp. corr,", "which can flip with the assay"])

    arrow(ax, 4.55, 5.15); arrow(ax, 9.45, 10.05); arrow(ax, 14.35, 14.95)

    # stage-4 mini motif: same effect, opposite label by platform
    minibars(ax, 16.4, 1.2, [0.5], [BLUE], w=0.5)      # RNA-seq -> NOVEL (blue)
    minibars(ax, 18.3, 1.2, [0.5], [RED], w=0.5)       # microarray -> REDUNDANT (red)
    ax.text(16.4, 1.0, "RNA-seq\nNOVEL", ha="center", va="top", fontsize=7.2, color=BLUE, fontweight="bold")
    ax.text(18.3, 1.0, "array\nREDUNDANT", ha="center", va="top", fontsize=7.2, color=RED, fontweight="bold")

    ax.text(10, 0.42, "Never below the best single view — a discovery engine that finds what the known biology misses, "
            "and says what is truly new vs already explained.", ha="center", va="center", fontsize=9.5,
            style="italic", color="#33485c")

    fig.savefig(os.path.join(FIGS, "graphical_abstract.png"), dpi=170, bbox_inches="tight", facecolor="white")
    print("wrote", os.path.join(FIGS, "graphical_abstract.png"))


if __name__ == "__main__":
    main()
