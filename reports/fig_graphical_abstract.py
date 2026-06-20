#!/usr/bin/env python3
"""Graphical abstract for the anchored multi-omics framework — a self-contained schematic (no data needed).

Four stages, left to right: (1) anchor on established knowledge (a zero-parameter textbook prior), (2) admit
genome-wide data only as a non-negative gated residual (never below the anchor), (3) mine the residual to
discover an anchor-orthogonal axis, (4) read transportability — the same hypothesis can be NOVEL or REDUNDANT
depending on the measurement platform. Writes reports/figs/graphical_abstract.png.

Run:  python reports/fig_graphical_abstract.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__)); FIGS = os.path.join(HERE, "figs"); os.makedirs(FIGS, exist_ok=True)
BLUE, RED, GREY, INK, SUB = "#2c7fb8", "#d95f5f", "#9aa7b2", "#1a2b3c", "#33485c"

W, H, YB = 2.75, 2.95, 1.15          # box width, height, bottom y
YT = YB + H                           # box top
XS = [0.25, 3.45, 6.65, 9.85]         # box left edges (gap 0.45 for arrows)


def draw_box(ax, x, color, title_lines, sublines):
    cx = x + W / 2
    ax.add_patch(FancyBboxPatch((x, YB), W, H, boxstyle="round,pad=0.015,rounding_size=0.10",
                                linewidth=1.6, edgecolor=color, facecolor=color + "1e"))
    for i, t in enumerate(title_lines):
        ax.text(cx, YT - 0.40 - i * 0.40, t, ha="center", va="center", fontsize=11, fontweight="bold", color=INK)
    y0 = YT - 0.40 - len(title_lines) * 0.40 - 0.18
    for i, s in enumerate(sublines):
        ax.text(cx, y0 - i * 0.30, s, ha="center", va="center", fontsize=7.7, color=SUB)
    return cx


def bars(ax, cx, base, heights, colors, bw=0.16, gap=0.07):
    n = len(heights); total = n * bw + (n - 1) * gap; x = cx - total / 2
    for h, c in zip(heights, colors):
        ax.add_patch(plt.Rectangle((x, base), bw, h, facecolor=c, edgecolor="none")); x += bw + gap


def arrow(ax, x0, x1):
    ax.add_patch(FancyArrowPatch((x0, YB + H / 2), (x1, YB + H / 2), arrowstyle="-|>",
                                 mutation_scale=16, linewidth=2, color="#8aa0b3"))


def main():
    fig, ax = plt.subplots(figsize=(12.8, 4.6))
    ax.set_xlim(0, 12.85); ax.set_ylim(0, 5.0); ax.axis("off")
    ax.text(6.4, 4.78, "Anchored multi-omics integration & knowledge-anchored residual discovery",
            ha="center", va="center", fontsize=13, fontweight="bold", color=INK)

    # 1 · anchor — a gene-signature glyph (several bars)
    cx = draw_box(ax, XS[0], BLUE, ["1 · Anchor on", "known biology"],
                  ["zero-parameter", "textbook prior"])
    bars(ax, cx, YB + 0.28, [0.40, 0.62, 0.30, 0.52, 0.36], [BLUE] * 5)

    # 2 · gate — anchor bar (tall) + small gated residual on top
    cx = draw_box(ax, XS[1], "#5b6b7a", ["2 · Gate the", "residual (β ≥ 0)"],
                  ["never below the", "best single view"])
    bars(ax, cx, YB + 0.28, [0.70, 0.16], ["#5b6b7a", GREY], bw=0.34, gap=0.10)

    # 3 · discover — anchor bar + orthogonal increment highlighted
    cx = draw_box(ax, XS[2], BLUE, ["3 · Discover", "orthogonal axis"],
                  ["mine what the", "anchor misses"])
    ax.add_patch(plt.Rectangle((cx - 0.17, YB + 0.28), 0.34, 0.55, facecolor=GREY, edgecolor="none"))
    ax.add_patch(plt.Rectangle((cx - 0.17, YB + 0.83), 0.34, 0.22, facecolor=BLUE, edgecolor="none"))

    # 4 · transportability — same effect, opposite label by platform (two labelled bars, inside the box)
    cx = draw_box(ax, XS[3], RED, ["4 · Read", "transportability"],
                  ["NOVEL / REDUNDANT,", "set by the assay"])
    bx = cx - 0.55
    ax.add_patch(plt.Rectangle((bx - 0.13, YB + 0.40), 0.26, 0.45, facecolor=BLUE, edgecolor="none"))
    ax.add_patch(plt.Rectangle((cx + 0.55 - 0.13, YB + 0.40), 0.26, 0.45, facecolor=RED, edgecolor="none"))
    ax.text(bx, YB + 0.26, "RNA-seq", ha="center", va="center", fontsize=6.6, color=BLUE, fontweight="bold")
    ax.text(cx + 0.55, YB + 0.26, "array", ha="center", va="center", fontsize=6.6, color=RED, fontweight="bold")

    for i in range(3):
        arrow(ax, XS[i] + W + 0.04, XS[i + 1] - 0.04)

    ax.text(6.4, 0.42, "Never below the best single view — a discovery engine that finds what the known biology "
            "misses, and says what is truly new vs already explained.", ha="center", va="center",
            fontsize=9.2, style="italic", color=SUB)

    fig.savefig(os.path.join(FIGS, "graphical_abstract.png"), dpi=170, bbox_inches="tight", facecolor="white")
    print("wrote", os.path.join(FIGS, "graphical_abstract.png"))


if __name__ == "__main__":
    main()
