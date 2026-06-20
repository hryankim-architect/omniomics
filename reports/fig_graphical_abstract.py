#!/usr/bin/env python3
"""Graphical abstract for the anchored multi-omics framework — a self-contained schematic (no data needed).

Four stages, left to right: (1) anchor on established knowledge (a zero-parameter textbook prior), (2) admit
genome-wide data only as a non-negative gated residual (never below the anchor), (3) mine the residual to
discover an anchor-orthogonal axis, (4) read transportability — the same hypothesis can be NOVEL or REDUNDANT
depending on the measurement platform. Colour semantics: blue = known/anchor, green = newly discovered,
red = redundant. Writes reports/figs/graphical_abstract.png.

Run:  python reports/fig_graphical_abstract.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__)); FIGS = os.path.join(HERE, "figs"); os.makedirs(FIGS, exist_ok=True)
BLUE, GREEN, RED, GREY, INK, SUB = "#2c7fb8", "#2ca25f", "#d95f5f", "#9aa7b2", "#1a2b3c", "#33485c"

W, H, YB = 2.75, 2.95, 1.35          # box width, height, bottom y
YT = YB + H
XS = [0.25, 3.45, 6.65, 9.85]


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


def glabel(ax, cx, text):
    ax.text(cx, YB + 0.16, text, ha="center", va="center", fontsize=6.8, color=SUB, style="italic")


def main():
    fig, ax = plt.subplots(figsize=(12.8, 4.9))
    ax.set_xlim(0, 12.85); ax.set_ylim(0, 5.3); ax.axis("off")
    ax.text(6.4, 5.08, "Anchored multi-omics integration & knowledge-anchored residual discovery",
            ha="center", va="center", fontsize=13, fontweight="bold", color=INK)

    # 1 · anchor — a gene-signature glyph (several known-biology bars)
    cx = draw_box(ax, XS[0], BLUE, ["1 · Anchor on", "known biology"], ["start from a fixed", "textbook signature"])
    x = cx - 0.42
    for hh in [0.34, 0.52, 0.26, 0.44, 0.30]:
        ax.add_patch(plt.Rectangle((x, YB + 0.42), 0.13, hh, facecolor=BLUE, edgecolor="none")); x += 0.205
    glabel(ax, cx, "known-biology prior")

    # 2 · gate — anchor bar with a small green gain stacked ON TOP of a "never below" baseline
    cx = draw_box(ax, XS[1], "#5b6b7a", ["2 · Gate the", "residual (β ≥ 0)"], ["add data only if it", "beats the anchor"])
    base = YB + 0.42
    ax.add_patch(plt.Rectangle((cx - 0.17, base), 0.34, 0.50, facecolor=BLUE, edgecolor="none"))        # anchor
    ax.add_patch(plt.Rectangle((cx - 0.17, base + 0.50), 0.34, 0.16, facecolor=GREEN, edgecolor="none"))  # gain
    ax.plot([cx - 0.5, cx + 0.5], [base + 0.50, base + 0.50], ls="--", lw=1.1, color="#7c8b99")           # never-below line
    ax.annotate("", xy=(cx + 0.34, base + 0.66), xytext=(cx + 0.34, base + 0.50),
                arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.4))
    ax.text(cx + 0.46, base + 0.50, "never\nbelow", ha="left", va="center", fontsize=6.0, color="#7c8b99")
    glabel(ax, cx, "anchor, plus gain only")

    # 3 · discover — two perpendicular arrows: anchor axis (blue) + orthogonal new axis (green)
    cx = draw_box(ax, XS[2], BLUE, ["3 · Discover", "orthogonal axis"], ["mine what the", "anchor misses"])
    ox, oy = cx - 0.30, YB + 0.42
    ax.annotate("", xy=(ox + 0.62, oy), xytext=(ox, oy), arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.2))
    ax.annotate("", xy=(ox, oy + 0.62), xytext=(ox, oy), arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=2.2))
    ax.text(ox + 0.64, oy - 0.02, "anchor", ha="left", va="center", fontsize=6.2, color=BLUE)
    ax.text(ox - 0.04, oy + 0.70, "new axis", ha="center", va="bottom", fontsize=6.2, color=GREEN, fontweight="bold")
    glabel(ax, cx, "orthogonal to the anchor")

    # 4 · transportability — same effect, NOVEL (green) on one platform, REDUNDANT (red) on another
    cx = draw_box(ax, XS[3], RED, ["4 · Read", "transportability"], ["is the hit truly new", "or already explained?"])
    lx, rx = cx - 0.62, cx + 0.62
    ax.add_patch(plt.Rectangle((lx - 0.14, YB + 0.62), 0.28, 0.42, facecolor=GREEN, edgecolor="none"))
    ax.add_patch(plt.Rectangle((rx - 0.14, YB + 0.62), 0.28, 0.42, facecolor=RED, edgecolor="none"))
    ax.text(lx, YB + 0.50, "NOVEL", ha="center", va="center", fontsize=6.4, color=GREEN, fontweight="bold")
    ax.text(rx, YB + 0.50, "REDUNDANT", ha="center", va="center", fontsize=6.0, color=RED, fontweight="bold")
    ax.text(lx, YB + 1.10, "RNA-seq", ha="center", va="center", fontsize=6.0, color=SUB)
    ax.text(rx, YB + 1.10, "array", ha="center", va="center", fontsize=6.0, color=SUB)
    ax.annotate("", xy=(rx - 0.22, YB + 0.83), xytext=(lx + 0.22, YB + 0.83),
                arrowprops=dict(arrowstyle="<|-|>", color="#7c8b99", lw=1.2))
    glabel(ax, cx, "same effect, flips by platform")

    for i in range(3):
        ax.add_patch(FancyArrowPatch((XS[i] + W + 0.04, YB + H / 2), (XS[i + 1] - 0.04, YB + H / 2),
                                     arrowstyle="-|>", mutation_scale=16, linewidth=2, color="#8aa0b3"))

    # colour key + tagline
    key = [(BLUE, "known anchor"), (GREEN, "newly discovered"), (RED, "redundant / collinear")]
    kx = 6.4 - 2.7
    for c, lab in key:
        ax.add_patch(plt.Rectangle((kx, 0.78), 0.22, 0.22, facecolor=c, edgecolor="none"))
        ax.text(kx + 0.30, 0.89, lab, ha="left", va="center", fontsize=8, color=SUB); kx += 1.95
    ax.text(6.4, 0.34, "Never below the best single view — a discovery engine that finds what the known biology "
            "misses, and says what is truly new vs already explained.", ha="center", va="center",
            fontsize=9.2, style="italic", color=SUB)

    fig.savefig(os.path.join(FIGS, "graphical_abstract.png"), dpi=170, bbox_inches="tight", facecolor="white")
    print("wrote", os.path.join(FIGS, "graphical_abstract.png"))


if __name__ == "__main__":
    main()
