"""
Publication-Quality Visualization for AFM Nanoindentation & Viscoelastic Fits.
Matches the official paper publication style: single-panel plot with experimental data
and analytical fitted curve without residual subplots.
"""

from typing import Optional
import os
import numpy as np
import matplotlib.pyplot as plt


def plot_force_fit(
    t: np.ndarray,
    y_exp: np.ndarray,
    y_fit: np.ndarray,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    show: bool = False
) -> plt.Figure:
    """
    Renders a publication figure matching the paper's exact visual style:
    single-panel force-versus-time curve with experimental scatter and fitted line.
    """
    t = np.asarray(t, dtype=float)
    y_exp = np.asarray(y_exp, dtype=float)
    y_fit = np.asarray(y_fit, dtype=float)

    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=300)

    # Plot experimental data (blue dots) and fitted curve (solid red line)
    ax.plot(t, y_exp, '.', color='#00008b', markersize=4, label='Experimental data', zorder=2)
    ax.plot(t, y_fit, '-', color='#c00000', linewidth=2.0, label='Fitted curve', zorder=3)

    ax.set_xlabel('time (s)', fontsize=12)
    ax.set_ylabel('Force (N)', fontsize=12)

    # Format scientific notation on Y axis
    ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
    ax.yaxis.get_offset_text().set_fontsize(11)

    ax.legend(loc='upper right', frameon=True, edgecolor='#333333', fontsize=11)
    ax.set_xlim(left=0, right=float(t[-1]))

    if title:
        ax.set_title(title, fontsize=12, pad=10)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure successfully saved to: {save_path}")

    if show:
        plt.show()

    return fig
