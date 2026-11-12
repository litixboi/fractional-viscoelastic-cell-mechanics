"""
Publication-Quality Visualization for AFM Nanoindentation & Viscoelastic Fits.
"""

from typing import Optional
import os
import numpy as np
import matplotlib.pyplot as plt


def plot_force_fit(
    t: np.ndarray,
    y_exp: np.ndarray,
    y_fit: np.ndarray,
    title: str = "AFM Nanoindentation & Fractional Kelvin-Voigt Fit",
    cell_condition: str = "PC-3 Cancer Cell",
    R2: Optional[float] = None,
    save_path: Optional[str] = None,
    show: bool = False
) -> plt.Figure:
    """
    Renders a dual-panel publication figure showing experimental force,
    theoretical fit, turnaround point, and fitting residuals.
    """
    t = np.asarray(t, dtype=float)
    y_exp_nN = np.asarray(y_exp, dtype=float) * 1e9  # Convert N to nN
    y_fit_nN = np.asarray(y_fit, dtype=float) * 1e9  # Convert N to nN
    residuals_nN = y_exp_nN - y_fit_nN
    Tr = t[-1] / 2.0

    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        figsize=(8, 6),
        dpi=300,
        gridspec_kw={'height_ratios': [2.5, 1.0]},
        sharex=True
    )

    # Panel 1: Primary fit
    ax1.scatter(t, y_exp_nN, s=12, color='#1f77b4', alpha=0.6, label='Experimental Force (AFM)', zorder=2)
    fit_label = f'Theoretical Fit ($R^2 = {R2:.4f}$)' if R2 is not None else 'Theoretical Fit'
    ax1.plot(t, y_fit_nN, color='#d62728', linewidth=2.0, label=fit_label, zorder=3)
    ax1.axvline(Tr, color='black', linestyle='--', linewidth=1.2, alpha=0.7, label=f'Turnaround $T_r = {Tr:.2f}\\,\\mathrm{{s}}$')

    ax1.set_ylabel('Force (nN)', fontsize=12, fontweight='bold')
    ax1.set_title(f"{title} - {cell_condition}", fontsize=13, fontweight='bold', pad=10)
    ax1.legend(loc='upper left', frameon=True, framealpha=0.9, fontsize=10)
    ax1.grid(True, linestyle=':', alpha=0.6)

    # Panel 2: Residuals
    ax2.plot(t, residuals_nN, color='#2ca02c', linewidth=1.0)
    ax2.axhline(0, color='red', linestyle='--', linewidth=1.0, alpha=0.8)
    ax2.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Residuals (nN)', fontsize=10, fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure successfully saved to: {save_path}")

    if show:
        plt.show()

    return fig
