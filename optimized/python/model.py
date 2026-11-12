"""
Analytical Contact Mechanics Model for Fractional Kelvin-Voigt Nanoindentation.
"""

from typing import Union, Tuple, Dict, Any, Optional
import numpy as np
from scipy.special import gamma


def fractional_kv_force(
    t: Union[np.ndarray, list],
    E0: float,
    tau: float,
    alpha: float,
    nu: float = 2.0e-6,
    theta_deg: float = 20.0,
    Tr: Optional[float] = None,
) -> np.ndarray:
    """
    Computes theoretical force F(t) for AFM conical/pyramidal nanoindentation
    on a fractional Kelvin-Voigt viscoelastic material.

    Parameters
    ----------
    t : array-like
        Time vector in seconds.
    E0 : float
        Instantaneous / relaxed Young's elastic modulus (Pa). Must be > 0.
    tau : float
        Characteristic viscoelastic relaxation time (seconds). Must be > 0.
    alpha : float
        Fractional derivative order (0 <= alpha <= 1, dimensionless).
    nu : float, optional
        Indentation velocity in m/s (default: 2.0e-6 m/s = 2 um/s).
    theta_deg : float, optional
        Half-opening angle of conical indenter in degrees (default: 20.0).
    Tr : float, optional
        Loading duration / turnaround time in seconds (default: t[-1] / 2).

    Returns
    -------
    np.ndarray
        Predicted contact force vector in Newtons (N).
    """
    t = np.asarray(t, dtype=float)
    if Tr is None:
        Tr = float(t[-1]) / 2.0

    theta_rad = np.radians(theta_deg)
    KG = 1.0 / (2.0 * np.tan(theta_rad))
    P_const = 2.0 * KG * E0 * (nu ** 2)

    out = np.zeros_like(t)
    mask_loading = (t <= Tr)
    mask_unloading = ~mask_loading

    gamma3_alpha = gamma(3.0 - alpha)
    gamma2_alpha = gamma(2.0 - alpha)

    # 1. Loading Phase: 0 <= t <= Tr
    if np.any(mask_loading):
        t1 = t[mask_loading]
        with np.errstate(divide='ignore', invalid='ignore'):
            term_visco = np.where(
                t1 > 0,
                ((t1 / tau) ** (-alpha)) / gamma3_alpha,
                0.0
            )
        out[mask_loading] = P_const * (t1 ** 2) * (0.5 + term_visco)

    # 2. Unloading Phase: t > Tr
    if np.any(mask_unloading):
        t2 = t[mask_unloading]
        with np.errstate(divide='ignore', invalid='ignore'):
            term1 = (t2 ** 2 / gamma3_alpha) * ((t2 / tau) ** (-alpha))
            term2 = (2.0 * Tr * (t2 - Tr) / gamma2_alpha) * (((t2 - Tr) / tau) ** (-alpha))
            term3 = ((2.0 * Tr - t2) ** 2) / 2.0
            out[mask_unloading] = P_const * (term1 - term2 + term3)

    # Sanitize floating point exceptions
    out = np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)
    return out


def compute_goodness_of_fit(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    num_params: int = 3
) -> Dict[str, float]:
    """
    Computes statistical goodness-of-fit metrics.

    Parameters
    ----------
    y_true : np.ndarray
        Experimental observed values.
    y_pred : np.ndarray
        Model predicted values.
    num_params : int, optional
        Number of fitted parameters (k=3 for E0, tau, alpha).

    Returns
    -------
    Dict[str, float]
        Dictionary with R2, adjusted_R2, SSE, RMSE, and MAE.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    n = len(y_true)

    SSres = float(np.sum((y_true - y_pred) ** 2))
    mean_y = float(np.mean(y_true))
    SStot = float(np.sum((y_true - mean_y) ** 2))
    if SStot < 1e-30:
        SStot = 1.0

    R2 = 1.0 - (SSres / SStot)
    if n > num_params + 1:
        adj_R2 = 1.0 - ((1.0 - R2) * (n - 1)) / (n - num_params - 1)
    else:
        adj_R2 = R2

    rmse = float(np.sqrt(SSres / n))
    mae = float(np.mean(np.abs(y_true - y_pred)))

    return {
        "R_squared": R2,
        "adjusted_R_squared": adj_R2,
        "SSE": SSres,
        "RMSE": rmse,
        "MAE": mae,
        "n_samples": n
    }


class FractionalKVModel:
    """
    Object-oriented wrapper for Fractional Kelvin-Voigt contact mechanics.
    """
    def __init__(self, nu: float = 2.0e-6, theta_deg: float = 20.0):
        self.nu = nu
        self.theta_deg = theta_deg
        self.E0: Optional[float] = None
        self.tau: Optional[float] = None
        self.alpha: Optional[float] = None
        self.metrics: Dict[str, float] = {}

    def predict(self, t: np.ndarray, Tr: Optional[float] = None) -> np.ndarray:
        if self.E0 is None or self.tau is None or self.alpha is None:
            raise ValueError("Model parameters are not set. Fit the model or assign parameters first.")
        return fractional_kv_force(t, self.E0, self.tau, self.alpha, nu=self.nu, theta_deg=self.theta_deg, Tr=Tr)

    def evaluate(self, t: np.ndarray, y_true: np.ndarray) -> Dict[str, float]:
        y_pred = self.predict(t)
        self.metrics = compute_goodness_of_fit(y_true, y_pred)
        return self.metrics
