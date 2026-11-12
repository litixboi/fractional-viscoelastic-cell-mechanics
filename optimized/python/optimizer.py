"""
High-Performance Parameter Optimizer for Fractional Kelvin-Voigt Model.
"""

from typing import Tuple, Dict, Any, Optional, List
import numpy as np
from scipy.optimize import minimize
from .model import fractional_kv_force, compute_goodness_of_fit


def fit_fractional_kv(
    t: np.ndarray,
    y: np.ndarray,
    nu: float = 2.0e-6,
    theta_deg: float = 20.0,
    p0: Optional[List[float]] = None,
    bounds: Optional[List[Tuple[float, float]]] = None,
    method: str = "Nelder-Mead",
    max_evals: int = 3000,
    multi_start: int = 1,
) -> Dict[str, Any]:
    """
    Fits constitutive parameters (E0, tau, alpha) to experimental AFM nanoindentation curves.

    Parameters
    ----------
    t : np.ndarray
        Time array in seconds.
    y : np.ndarray
        Experimental measured force in Newtons.
    nu : float, optional
        Indentation velocity in m/s (default: 2.0e-6).
    theta_deg : float, optional
        Half-opening indenter angle in degrees (default: 20.0).
    p0 : list of float, optional
        Initial parameter guess [E0, tau, alpha]. Default: [5000.0, 0.5, 0.05].
    bounds : list of (lower, upper), optional
        Parameter bounds. Default: [(0.0, None), (1e-8, None), (0.0, 1.0)].
    method : str, optional
        Optimization algorithm ('Nelder-Mead' or 'L-BFGS-B').
    max_evals : int, optional
        Maximum function evaluations.
    multi_start : int, optional
        Number of random starts around p0 to guarantee global convergence.

    Returns
    -------
    Dict[str, Any]
        Optimal parameters, metrics, predicted curve, and solver metadata.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    Tr = float(t[-1]) / 2.0

    if p0 is None:
        p0 = [5000.0, 0.5, 0.05]
    if bounds is None:
        bounds = [(0.0, None), (1e-8, None), (0.0, 1.0)]

    lb = np.array([b[0] if b[0] is not None else 0.0 for b in bounds])
    ub = np.array([b[1] if b[1] is not None else np.inf for b in bounds])

    # Objective function with exterior penalty
    def objective(p):
        p_arr = np.asarray(p, dtype=float)
        viol_lb = np.maximum(0.0, lb - p_arr)
        viol_ub = np.maximum(0.0, p_arr - ub)
        penalty = 1e6 * np.sum((viol_lb + viol_ub) ** 2)
        
        y_pred = fractional_kv_force(t, p_arr[0], p_arr[1], p_arr[2], nu=nu, theta_deg=theta_deg, Tr=Tr)
        sse = np.sum((y - y_pred) ** 2)
        return sse + penalty

    best_res = None
    best_cost = float("inf")

    # Multi-start strategy
    starts = [np.array(p0, dtype=float)]
    if multi_start > 1:
        np.random.seed(42)
        for _ in range(multi_start - 1):
            perturb = 1.0 + 0.5 * (np.random.rand(3) - 0.5)
            s = np.clip(starts[0] * perturb, lb + 1e-6, np.minimum(ub, 1e7))
            starts.append(s)

    for start_point in starts:
        if method == "Nelder-Mead":
            res = minimize(
                objective,
                start_point,
                method="Nelder-Mead",
                options={"maxiter": max_evals, "xatol": 1e-6, "fatol": 1e-20}
            )
        else:
            res = minimize(
                objective,
                start_point,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxfun": max_evals, "ftol": 1e-20, "gtol": 1e-12}
            )

        if res.fun < best_cost:
            best_cost = res.fun
            best_res = res

    p_opt = np.clip(best_res.x, lb, ub)
    E0_opt, tau_opt, alpha_opt = float(p_opt[0]), float(p_opt[1]), float(p_opt[2])
    y_pred = fractional_kv_force(t, E0_opt, tau_opt, alpha_opt, nu=nu, theta_deg=theta_deg, Tr=Tr)
    metrics = compute_goodness_of_fit(y, y_pred, num_params=3)

    return {
        "E0": E0_opt,
        "tau": tau_opt,
        "alpha": alpha_opt,
        "p_opt": [E0_opt, tau_opt, alpha_opt],
        "metrics": metrics,
        "R_squared": metrics["R_squared"],
        "adjusted_R_squared": metrics["adjusted_R_squared"],
        "SSE": metrics["SSE"],
        "RMSE": metrics["RMSE"],
        "y_pred": y_pred,
        "iterations": best_res.nit if hasattr(best_res, "nit") else 0,
        "nfev": best_res.nfev,
        "success": best_res.success
    }
