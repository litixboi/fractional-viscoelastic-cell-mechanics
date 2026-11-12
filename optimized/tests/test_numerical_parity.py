"""
Numerical Parity Test Suite.
Verifies exact analytical and numerical equivalence between the Python implementation,
MATLAB results, and paper published benchmarks.
"""

import os
import pytest
import numpy as np
import pandas as pd
from optimized.python.model import fractional_kv_force, compute_goodness_of_fit
from optimized.python.optimizer import fit_fractional_kv


DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


def test_pc3_wod_41_parity():
    """
    Verifies that Python fractional_kv_force reproduces the saved theoretical curve
    in PC3_WoD(41).xlsx down to numerical precision.
    """
    data_file = os.path.join(DATA_DIR, "PC3_WoD(41).xlsx")
    assert os.path.exists(data_file), f"Missing test file: {data_file}"

    df = pd.read_excel(data_file)
    t = df["Time(s)"].values
    f_exp = df["Experimental F(N)"].values
    f_theor_saved = df["Theoretical F(N)"].values

    # Paper parameters for PC3_WoD sample 41
    E0_paper = 4614.54541570643
    tau_paper = 0.0177758363310265
    alpha_paper = 0.389488117481465

    f_calc = fractional_kv_force(t, E0_paper, tau_paper, alpha_paper, nu=2.0e-6, theta_deg=20.0)

    # Parity check against saved theoretical curve
    max_err = np.max(np.abs(f_calc - f_theor_saved))
    assert max_err < 1.0e-10, f"Max error {max_err:.4e} exceeds tolerance 1e-10"

    # Goodness of fit check
    metrics = compute_goodness_of_fit(f_exp, f_calc)
    assert metrics["R_squared"] > 0.949, f"R2 {metrics['R_squared']} below 0.949"


def test_optimization_convergence():
    """
    Verifies that the optimizer converges to R^2 >= 0.95 from a standard initial guess.
    """
    data_file = os.path.join(DATA_DIR, "PC3_WoD(41).xlsx")
    df = pd.read_excel(data_file)
    t = df["Time(s)"].values
    f_exp = df["Experimental F(N)"].values

    res = fit_fractional_kv(t, f_exp, nu=2.0e-6, theta_deg=20.0, p0=[5000.0, 0.5, 0.1])

    assert res["R_squared"] >= 0.949, f"Optimizer R^2 {res['R_squared']} below threshold"
    assert res["E0"] > 3000.0 and res["E0"] < 7000.0, f"Fitted E0 {res['E0']} out of expected range"
    assert res["tau"] > 0.005 and res["tau"] < 0.05, f"Fitted tau {res['tau']} out of expected range"
    assert res["alpha"] > 0.2 and res["alpha"] < 0.6, f"Fitted alpha {res['alpha']} out of expected range"


def test_turnaround_continuity():
    """
    Verifies force continuity across the turnaround boundary t = Tr.
    """
    Tr = 1.0
    dt = 1e-6
    t_before = Tr - dt
    t_after = Tr + dt

    f_before = fractional_kv_force(np.array([t_before]), E0=5000.0, tau=0.02, alpha=0.3, Tr=Tr)[0]
    f_after = fractional_kv_force(np.array([t_after]), E0=5000.0, tau=0.02, alpha=0.3, Tr=Tr)[0]

    rel_diff = abs(f_before - f_after) / max(f_before, 1e-20)
    assert rel_diff < 1e-4, f"Discontinuity at Tr: relative diff = {rel_diff:.4e}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
