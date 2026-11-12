"""
Unit Tests for Fractional Kelvin-Voigt Model Edge Cases and Mechanics.
"""

import pytest
import numpy as np
from optimized.python.model import fractional_kv_force, compute_goodness_of_fit, FractionalKVModel


def test_zero_time():
    """At t = 0, force should be exactly 0."""
    f = fractional_kv_force(np.array([0.0]), E0=5000.0, tau=0.01, alpha=0.3)
    assert np.isclose(f[0], 0.0), f"Force at t=0 should be 0, got {f[0]}"


def test_positive_force_monotonic_loading():
    """During loading, force should be monotonically increasing."""
    t = np.linspace(0.01, 1.0, 100)
    f = fractional_kv_force(t, E0=5000.0, tau=0.02, alpha=0.3, Tr=1.0)
    diff = np.diff(f)
    assert np.all(diff > 0), "Force during loading phase must be monotonically increasing"


def test_oo_model_wrapper():
    """Tests the object-oriented FractionalKVModel interface."""
    model = FractionalKVModel(nu=2.0e-6, theta_deg=20.0)
    with pytest.raises(ValueError):
        model.predict(np.array([0.1, 0.2]))

    model.E0 = 5000.0
    model.tau = 0.02
    model.alpha = 0.3
    t = np.linspace(0.01, 2.0, 50)
    pred = model.predict(t)
    assert len(pred) == 50
    assert np.all(pred[t <= 1.0] > 0), "Loading phase force must be strictly positive"

    # Test evaluate
    y_synth = pred * (1 + 0.01 * np.random.randn(len(pred)))
    metrics = model.evaluate(t, y_synth)
    assert metrics["R_squared"] > 0.99
