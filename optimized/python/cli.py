"""
Command-Line Interface (CLI) for Fractional Kelvin-Voigt Fitting.
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np

from .model import fractional_kv_force
from .optimizer import fit_fractional_kv
from .plot import plot_force_fit


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit Fractional Kelvin-Voigt contact mechanics model to AFM nanoindentation data."
    )
    parser.add_argument(
        "--file", "-f", required=True,
        help="Path to Excel (.xlsx) or CSV file containing experimental data (Col 1: Time(s), Col 2: Force(N))."
    )
    parser.add_argument(
        "--velocity", "-v", type=float, default=2.0e-6,
        help="Indentation velocity in m/s (default: 2e-6)."
    )
    parser.add_argument(
        "--theta", "-t", type=float, default=20.0,
        help="Half-opening angle of indenter in degrees (default: 20.0)."
    )
    parser.add_argument(
        "--p0", nargs=3, type=float, default=[5000.0, 0.5, 0.05],
        help="Initial parameter guess: E0 (Pa), tau (s), alpha (dimensionless)."
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="Path to save output results CSV or Excel."
    )
    parser.add_argument(
        "--plot", "-p", default=None,
        help="Path to save fit plot image (e.g. results/fit.png)."
    )
    parser.add_argument(
        "--title", default="AFM Nanoindentation Fit",
        help="Title for the output figure."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.file):
        print(f"Error: File not found: {args.file}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading data from: {args.file} ...")
    if args.file.endswith('.csv'):
        df = pd.read_csv(args.file)
    else:
        df = pd.read_excel(args.file)

    t = df.iloc[:, 0].values
    y = df.iloc[:, 1].values

    print(f"Running optimization with velocity = {args.velocity:.2e} m/s, theta = {args.theta} deg ...")
    res = fit_fractional_kv(t, y, nu=args.velocity, theta_deg=args.theta, p0=args.p0)

    print("\n=== Optimization Results ===")
    print(f"  E0 (Elastic Modulus)     : {res['E0']:.4f} Pa")
    print(f"  tau (Relaxation Time)    : {res['tau']:.6f} s")
    print(f"  alpha (Fractional Order) : {res['alpha']:.6f}")
    print(f"  R-squared                : {res['R_squared']:.6f}")
    print(f"  Adjusted R-squared       : {res['adjusted_R_squared']:.6f}")
    print(f"  SSE                      : {res['SSE']:.4e} N^2")
    print(f"  RMSE                     : {res['RMSE']:.4e} N")
    print(f"  Function Evaluations     : {res['nfev']}")

    out_df = pd.DataFrame({
        "Time(s)": t,
        "Experimental_Force(N)": y,
        "Theoretical_Force(N)": res["y_pred"],
        "Residuals(N)": y - res["y_pred"]
    })

    if args.output:
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        if args.output.endswith('.csv'):
            out_df.to_csv(args.output, index=False)
        else:
            out_df.to_excel(args.output, index=False)
        print(f"Fitted data saved to: {args.output}")

    if args.plot:
        plot_force_fit(
            t, y, res["y_pred"],
            title=args.title,
            save_path=args.plot
        )


if __name__ == "__main__":
    main()
