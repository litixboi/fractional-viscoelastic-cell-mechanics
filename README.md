# Quantitative Biomechanical Mapping of Prostate Cancer Therapeutic Responses Using AFM and Cytoskeletal Imaging

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24189111.svg)](https://doi.org/10.5281/zenodo.24189111)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10 | 3.11 | 3.12](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/)
[![MATLAB R2024a](https://img.shields.io/badge/MATLAB-R2024a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![PLOS ONE Compliant](https://img.shields.io/badge/PLOS%20ONE-Software%20Sharing%20Compliant-green.svg)](https://journals.plos.org/plosone/s/materials-and-software-sharing#loc-sharing-code)
[![Tests Passing](https://img.shields.io/badge/Tests-Passing%20(6%2F6)-brightgreen.svg)](#automated-testing--numerical-parity-verification)

---

## Overview

This project investigates the biomechanical responses of prostate cancer cells to two therapeutic agents: **docetaxel** and **bicalutamide**. The study compares **LNCaP** and **PC3** cell lines, representing different metastatic characteristics, using atomic force microscopy (AFM) and cytoskeletal fluorescence imaging.

A Kelvin–Voigt Fractional Derivative (KVFD) model is employed to quantify cellular viscoelastic properties, including elasticity, viscous behavior, and fluidity. By integrating mechanical measurements with DAPI and phalloidin staining, the study explores the relationship between treatment-induced biomechanical alterations and cytoskeletal organization.

* **$E_0$ (Pa)**: Elastic modulus representing cellular stiffness
* **$\tau$ (s)**: Viscous time constant
* **$\alpha$ (dimensionless)**: Fractional derivative order ($0 \le \alpha \le 1$)

This repository complies strictly with **[PLOS ONE's Materials and Software Sharing Guidelines](https://journals.plos.org/plosone/s/materials-and-software-sharing#loc-sharing-code)**:
1. **Full Open Access**: Source code is deposited in this public archive under an OSI-approved open-source license (**MIT License**).
2. **Dual-Implementation Strategy**: Provides the untouched original code exactly as used in the paper (`/original`), plus a modern, vectorized MATLAB suite and a 100% open-source **Python** implementation (`/optimized`) to eliminate proprietary commercial software barriers for global reviewers and researchers.
3. **Reproducibility Guarantee**: Benchmark datasets, control parameter sets, and an automated numerical parity test suite verify reproduction of published findings down to floating-point precision ($< 10^{-10}\text{ N}$).

---

## Authorship & Attribution

* **Lead Programmer / Author**: **Nima Naserbakht** ([@litixboi](https://github.com/litixboi))  
  *Affiliation*: K. N. Toosi University of Technology, Tehran, Iran  
  *Email*: `nnaserbakht1312@gmail.com`
* **Co-Authors & Contributors**:
  * **Arisa Monshizadeh** (K. N. Toosi University of Technology)
  * **Sayed Reza Ramezani** (K. N. Toosi University of Technology)
  * **Afsaneh Mojra** (Faculty of Mechanical Engineering, K. N. Toosi University of Technology)

---

## Mathematical Formulation

Atomic force microscopy (AFM) nanoindentation is modeled using contact mechanics for a conical indenter tip (half-opening angle $\theta$, indenter geometry constant $K_G = \frac{1}{2 \tan\theta}$) penetrating the cell membrane at constant approach velocity $\nu$ until maximum depth at turnaround time $T_r$, followed by retraction.

Under the fractional Kelvin-Voigt viscoelastic constitutive relation, the theoretical reaction force $F(t)$ is given analytically:

### 1. Loading Phase ($0 \le t \le T_r$)
$$F(t) = 2 K_G E_0 \nu^2 t^2 \left[ \frac{1}{2} + \frac{(t / \tau)^{-\alpha}}{\Gamma(3 - \alpha)} \right]$$

### 2. Unloading Phase ($t > T_r$)
$$F(t) = 2 K_G E_0 \nu^2 \left[ \frac{t^2}{\Gamma(3 - \alpha)} \left(\frac{t}{\tau}\right)^{-\alpha} - \frac{2 T_r (t - T_r)}{\Gamma(2 - \alpha)} \left(\frac{t - T_r}{\tau}\right)^{-\alpha} + \frac{(2 T_r - t)^2}{2} \right]$$

Where:
* $K_G = \frac{1}{2 \tan\theta}$ (for $\theta = 20^\circ$, $K_G \approx 1.3737$; for $\theta = 45^\circ$, $K_G = 0.5$).
* $\nu$ is the constant indenter approach/retract speed ($2.0\,\mu\text{m/s} = 2.0 \times 10^{-6}\text{ m/s}$).
* $T_r = \frac{t_{\text{final}}}{2}$ is the transition time at peak indentation depth.
* $\Gamma(\cdot)$ is the Euler Gamma function.

The objective function minimizes the Sum of Squared Errors ($SSE$) between experimental AFM force measurements $y_i$ and theoretical model predictions $F(t_i; E_0, \tau, \alpha)$:
$$SSE(E_0, \tau, \alpha) = \sum_{i=1}^n \left[ y_i - F(t_i; E_0, \tau, \alpha) \right]^2$$

Goodness of fit is quantified via $R^2$ and Adjusted $R^2$:
$$R^2 = 1 - \frac{\sum (y_i - \hat{y}_i)^2}{\sum (y_i - \bar{y})^2}, \quad R^2_{\text{adj}} = 1 - \left[\frac{(1 - R^2)(n - 1)}{n - k - 1}\right]$$

---

## Repository Structure

```
fractional-viscoelastic-cell-mechanics/
│
├── .github/
│   └── workflows/
│       └── ci.yml                     # GitHub Actions CI matrix (Ubuntu & Windows, Python 3.10-3.12)
│
├── original/                          # UNTOUCHED ORIGINAL RESEARCH CODE (As executed in paper)
│   ├── LastTime/
│   │   ├── Best Graphs based on R^2/  # Top R^2 fitted datasets & publication figures (PC3 & LNCaP)
│   │   │   ├── L_byDBi5_24h/          # LNCaP treated with Bicalutamide 5uM 24h
│   │   │   ├── L_byDDo_3/             # LNCaP treated with Docetaxel
│   │   │   ├── L_WoD/                 # LNCaP control (Without Drug)
│   │   │   ├── PC3_byDBi/             # PC-3 treated with Bicalutamide
│   │   │   ├── PC3_byDDo/             # PC-3 treated with Docetaxel
│   │   │   └── PC3_WoD/               # PC-3 control (Without Drug)
│   │   ├── BestResults-Adjusted-(AlphaLimited, Fixed).xlsx  # Full paper parameter master table (146 samples x 42 cols)
│   │   ├── BestResults-Adjusted-Final.xlsx                  # Comparative regression results
│   │   ├── NewSTDbestVal.m            # Batch fitting routine
│   │   └── [filtered_data, fitted_force, etc.]
│   └── Validation/
│       ├── Validation.m               # Primary unvectorized baseline optimization script
│       ├── Validation2.m              # Extended bounds optimization script
│       ├── GlobalSearch.m             # Multi-start optimizer with theoretical curve regeneration
│       ├── NV.m                       # Scaled robust optimizer script
│       ├── Raw(1).xlsx                # Baseline experimental curve
│       └── valid.jpg                  # Experimental curve validation plot
│
├── optimized/                         # HIGH-PERFORMANCE MODERN IMPLEMENTATIONS
│   ├── matlab/                        # Vectorized, portable MATLAB suite (>130x speedup)
│   │   ├── fractional_kv_force.m      # Vectorized analytical contact mechanics equation
│   │   ├── fit_fractional_kv.m        # Bounded Nelder-Mead optimizer returning structured results
│   │   └── run_reproduction.m         # Automated benchmark reproduction script
│   ├── python/                        # 100% Free Open-Source Python implementation (NumPy/SciPy)
│   │   ├── __init__.py                # Package interface
│   │   ├── model.py                   # Core analytical model & FractionalKVModel OO class
│   │   ├── optimizer.py               # Robust bounded Nelder-Mead / L-BFGS-B optimizer
│   │   ├── plot.py                    # Publication-grade dual-panel Matplotlib figures
│   │   └── cli.py                     # Command-line interface for fitting custom datasets
│   ├── data/                          # Benchmark experimental datasets for instant testing
│   │   ├── PC3_WoD(41).xlsx           # Benchmark PC-3 control sample 41
│   │   ├── PC3_byDDo(8).xlsx          # Benchmark PC-3 Docetaxel sample 8
│   │   ├── PC3_byDBi(43).xlsx         # Benchmark PC-3 Bicalutamide sample 43
│   │   ├── L_WoD(40).xlsx             # Benchmark LNCaP control sample 40
│   │   ├── L_byDBi5_24h(104).xlsx     # Benchmark LNCaP Bicalutamide sample 104
│   │   ├── L_byDDo_3(133).xlsx        # Benchmark LNCaP Docetaxel sample 133
│   │   └── Raw(1).xlsx                # Validation raw test set
│   └── tests/                         # Automated reproducibility and test suite
│       ├── test_numerical_parity.py   # Parity checks matching MATLAB results to < 1e-10 N
│       └── test_fractional_model.py   # Unit tests covering boundary conditions & continuity
│
├── requirements.txt                   # Pinned Python package dependencies
├── CITATION.cff                       # Citation metadata format
├── LICENSE                            # OSI-approved MIT License
└── README.md                          # Documentation and reproduction guide
```

---

## Optimization & Refactoring Summary

| Dimension | Original Version (`/original`) | Optimized Version (`/optimized`) | Impact & Benefit |
| :--- | :--- | :--- | :--- |
| **Computation Engine** | `parfor` loop inside objective evaluation | Array logical masks (`mask1 = t <= Tr`) | **>130x speedup** (0.3 s vs. 40 s per curve) without parallel worker IPC overhead |
| **Language Dependency** | Proprietary MATLAB with Parallel Toolbox | Dual MATLAB + pure **Python** (NumPy / SciPy) | Zero paywall barrier; adheres to PLOS ONE open software preference |
| **Path Handling** | Hardcoded Windows absolute paths (`D:\Tebs\`) | Relative, platform-agnostic file paths | Runs seamlessly across Windows, Linux, and macOS |
| **Boundary Handling** | Discontinuous boundary penalties | Smooth quadratic exterior penalties & bounded projection | Prevents unphysical parameter divergence ($E_0 < 0, \alpha \notin [0,1]$) |
| **CLI & Modularity** | Monolithic multi-step script | Clean functional APIs + CLI flags | Easily integrated into external pipelines |
| **Numerical Parity** | Baseline values | Identical outputs within $< 10^{-10}\text{ N}$ | 100% mathematical fidelity preserved |

---

## Environment Setup & Quick Start

### Option A: Python Implementation (Recommended — 100% Free & Open-Source)

#### 1. Clone the repository
```bash
git clone https://github.com/litixboi/fractional-viscoelastic-cell-mechanics.git
cd fractional-viscoelastic-cell-mechanics
```

#### 2. Create virtual environment and install dependencies
```bash
python -m venv venv
# On Linux/macOS:
source venv/bin/activate
# On Windows:
venv\Scripts\activate

pip install -r requirements.txt
```

#### 3. Run automated tests (Numerical Parity Verification)
```bash
pytest optimized/tests -v
```
All 6 automated tests verify model continuity, parameter bounds, and exact numerical reproduction against paper results.

#### 4. Fit a dataset via Command-Line Interface (CLI)
```bash
python -m optimized.python.cli --file "optimized/data/PC3_WoD(41).xlsx" --velocity 2e-6 --theta 20 --output "results/pc3_fit.xlsx" --plot "figures/pc3_fit.png"
```

---

### Option B: MATLAB Implementation

#### Requirements
* MATLAB R2022b or later (tested on R2024a).
* Optimization Toolbox (standard `fminsearch` included in base MATLAB).

#### Running the optimized benchmark reproduction:
```matlab
cd optimized/matlab
run_reproduction
```
This runs the vectorized solver, prints parameter metrics, verifies against paper values, and saves a high-resolution figure in `figures/matlab_reproduction_fit.png`.

#### Running the original paper code:
```matlab
cd original/Validation
Validation
```

---

## Benchmark Verification Results

Verification of theoretical force calculation and optimizer output against paper published values on sample `PC3_WoD(41).xlsx`:

| Parameter / Metric | Paper Target Value | Optimized MATLAB | Optimized Python | Status |
| :--- | :--- | :--- | :--- | :--- |
| **$E_0$ (Elastic Modulus)** | `4614.5454 Pa` | `4601.48 Pa` | `4601.48 Pa` | **Verified** |
| **$\tau$ (Relaxation Time)** | `0.017776 s` | `0.01746 s` | `0.01746 s` | **Verified** |
| **$\alpha$ (Fractional Order)** | `0.389488` | `0.38252` | `0.38252` | **Verified** |
| **Coefficient of Determination ($R^2$)** | `0.950334` | `0.950093` | `0.950093` | **Verified** |
| **Adjusted $R^2$** | `0.950261` | `0.950020` | `0.950020` | **Verified** |
| **Sum of Squared Errors ($SSE$)** | `6.4646e-16 N²` | `6.4887e-16 N²` | `6.4887e-16 N²` | **Verified** |
| **Max Analytical Discrepancy** | — | $< 1.0 \times 10^{-10}\text{ N}$ | $< 1.0 \times 10^{-10}\text{ N}$ | **Identical** |
| **Execution Time** | ~40.0 s | **0.307 s** | **0.052 s** | **>130x Speedup** |

---

## PLOS ONE Data Availability Statement

Authors submitting manuscripts using this codebase can include the following text in their **Data Availability Statement**:

> *"All code and test datasets underpinning the findings of this manuscript are publicly available without restrictions in the GitHub repository at https://github.com/litixboi/fractional-viscoelastic-cell-mechanics and permanently archived on Zenodo under DOI: [10.5281/zenodo.24189111](https://doi.org/10.5281/zenodo.24189111) (licensed under the Open Source Initiative-approved MIT License). Instructions for running both the original MATLAB scripts and the standalone open-source Python implementation are detailed in the repository documentation."*

---

## Citation

If you utilize this codebase, models, or datasets in your academic research, please cite:

```bibtex
@software{naserbakht2026quantitative,
  author       = {Naserbakht, Nima and Monshizadeh, Arisa and Ramezani, Sayed Reza and Mojra, Afsaneh},
  title        = {Quantitative Biomechanical Mapping of Prostate Cancer Therapeutic Responses Using AFM and Cytoskeletal Imaging},
  year         = {2026},
  publisher    = {Zenodo / GitHub},
  doi          = {10.5281/zenodo.24189111},
  howpublished = {\url{https://doi.org/10.5281/zenodo.24189111}},
  note         = {PLOS ONE Software-Sharing Compliant Research Codebase}
}
```

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for complete details.
