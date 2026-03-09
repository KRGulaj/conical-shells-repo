# conical-shells-repo
Aerodynamic Characterization of Hollow Conical Shells in Free Fall - A Comparative Experimental and Numerical Study on the Influence of Apex Angle on Subsonic Drag Coefficients
# Abstract
This repository contains the experimental data, computational models, and source code for an independent research project investigating the aerodynamic behavior of hollow conical shells. The primary objective of this study is to bridge theoretical fluid dynamics with empirical validation by systematically correlating free-fall kinematics with computational fluid dynamics (CFD) predictions.

# Methodology
The research methodology is structured around a dual-approach validation process to accurately determine drag coefficients across varying bluff-body geometries.

Initially, a physical experimental framework was engineered to capture the free-fall kinematics of the conical shells using high-framerate video tracking. To guarantee the mathematical rigor of the empirical data, a comprehensive measurement uncertainty analysis was executed utilizing the Root Sum Square (RSS) framework, tightly quantifying all experimental error margins.

In parallel with the physical testing, the computational phase involved steady-state RANS simulations executed in OpenFOAM. Utilizing the k-omega SST turbulence model, this phase focused on extracting complex wake topologies and computationally deriving the drag coefficients for the corresponding geometries.

To synthesize the empirical and numerical data streams, custom Python routines were developed. These scripts automate the data processing and apply Weighted Least Squares (WLS) regression, providing a critical, quantitative synthesis of the physical test results against the CFD predictions.

## Computational & Analytical Stack
* **CFD Solvers & Meshing:** OpenFOAM
* **Data Processing & Regression:** Python (NumPy, SciPy, Matplotlib)
* **Computer-Aided Design (CAD):** CATIA
* **Typesetting & Documentation:** LaTeX

# Key Results
The integration of the experimental kinematics with the numerical CFD models yielded a strong correlation in predicting both terminal velocities and aerodynamic drag coefficients as a function of the conical apex angle. A central output of this comparative analysis is visualized below:

![Correlation between Experimental and Computational Drag Coefficients](results/figure_6_7.png)
*(Figure 6.7: Correlation between empirically derived and computationally predicted drag coefficients across tested geometries)*

> **Note:** For a comprehensive breakdown of the mathematical derivations, boundary conditions, and the exact physical experimental setup, please refer to the complete technical report located in the `docs/` directory.

# Repository Structure
* `/docs` - The final technical report (PDF) detailing the theoretical background and findings.
* `/src` - Python source code for data analysis, plotting, and WLS regression modeling.
* `/data` - Experimental kinematics datasets (categorized into raw and processed) along with uncertainty calculation matrices.
* `/cfd` - Core configuration files, mesh setups, and solver logs for the OpenFOAM/SimScale simulations.
* `/results` - Post-processing visualizations, including wake topologies, pressure gradients, and synthesized data plots.

# License
This project is open-source and available under the MIT License.
