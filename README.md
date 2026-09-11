# kb_fluct: Kirkwood-Buff Density Fluctuation Analysis

`kb_fluct` is an object-oriented Python package designed for the thermodynamic analysis of Molecular Dynamics (MD) trajectories, specifically focusing on density fluctuations derived from Kirkwood-Buff Integrals (KBI).

It extracts structural geometry (such as tilt and eccentricity) from the density fluctuation distribution, allowing researchers to seamlessly detect force field artifacts and clustering behaviors.

## Installation

```bash
git clone https://github.com/jiyoung-yang/kb_fluct.git
cd kb_fluct
pip install .
```

## Quick Start

```python
from kb_fluct import FluctuationMatrix, EllipseAnalyzer

# Build S-matrix from KBIs
fluct = FluctuationMatrix(rho_1=33.3, rho_2=1.5).fit(G11=-0.02, G22=0.5, G12=-0.1)

# Analyze geometric properties of the density fluctuation ellipse
analyzer = EllipseAnalyzer().fit(fluct.S_matrix_)

print(f"Tilt Angle: {analyzer.tilt_angle_degrees_:.1f}°")
print(f"Eccentricity: {analyzer.eccentricity_:.2f}")
```

## Architecture
This package is divided into modular components:
- `kb_fluct.io`: Parsing GROMACS `.xvg` outputs.
- `kb_fluct.core`: Thermodynamic Simpson integration and `GangulyCorrection` for finite-size (closed NpT ensemble) effects.
- `kb_fluct.analysis`: Building the Fluctuation Matrix (`S-matrix`) and extracting its physical geometry.

## Example Notebooks

Comprehensive tutorials are available in the `examples/` directory:
- [a_kD_KBI_derivation.ipynb](examples/a_kD_KBI_derivation.ipynb): **a.** Mathematical Derivation of $k_D$ and 3D Schur Complement
- [b_density_fluctuation_NaCl_solution.ipynb](examples/b_density_fluctuation_NaCl_solution.ipynb): **b.** 2D Density Fluctuation Analysis in NaCl Solution
- [c_density_fluctuation_fg_cg_compare.ipynb](examples/c_density_fluctuation_fg_cg_compare.ipynb): **c.** Fine-Grained vs. Coarse-Grained Mixture Model Comparison



## References & Citations

If you use this package in your research, please ensure you cite the underlying theoretical frameworks:

1. **Kirkwood-Buff Theory (KBI):**
   > Kirkwood, J. G.; Buff, F. P. The Statistical Mechanical Theory of Solutions. I. *J. Chem. Phys.* **1951**, 19 (6), 774–777. [doi:10.1063/1.1748352](https://doi.org/10.1063/1.1748352)
   
2. **Ganguly Finite-Size Correction (Closed NpT Ensemble):**
   > Ganguly, P.; van der Vegt, N. F. A. Convergence of Sampling Kirkwood–Buff Integrals of Aqueous Solutions with Molecular Dynamics Simulations. *J. Chem. Theory Comput.* **2013**, 9 (3), 1347–1355. [doi:10.1021/ct301017q](https://doi.org/10.1021/ct301017q)

3. **Diffusion Interaction Parameter via Kirkwood-Buff Theory:**
   > Yang, J.; Brosz, M.; Adebar, N.; Burkert, O.; Schulze, M.; Smiatek, J. Influence of Co-Solutes and Solvents on Diffusion Interaction Parameters in Multicomponent Solutions: New Insights through the Kirkwood–Buff Theory. *J. Phys. Chem. B* **2025**, 129, 10381–10391. [doi:10.1021/acs.jpcb.5c03102](https://doi.org/10.1021/acs.jpcb.5c03102)
