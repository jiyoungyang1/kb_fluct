import numpy as np
try:
    from scipy.integrate import simpson
except ImportError:
    from scipy.integrate import simps as simpson

def calc_kbi(r: np.ndarray, g: np.ndarray) -> float:
    """
    Calculates the Kirkwood-Buff Integral (G_ij) from a radial distribution function.
    Formula: G_ij = integral_0^infty 4*pi*r^2 * (g(r) - 1) dr
    """
    integrand = 4.0 * np.pi * (r**2) * (g - 1.0)
    return float(simpson(integrand, r))
