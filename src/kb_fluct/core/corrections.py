import numpy as np

class GangulyCorrection:
    """
    Applies the Ganguly finite-size correction to an RDF computed from a closed (NpT) ensemble.
    
    Based on the methodology described in:
    Ganguly, P.; van der Vegt, N. F. A. Convergence of Sampling Kirkwood-Buff Integrals 
    of Aqueous Solutions with Molecular Dynamics Simulations. 
    J. Chem. Theory Comput. 2013, 9 (3), 1347-1355.
    """
    def __init__(self, N: float, V: float, is_self_pair: bool):
        """
        Args:
            N: Number of particles of the TARGET species (the "sel" species in gmx rdf).
            V: Average volume of the simulation box.
            is_self_pair: True if calculating interactions between the same species (e.g., water-water, ion-ion).
        """
        self.N = float(N)
        self.V = float(V)
        self.delta = 1 if is_self_pair else 0
        self.rho_j = self.N / self.V

    def transform(self, r: np.ndarray, g: np.ndarray) -> np.ndarray:
        """
        Applies the correction to the raw g(r).
        Returns:
            Corrected g(r) array.
        """
        dr = r[1] - r[0]
        fr = 1.0 - (4.0/3.0) * (np.pi/self.V) * (r**3)
        
        f_r = 4.0 * np.pi * (g - 1.0) * (r**2)
        rkbi = np.cumsum(f_r * dr)
        
        delta_N = self.rho_j * rkbi
        corr_numer = self.N * fr
        corr_denom = corr_numer - delta_N - self.delta
        
        with np.errstate(divide="ignore", invalid="ignore"):
            corr = np.where(corr_denom != 0, corr_numer / corr_denom, 1.0)
            
        return g * corr
