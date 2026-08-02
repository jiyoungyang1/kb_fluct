import numpy as np

class BinarySystem:
    """
    Represents a binary solution system for Kirkwood-Buff Theory calculations.
    """
    def __init__(self, rho1, rho2, G11, G22, G12, eta0=None):
        """
        Initialize the binary system.
        
        Args:
            rho1 (float): Density of solvent.
            rho2 (float): Density of solute.
            G11 (float): KBI between solvent molecules.
            G22 (float): KBI between solute molecules.
            G12 (float): Cross KBI between solvent and solute.
            eta0 (float, optional): Viscosity of the pure solvent. Required for approximations.
        """
        self.rho1 = rho1
        self.rho2 = rho2
        self.G11 = G11
        self.G22 = G22
        self.G12 = G12
        self.eta0 = eta0

    def calculate_Gamma(self):
        """
        Calculate the exact thermodynamic factor (Gamma).
        """
        detB = (self.rho1 + self.rho1**2 * self.G11) * (self.rho2 + self.rho2**2 * self.G22) - (self.rho1**2 * self.rho2**2 * self.G12**2)
        A22 = (self.rho1 * (1 + self.rho1 * self.G11)) / detB
        return self.rho2 * A22

    def calculate_kD_exact(self):
        """
        Calculate exact diffusion interaction parameter (k_D).
        """
        Gamma = self.calculate_Gamma()
        return (Gamma - 1) / self.rho2

    def calculate_kD_approx(self):
        """
        Calculate approximate k_D (Eq 23).
        Requires eta0 to be initialized.
        """
        if self.eta0 is None:
            raise ValueError("eta0 (solvent viscosity) must be provided for approximate calculations.")
        num = (self.rho1 * self.G12**2 / self.eta0) - self.G22
        den = 1 + self.rho2 * (self.G22 - (self.rho1 * self.G12**2 / self.eta0))
        return num / den

    def get_B_matrix(self):
        """
        Calculate the cofactor matrix B for density fluctuations.
        """
        B11 = self.rho1 + self.rho1**2 * self.G11
        B22 = self.rho2 + self.rho2**2 * self.G22
        B12 = self.rho1 * self.rho2 * self.G12
        B21 = B12
        return np.array([[B11, B12], [B21, B22]])
