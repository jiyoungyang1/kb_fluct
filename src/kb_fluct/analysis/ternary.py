import numpy as np

class TernarySystem:
    """
    Represents a ternary solution system for Kirkwood-Buff Theory calculations.
    """
    def __init__(self, rho1, rho2, rho3, G11, G22, G33, G12, G13, G23, eta0=None):
        """
        Initialize the ternary system.
        
        Args:
            rho1 (float): Density of solvent.
            rho2 (float): Density of primary solute.
            rho3 (float): Density of co-solute.
            G11, G22, G33 (float): Self KBIs.
            G12, G13, G23 (float): Cross KBIs.
            eta0 (float, optional): Viscosity of pure solvent. Required for approximations.
        """
        self.rho1 = rho1
        self.rho2 = rho2
        self.rho3 = rho3
        self.G11 = G11
        self.G22 = G22
        self.G33 = G33
        self.G12 = G12
        self.G13 = G13
        self.G23 = G23
        self.eta0 = eta0

    def calculate_kD_exact(self):
        """
        Calculate exact diffusion interaction parameter (k_D) for a ternary system (Eq 37).
        """
        term = self.rho3 + self.rho3**2 * self.G33
        M11 = self.rho1 + self.rho1**2 * self.G11 - (self.rho1 * self.rho3 * self.G13)**2 / term
        M22 = self.rho2 + self.rho2**2 * self.G22 - (self.rho2 * self.rho3 * self.G23)**2 / term
        M12 = self.rho1 * self.rho2 * self.G12 - (self.rho1 * self.rho3 * self.G13) * (self.rho2 * self.rho3 * self.G23) / term
        
        num = M11 * (1 - M22) + M12**2
        den = M11 * M22 - M12**2
        return (1 / self.rho2) * (num / den)

    def calculate_kD_approx(self):
        """
        Calculate k_D for a ternary solution with low cosolute concentrations (Eq 43).
        Requires eta0 to be initialized.
        """
        if self.eta0 is None:
            raise ValueError("eta0 (solvent viscosity) must be provided for approximate calculations.")
            
        M_bar_11 = self.rho1 * self.eta0 - (self.rho1 * self.rho3 * self.G13)**2 / self.rho3
        M_bar_22 = self.rho2 * (1 + self.rho2 * self.G22) - (self.rho2 * self.rho3 * self.G23)**2 / self.rho3
        M_bar_12 = self.rho1 * self.rho2 * self.G12 - (self.rho1 * self.rho3 * self.G13) * (self.rho2 * self.rho3 * self.G23) / self.rho3
        
        num = M_bar_11 * (1 - M_bar_22) + M_bar_12**2
        den = M_bar_11 * M_bar_22 - M_bar_12**2
        return (1 / self.rho2) * (num / den)
