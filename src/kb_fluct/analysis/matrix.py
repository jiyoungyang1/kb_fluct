import numpy as np

class FluctuationMatrix:
    """
    Constructs the density fluctuation covariance matrix (S matrix) from Kirkwood-Buff Integrals.
    
    Based on the theoretical framework defined in:
    Yang, J.; Brosz, M.; Adebar, N.; Burkert, O.; Schulze, M.; Smiatek, J. 
    Influence of Co-Solutes and Solvents on Diffusion Interaction Parameters 
    in Multicomponent Solutions: New Insights through the Kirkwood-Buff Theory. 
    J. Phys. Chem. B 2025, 129, 10381-10391.
    """
    def __init__(self, rho_1: float, rho_2: float):
        """
        Args:
            rho_1: Number density of species 1 (solvent).
            rho_2: Number density of species 2 (solute).
        """
        self.rho_1 = float(rho_1)
        self.rho_2 = float(rho_2)
        self.S_matrix_ = None
        self.S11_ = None
        self.S22_ = None
        self.S12_ = None
        self.results = {}

    def fit(self, G11: float, G22: float, G12: float):
        """
        Calculates the elements of the S matrix.
        Args:
            G11: KBI for solvent-solvent.
            G22: KBI for solute-solute.
            G12: KBI for solvent-solute.
        """
        self.S11_ = self.rho_1 + (self.rho_1**2) * G11
        self.S22_ = self.rho_2 + (self.rho_2**2) * G22
        self.S12_ = self.rho_1 * self.rho_2 * G12
        
        self.S_matrix_ = np.array([
            [self.S11_, self.S12_],
            [self.S12_, self.S22_]
        ])
        
        self.results['S11'] = self.S11_
        self.results['S22'] = self.S22_
        self.results['S12'] = self.S12_
        self.results['S_matrix'] = self.S_matrix_
        
        return self
