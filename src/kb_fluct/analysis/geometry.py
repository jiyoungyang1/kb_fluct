import numpy as np

class EllipseAnalyzer:
    """
    Extracts geometric and mathematical properties (eigenvalues, eigenvectors, tilt, eccentricity, area) 
    from a 2x2 Covariance Matrix (S matrix).
    Follows scikit-learn style API.
    """
    def __init__(self):
        self.eigenvalues_ = None
        self.eigenvectors_ = None
        self.tilt_angle_degrees_ = None
        self.eccentricity_ = None
        self.area_ = None
        self.is_positive_definite_ = None

    def fit(self, S_matrix: np.ndarray):
        """
        Computes properties from the S matrix.
        """
        if S_matrix.shape != (2, 2):
            raise ValueError("S_matrix must be a 2x2 numpy array.")
            
        # Check if matrix is positive semi-definite (physical requirement for covariance)
        try:
            np.linalg.cholesky(S_matrix)
            self.is_positive_definite_ = True
        except np.linalg.LinAlgError:
            self.is_positive_definite_ = False

        # Eigen decomposition
        eigvals, eigvecs = np.linalg.eigh(S_matrix)
        
        # Sort eigenvalues and eigenvectors in descending order
        idx = eigvals.argsort()[::-1]
        self.eigenvalues_ = eigvals[idx]
        self.eigenvectors_ = eigvecs[:, idx]
        
        l1, l2 = self.eigenvalues_ # l1 >= l2
        
        # Tilt angle of the major axis (v1)
        v1 = self.eigenvectors_[:, 0]
        angle_rad = np.arctan2(v1[1], v1[0])
        self.tilt_angle_degrees_ = np.degrees(angle_rad)
        
        # Eccentricity e = sqrt(1 - (l2/l1)^2) if we consider axes length proportional to sqrt(eigenval)
        # Note: Major axis a ~ sqrt(l1), minor axis b ~ sqrt(l2)
        if l1 <= 0 or l2 < 0:
            self.eccentricity_ = float("nan") # Unphysical
        else:
            self.eccentricity_ = np.sqrt(1.0 - (l2/l1))
            
        # Area = pi * a * b
        if l1 >= 0 and l2 >= 0:
            self.area_ = np.pi * np.sqrt(l1 * l2)
        else:
            self.area_ = float("nan")
            
        return self
