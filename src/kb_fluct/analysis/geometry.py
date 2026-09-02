import numpy as np
import scipy.stats as stats

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
        self.results = {}

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
        self.results['is_positive_definite'] = self.is_positive_definite_

        # Eigen decomposition
        eigvals, eigvecs = np.linalg.eigh(S_matrix)
        
        # Sort eigenvalues and eigenvectors in descending order
        idx = eigvals.argsort()[::-1]
        self.eigenvalues_ = eigvals[idx]
        self.eigenvectors_ = eigvecs[:, idx]
        self.results['eigenvalues'] = self.eigenvalues_
        self.results['eigenvectors'] = self.eigenvectors_
        
        l1, l2 = self.eigenvalues_ # l1 >= l2
        
        # Tilt angle of the major axis (v1)
        v1 = self.eigenvectors_[:, 0]
        angle_rad = np.arctan2(v1[1], v1[0])
        self.tilt_angle_degrees_ = np.degrees(angle_rad)
        self.results['tilt_angle_degrees'] = self.tilt_angle_degrees_
        
        # Eccentricity e = sqrt(1 - (l2/l1)^2) if we consider axes length proportional to sqrt(eigenval)
        # Note: Major axis a ~ sqrt(l1), minor axis b ~ sqrt(l2)
        if l1 <= 0 or l2 < 0:
            self.eccentricity_ = float("nan") # Unphysical
        else:
            self.eccentricity_ = np.sqrt(1.0 - (l2/l1))
        self.results['eccentricity'] = self.eccentricity_
            
        # Area = pi * a * b
        if l1 >= 0 and l2 >= 0:
            self.area_ = np.pi * np.sqrt(l1 * l2)
        else:
            self.area_ = float("nan")
        self.results['area'] = self.area_
            
        return self

class DistributionAnalyzer:
    """
    Extracts statistical metrics (skewness, kurtosis) from a set of 2D data points.
    Follows scikit-learn style API.
    """
    def __init__(self):
        self.results = {}

    def fit(self, points: np.ndarray):
        """
        Computes skewness and kurtosis for a set of 2D points (shape: N x 2).
        points[:, 0] corresponds to the first dimension (e.g., solvent),
        points[:, 1] corresponds to the second dimension (e.g., solute).
        """
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError("points must be a 2D numpy array of shape (N, 2).")

        # Calculate skewness and Fisher (excess) kurtosis for each dimension
        # skewness = 0 and kurtosis = 0 indicates a perfect normal distribution.
        self.results['skewness'] = stats.skew(points, axis=0)
        self.results['kurtosis'] = stats.kurtosis(points, axis=0, fisher=True)
        
        return self
