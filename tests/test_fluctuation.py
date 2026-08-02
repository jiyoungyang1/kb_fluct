import pytest
import numpy as np
from kb_fluct import FluctuationMatrix, EllipseAnalyzer

def test_fluctuation_matrix_and_ellipse():
    # Mock KBI values
    rho_1, rho_2 = 33.3, 1.5
    G11, G22, G12 = -0.02, 0.5, -0.1

    # Fit FluctuationMatrix
    fluct = FluctuationMatrix(rho_1, rho_2).fit(G11, G22, G12)
    
    assert fluct.S_matrix_.shape == (2, 2)
    assert np.isclose(fluct.S11_, rho_1 + (rho_1**2)*G11)
    assert np.isclose(fluct.S12_, rho_1*rho_2*G12)

    # Analyze Ellipse
    analyzer = EllipseAnalyzer().fit(fluct.S_matrix_)
    
    assert analyzer.is_positive_definite_ is True
    assert len(analyzer.eigenvalues_) == 2
    assert analyzer.tilt_angle_degrees_ is not None
    assert analyzer.eccentricity_ >= 0.0 and analyzer.eccentricity_ <= 1.0
