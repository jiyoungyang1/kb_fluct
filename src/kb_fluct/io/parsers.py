import numpy as np

def read_xvg(filepath: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Reads a GROMACS .xvg file and returns the first two columns (e.g., r and g(r)).
    Ignores comments starting with @ or #.
    """
    x, y = [], []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("@") or line.startswith("#"): 
                continue
            parts = line.split()
            if len(parts) >= 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
    return np.array(x), np.array(y)
