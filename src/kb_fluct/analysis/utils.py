import pandas as pd
import numpy as np

def summarize_results(*analyzers, names=None):
    """
    Takes one or more analyzer objects (e.g., FluctuationMatrix, EllipseAnalyzer)
    and formats their `results` dictionaries into a neat pandas DataFrame for display.
    Benchmarked against MDAnalysis/RDKit patterns of using pandas for tabular result viewing.
    """
    summaries = []
    
    for i, analyzer in enumerate(analyzers):
        if not hasattr(analyzer, 'results') or not analyzer.results:
            continue
            
        flat_results = {}
        for key, value in analyzer.results.items():
            # Skip 2D matrices for a clean summary table, unless they are small
            if isinstance(value, np.ndarray):
                if value.ndim == 1:
                    # Flatten 1D arrays like skewness and kurtosis
                    for j, val in enumerate(value):
                        comp_name = f"comp_{j+1}"
                        flat_results[f"{key}_{comp_name}"] = round(val, 4)
                elif value.ndim == 2:
                    # Skip full 2D matrices to keep the table clean
                    flat_results[key] = "<matrix>"
            elif isinstance(value, float):
                flat_results[key] = round(value, 4)
            else:
                flat_results[key] = value
                
        name = names[i] if names and i < len(names) else analyzer.__class__.__name__
        flat_results['Analyzer'] = name
        summaries.append(flat_results)
        
    if not summaries:
        return pd.DataFrame()
        
    df = pd.DataFrame(summaries)
    # Move 'Analyzer' column to the front
    cols = ['Analyzer'] + [c for c in df.columns if c != 'Analyzer']
    return df[cols].fillna('-')
