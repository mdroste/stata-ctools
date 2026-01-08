#!/usr/bin/env python3
"""
Benchmark quantile regression in Python vs Stata cqreg/qreg
Uses same data generation as Stata benchmark: N=5M, y = 2*x + noise
"""

import numpy as np
import time

# Try to import statsmodels for quantile regression
try:
    import statsmodels.api as sm
    from statsmodels.regression.quantile_regression import QuantReg
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("statsmodels not installed. Install with: pip install statsmodels")

def generate_data(n=5_000_000, seed=12345):
    """Generate same data as Stata benchmark"""
    np.random.seed(seed)
    x = np.random.randn(n)
    y = 2 * x + np.random.randn(n) * 5
    return y, x

def run_statsmodels_qreg(y, x, quantile=0.5):
    """Run quantile regression using statsmodels"""
    # Add constant
    X = sm.add_constant(x)
    
    # Fit model
    model = QuantReg(y, X)
    result = model.fit(q=quantile)
    
    return result

def main():
    print("=" * 60)
    print("Python Quantile Regression Benchmark")
    print("N = 5,000,000, K = 2 (constant + x)")
    print("=" * 60)
    print()
    
    # Generate data
    print("Generating data...")
    t0 = time.perf_counter()
    y, x = generate_data()
    t_gen = time.perf_counter() - t0
    print(f"Data generation: {t_gen:.3f} seconds")
    print()
    
    if HAS_STATSMODELS:
        print("Running statsmodels QuantReg (interior point method)...")
        t0 = time.perf_counter()
        result = run_statsmodels_qreg(y, x, quantile=0.5)
        t_qreg = time.perf_counter() - t0
        
        print(f"\nstatsmodels QuantReg time: {t_qreg:.3f} seconds")
        print()
        print("Results:")
        print("-" * 40)
        print(f"{'Variable':<12} {'Coef':>12} {'Std.Err':>12} {'t':>10} {'P>|t|':>10}")
        print("-" * 40)
        
        # Get results
        params = result.params
        bse = result.bse
        tvalues = result.tvalues
        pvalues = result.pvalues
        
        var_names = ['const', 'x']
        for i, name in enumerate(var_names):
            print(f"{name:<12} {params[i]:>12.6f} {bse[i]:>12.6f} {tvalues[i]:>10.2f} {pvalues[i]:>10.4f}")
        
        print("-" * 40)
        print()
        print("Comparison with Stata (expected: const ≈ 0, x ≈ 2):")
        print(f"  Python const: {params[0]:.6f}")
        print(f"  Python x:     {params[1]:.6f}")
    
    print()
    print("=" * 60)
    print("Summary:")
    print("=" * 60)
    if HAS_STATSMODELS:
        print(f"  statsmodels QuantReg: {t_qreg:.3f} seconds")
    print()
    print("For comparison, Stata times on same data:")
    print("  qreg:             ~10 seconds")
    print("  cqreg (residual): ~2.2 seconds")
    print("  cqreg (fitted):   ~4.6 seconds")

if __name__ == "__main__":
    main()
