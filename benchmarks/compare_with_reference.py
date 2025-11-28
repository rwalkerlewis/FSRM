#!/usr/bin/env python3
"""
SCEC Benchmark Verification Tool

Compares FSRM results with reference solutions from:
- SCEC website
- SeisSol
- Other verified codes

Usage:
    python compare_with_reference.py <benchmark_name>
    python compare_with_reference.py tpv5
    python compare_with_reference.py --all
"""

import sys
import os
import numpy as np
import h5py
import argparse
from pathlib import Path
import json

# Verification tolerances (stricter = better)
TOLERANCES = {
    'rupture_time': 0.05,      # 5% error allowed
    'slip': 0.05,              # 5% error
    'slip_rate': 0.10,         # 10% error (more variable)
    'rupture_velocity': 0.05,  # 5% error
    'energy': 0.01,            # 1% error (should be very accurate)
}

# Benchmark-specific tolerances (some are more challenging)
BENCHMARK_TOLERANCES = {
    'tpv13': {'slip_rate': 0.15},  # Plasticity adds variability
    'tpv34': {'slip_rate': 0.15, 'rupture_time': 0.10},  # TP is sensitive
    'tpv104': {'slip_rate': 0.15, 'rupture_time': 0.10}, # Strong VW is stiff
}

class BenchmarkVerifier:
    def __init__(self, benchmark_name, reference_dir='benchmarks/reference'):
        self.benchmark_name = benchmark_name
        self.reference_dir = Path(reference_dir) / benchmark_name
        self.output_dir = Path('output') / benchmark_name
        
        # Results storage
        self.metrics = {}
        self.passed = []
        self.failed = []
        
    def load_fsrm_data(self):
        """Load FSRM output data"""
        print(f"Loading FSRM data for {self.benchmark_name}...")
        
        # Try HDF5 first
        fault_file = self.output_dir / 'fault_output.h5'
        if fault_file.exists():
            with h5py.File(fault_file, 'r') as f:
                self.fsrm_data = {
                    'rupture_time': np.array(f['rupture_time']),
                    'slip': np.array(f['slip']),
                    'slip_rate': np.array(f['slip_rate']),
                    'shear_stress': np.array(f['shear_stress']),
                    'normal_stress': np.array(f['normal_stress']),
                }
                if 'state_variable' in f:
                    self.fsrm_data['state_variable'] = np.array(f['state_variable'])
                if 'pore_pressure' in f:
                    self.fsrm_data['pore_pressure'] = np.array(f['pore_pressure'])
            return True
        
        print(f"  ✗ FSRM output not found at {fault_file}")
        return False
    
    def load_reference_data(self):
        """Load reference solution"""
        print(f"Loading reference data...")
        
        ref_file = self.reference_dir / 'reference_solution.h5'
        if not ref_file.exists():
            print(f"  ⚠ Reference solution not found at {ref_file}")
            return False
        
        with h5py.File(ref_file, 'r') as f:
            self.reference_data = {
                'rupture_time': np.array(f['rupture_time']),
                'slip': np.array(f['slip']),
                'slip_rate': np.array(f['slip_rate']),
            }
            
            # Optional fields
            for field in ['shear_stress', 'normal_stress', 'state_variable', 'pore_pressure']:
                if field in f:
                    self.reference_data[field] = np.array(f[field])
        
        return True
    
    def compute_error_metric(self, field, metric_type='l2_relative'):
        """Compute error between FSRM and reference"""
        
        fsrm = self.fsrm_data[field]
        ref = self.reference_data[field]
        
        # Ensure same shape (interpolate if needed)
        if fsrm.shape != ref.shape:
            print(f"  ⚠ Shape mismatch for {field}: {fsrm.shape} vs {ref.shape}")
            # TODO: Implement interpolation
            return None
        
        if metric_type == 'l2_relative':
            # L2 relative error
            diff = fsrm - ref
            error = np.sqrt(np.mean(diff**2)) / np.sqrt(np.mean(ref**2))
        
        elif metric_type == 'max_relative':
            # Maximum relative error
            error = np.max(np.abs(fsrm - ref) / np.abs(ref))
        
        elif metric_type == 'rmse':
            # Root mean square error
            error = np.sqrt(np.mean((fsrm - ref)**2))
        
        else:
            raise ValueError(f"Unknown metric type: {metric_type}")
        
        return error
    
    def verify_rupture_time(self):
        """Verify rupture arrival time"""
        print("\n  Checking rupture time...")
        
        error = self.compute_error_metric('rupture_time', 'l2_relative')
        if error is None:
            return False
        
        # Get tolerance
        tol = BENCHMARK_TOLERANCES.get(self.benchmark_name, {}).get('rupture_time',
                                                                     TOLERANCES['rupture_time'])
        
        self.metrics['rupture_time_error'] = error
        
        if error < tol:
            print(f"    ✓ Rupture time error: {error:.4f} < {tol:.4f}")
            self.passed.append('rupture_time')
            return True
        else:
            print(f"    ✗ Rupture time error: {error:.4f} >= {tol:.4f}")
            self.failed.append('rupture_time')
            return False
    
    def verify_slip_distribution(self):
        """Verify final slip distribution"""
        print("\n  Checking slip distribution...")
        
        error = self.compute_error_metric('slip', 'l2_relative')
        if error is None:
            return False
        
        tol = BENCHMARK_TOLERANCES.get(self.benchmark_name, {}).get('slip',
                                                                     TOLERANCES['slip'])
        
        self.metrics['slip_error'] = error
        
        if error < tol:
            print(f"    ✓ Slip error: {error:.4f} < {tol:.4f}")
            self.passed.append('slip')
            return True
        else:
            print(f"    ✗ Slip error: {error:.4f} >= {tol:.4f}")
            self.failed.append('slip')
            return False
    
    def verify_slip_rate(self):
        """Verify peak slip rate"""
        print("\n  Checking slip rate...")
        
        # Compare peak slip rates
        fsrm_peak = np.max(self.fsrm_data['slip_rate'])
        ref_peak = np.max(self.reference_data['slip_rate'])
        
        error = np.abs(fsrm_peak - ref_peak) / ref_peak
        
        tol = BENCHMARK_TOLERANCES.get(self.benchmark_name, {}).get('slip_rate',
                                                                     TOLERANCES['slip_rate'])
        
        self.metrics['slip_rate_peak_error'] = error
        self.metrics['fsrm_peak_slip_rate'] = fsrm_peak
        self.metrics['ref_peak_slip_rate'] = ref_peak
        
        if error < tol:
            print(f"    ✓ Peak slip rate: {fsrm_peak:.3f} m/s vs {ref_peak:.3f} m/s")
            print(f"      Error: {error:.4f} < {tol:.4f}")
            self.passed.append('slip_rate')
            return True
        else:
            print(f"    ✗ Peak slip rate: {fsrm_peak:.3f} m/s vs {ref_peak:.3f} m/s")
            print(f"      Error: {error:.4f} >= {tol:.4f}")
            self.failed.append('slip_rate')
            return False
    
    def verify_rupture_velocity(self):
        """Estimate and verify average rupture velocity"""
        print("\n  Checking rupture velocity...")
        
        # Compute rupture velocity from rupture time
        rt = self.fsrm_data['rupture_time']
        
        # Simple estimate: distance / time
        # (This is approximate - proper method uses contours)
        
        # For now, use simple metric
        # TODO: Implement proper rupture velocity computation
        
        print("    ⚠ Rupture velocity check not yet implemented")
        return True
    
    def verify_energy_balance(self):
        """Check energy conservation"""
        print("\n  Checking energy balance...")
        
        # Load energy data if available
        energy_file = self.output_dir / 'energy.txt'
        if not energy_file.exists():
            print("    ⚠ Energy output not found")
            return True
        
        # Read energy data
        data = np.loadtxt(energy_file)
        time = data[:, 0]
        kinetic = data[:, 1]
        elastic = data[:, 2]
        fracture = data[:, 3]
        
        total = kinetic + elastic + fracture
        
        # Check conservation (should be constant or monotonic)
        energy_change = np.abs(total[-1] - total[0]) / total[0]
        
        tol = TOLERANCES['energy']
        
        self.metrics['energy_balance_error'] = energy_change
        
        if energy_change < tol:
            print(f"    ✓ Energy balance: {energy_change:.6f} < {tol:.6f}")
            self.passed.append('energy')
            return True
        else:
            print(f"    ✗ Energy balance: {energy_change:.6f} >= {tol:.6f}")
            self.failed.append('energy')
            return False
    
    def run_verification(self):
        """Run all verification checks"""
        print(f"\n{'='*60}")
        print(f"Verifying SCEC Benchmark: {self.benchmark_name}")
        print(f"{'='*60}")
        
        # Load data
        if not self.load_fsrm_data():
            print("\n✗ Cannot load FSRM data")
            return False
        
        if not self.load_reference_data():
            print("\n⚠ Cannot load reference data (will skip comparison)")
            # Can still check internal consistency
        
        # Run checks
        checks = []
        
        if 'rupture_time' in self.reference_data:
            checks.append(self.verify_rupture_time())
        
        if 'slip' in self.reference_data:
            checks.append(self.verify_slip_distribution())
        
        if 'slip_rate' in self.reference_data:
            checks.append(self.verify_slip_rate())
        
        checks.append(self.verify_rupture_velocity())
        checks.append(self.verify_energy_balance())
        
        # Summary
        print(f"\n{'='*60}")
        print("Verification Summary")
        print(f"{'='*60}")
        print(f"Passed: {len(self.passed)}/{len(self.passed) + len(self.failed)}")
        
        if self.passed:
            print("\nPassed checks:")
            for check in self.passed:
                print(f"  ✓ {check}")
        
        if self.failed:
            print("\nFailed checks:")
            for check in self.failed:
                print(f"  ✗ {check}")
        
        # Save metrics
        metrics_file = self.output_dir / 'verification_metrics.json'
        with open(metrics_file, 'w') as f:
            json.dump(self.metrics, f, indent=2)
        print(f"\nMetrics saved to {metrics_file}")
        
        # Overall result
        if len(self.failed) == 0:
            print(f"\n✓ {self.benchmark_name} PASSED verification")
            return True
        else:
            print(f"\n✗ {self.benchmark_name} FAILED verification")
            return False

def main():
    parser = argparse.ArgumentParser(description='Verify SCEC benchmark results')
    parser.add_argument('benchmark', nargs='?', help='Benchmark name (e.g., tpv5)')
    parser.add_argument('--all', action='store_true', help='Verify all benchmarks')
    parser.add_argument('--reference-dir', default='benchmarks/reference',
                       help='Directory containing reference solutions')
    
    args = parser.parse_args()
    
    if args.all:
        # Verify all available benchmarks
        benchmarks = ['tpv5', 'tpv10', 'tpv13', 'tpv16', 'tpv34', 'tpv101', 'tpv104']
        
        results = {}
        for bm in benchmarks:
            verifier = BenchmarkVerifier(bm, args.reference_dir)
            results[bm] = verifier.run_verification()
        
        # Overall summary
        print(f"\n\n{'='*60}")
        print("Overall Summary")
        print(f"{'='*60}")
        for bm, passed in results.items():
            status = "✓ PASS" if passed else "✗ FAIL"
            print(f"{bm:10s}: {status}")
        
        # Exit code
        if all(results.values()):
            sys.exit(0)
        else:
            sys.exit(1)
    
    elif args.benchmark:
        verifier = BenchmarkVerifier(args.benchmark, args.reference_dir)
        success = verifier.run_verification()
        sys.exit(0 if success else 1)
    
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == '__main__':
    main()
