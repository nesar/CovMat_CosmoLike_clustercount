#!/usr/bin/env python3
"""
Validation and Testing Script for CosmoLike-Compatible Implementations

This script tests the existing implementations against CosmoLike standards
and the paper equations.
"""

import numpy as np
import sys
from pathlib import Path

# Import the implementations
try:
    from cosmolike_compatible_cluster_covariance_Y1 import CosmoLikeClusterCovariance
    from cov_mat_cosmolike import ClusterNumberCountCovariance
except ImportError as e:
    print(f"Error importing modules: {e}")
    print("Make sure you're in the correct directory")
    sys.exit(1)

class CosmoLikeValidator:
    """Validate implementations against CosmoLike standards"""

    # CosmoLike official specifications
    COSMOLIKE_RICHNESS_BINS = [(20, 30), (30, 45), (45, 70), (70, 120), (120, 220)]

    COSMOLIKE_Y1 = {
        'survey_area': 12300.0,
        'n_z_bins': 3,
        'z_bins': [(0.2, 0.5), (0.5, 0.8), (0.8, 1.1)],
        'n_data_points': 15  # 3 z-bins × 5 richness bins
    }

    COSMOLIKE_Y10 = {
        'survey_area': 16500.0,
        'n_z_bins': 4,
        'z_bins': [(0.2, 0.45), (0.45, 0.7), (0.7, 0.95), (0.95, 1.2)],
        'n_data_points': 20  # 4 z-bins × 5 richness bins
    }

    def __init__(self):
        self.results = {}

    def test_binning_scheme(self, calculator):
        """Test if richness binning matches CosmoLike standard"""
        print("\n" + "="*70)
        print("TEST 1: Binning Scheme Validation")
        print("="*70)

        richness_bins = calculator.lambda_bins
        print(f"\nImplementation richness bins:")
        for i, (lmin, lmax) in enumerate(richness_bins):
            print(f"  Bin {i+1}: {lmin}-{lmax}")

        print(f"\nCosmoLike standard richness bins:")
        for i, (lmin, lmax) in enumerate(self.COSMOLIKE_RICHNESS_BINS):
            print(f"  Bin {i+1}: {lmin}-{lmax}")

        # Check if bins match
        bins_match = richness_bins == self.COSMOLIKE_RICHNESS_BINS

        if bins_match:
            print(f"\n✅ PASS: Richness bins EXACTLY match CosmoLike standard")
            self.results['binning'] = 'PASS'
        else:
            print(f"\n❌ FAIL: Richness bins DO NOT match CosmoLike standard")
            print(f"   Expected: {self.COSMOLIKE_RICHNESS_BINS}")
            print(f"   Got: {richness_bins}")
            self.results['binning'] = 'FAIL'

        return bins_match

    def test_survey_configuration(self, calculator, config_name):
        """Test if survey configuration matches CosmoLike Y1/Y10"""
        print("\n" + "="*70)
        print(f"TEST 2: Survey Configuration Validation ({config_name})")
        print("="*70)

        if config_name == "Y1":
            expected = self.COSMOLIKE_Y1
        elif config_name == "Y10":
            expected = self.COSMOLIKE_Y10
        else:
            print(f"Unknown configuration: {config_name}")
            return False

        # Check survey area
        area_match = abs(calculator.survey_area - expected['survey_area']) < 1.0
        print(f"\nSurvey area:")
        print(f"  Expected: {expected['survey_area']} deg²")
        print(f"  Got: {calculator.survey_area} deg²")
        print(f"  Match: {'✅ YES' if area_match else '❌ NO'}")

        # Check number of redshift bins
        z_bins_match = calculator.n_z_bins == expected['n_z_bins']
        print(f"\nNumber of redshift bins:")
        print(f"  Expected: {expected['n_z_bins']}")
        print(f"  Got: {calculator.n_z_bins}")
        print(f"  Match: {'✅ YES' if z_bins_match else '❌ NO'}")

        # Check total data points
        n_data = calculator.n_z_bins * calculator.n_lambda_bins
        data_match = n_data == expected['n_data_points']
        print(f"\nTotal cluster data points:")
        print(f"  Expected: {expected['n_data_points']}")
        print(f"  Got: {n_data}")
        print(f"  Match: {'✅ YES' if data_match else '❌ NO'}")

        all_match = area_match and z_bins_match and data_match

        if all_match:
            print(f"\n✅ PASS: Survey configuration matches CosmoLike {config_name}")
            self.results[f'survey_{config_name}'] = 'PASS'
        else:
            print(f"\n❌ FAIL: Survey configuration does NOT match CosmoLike {config_name}")
            self.results[f'survey_{config_name}'] = 'FAIL'

        return all_match

    def test_cluster_counts_realism(self, calculator):
        """Test if cluster counts are in realistic range"""
        print("\n" + "="*70)
        print("TEST 3: Cluster Count Realism")
        print("="*70)

        # Calculate a few representative bins
        print("\nCalculating sample cluster counts...")

        # Low-z, low-richness (should be highest)
        N_low_z_low_rich = calculator.cluster_number_count(0.2, 0.4, 20, 30)

        # Mid-z, mid-richness
        N_mid_z_mid_rich = calculator.cluster_number_count(0.5, 0.7, 45, 70)

        # High-z, high-richness (should be lowest)
        N_high_z_high_rich = calculator.cluster_number_count(0.8, 1.0, 120, 220)

        print(f"\nSample counts:")
        print(f"  Low-z, low-richness (z=0.2-0.4, λ=20-30): {N_low_z_low_rich:.1f}")
        print(f"  Mid-z, mid-richness (z=0.5-0.7, λ=45-70): {N_mid_z_mid_rich:.1f}")
        print(f"  High-z, high-richness (z=0.8-1.0, λ=120-220): {N_high_z_high_rich:.1f}")

        # Realistic ranges (order of magnitude checks)
        # For LSST-like survey: few hundred to few thousand per bin
        realistic_min = 10  # At least 10 clusters per bin
        realistic_max = 1e6  # Not more than a million per bin

        counts_realistic = (
            realistic_min < N_low_z_low_rich < realistic_max and
            realistic_min < N_mid_z_mid_rich < realistic_max and
            realistic_min < N_high_z_high_rich < realistic_max
        )

        # Check trend: counts should decrease with richness and redshift
        trend_correct = N_low_z_low_rich > N_mid_z_mid_rich > N_high_z_high_rich

        print(f"\nRealism check:")
        print(f"  Counts in realistic range [{realistic_min:.0f}, {realistic_max:.0e}]: {'✅ YES' if counts_realistic else '❌ NO'}")
        print(f"  Correct trend (decreasing with z and λ): {'✅ YES' if trend_correct else '❌ NO'}")

        if counts_realistic and trend_correct:
            print(f"\n✅ PASS: Cluster counts are realistic")
            self.results['counts'] = 'PASS'
        else:
            print(f"\n⚠️  WARNING: Cluster counts may not be realistic")
            self.results['counts'] = 'WARNING'

        return counts_realistic and trend_correct

    def test_bias_values(self, calculator):
        """Test if bias values are physically reasonable"""
        print("\n" + "="*70)
        print("TEST 4: Cluster Bias Validation")
        print("="*70)

        z_test = 0.3

        # Calculate bias for different richness bins
        print(f"\nCluster bias at z={z_test}:")
        biases = []
        for i, (lmin, lmax) in enumerate(calculator.lambda_bins):
            b = calculator.cluster_bias(z_test, lmin, lmax)
            biases.append(b)
            print(f"  λ = {lmin:3d}-{lmax:3d}: b = {b:.3f}")

        # Physical checks
        # 1. All biases should be > 1 (clusters are biased tracers)
        all_above_one = all(b > 1.0 for b in biases)

        # 2. Bias should increase with richness (higher mass → higher bias)
        increasing_trend = all(biases[i] < biases[i+1] for i in range(len(biases)-1))

        # 3. Bias should be in reasonable range (1.5 to 5.0 for clusters)
        reasonable_range = all(1.5 < b < 5.0 for b in biases)

        print(f"\nPhysical checks:")
        print(f"  All b > 1.0 (biased tracers): {'✅ YES' if all_above_one else '❌ NO'}")
        print(f"  Increasing with richness: {'✅ YES' if increasing_trend else '❌ NO'}")
        print(f"  In reasonable range [1.5, 5.0]: {'✅ YES' if reasonable_range else '❌ NO'}")

        if all_above_one and increasing_trend and reasonable_range:
            print(f"\n✅ PASS: Cluster biases are physically reasonable")
            self.results['bias'] = 'PASS'
        else:
            print(f"\n❌ FAIL: Cluster biases have issues")
            self.results['bias'] = 'FAIL'

        return all_above_one and increasing_trend and reasonable_range

    def test_covariance_properties(self, calculator):
        """Test covariance matrix properties"""
        print("\n" + "="*70)
        print("TEST 5: Covariance Matrix Properties")
        print("="*70)

        print("\nThis test will calculate a small covariance matrix (2x2 bins)...")
        print("(Full matrix calculation would be very slow)\n")

        # Calculate for just 2 bins in same z-bin
        z_min, z_max = 0.2, 0.4
        lambda1 = (20, 30)
        lambda2 = (30, 45)

        # Get counts
        N1 = calculator.cluster_number_count(z_min, z_max, *lambda1)
        N2 = calculator.cluster_number_count(z_min, z_max, *lambda2)

        # Get biases
        z_center = 0.5 * (z_min + z_max)
        b1 = calculator.cluster_bias(z_center, *lambda1)
        b2 = calculator.cluster_bias(z_center, *lambda2)

        # Get SSV
        sigma_b_sq = calculator.super_sample_variance(z_center)

        # Calculate covariance elements
        vol_weight = calculator.comoving_volume_element(z_center) * (z_max - z_min)

        # Diagonal elements
        SSV_11 = calculator.Omega_s**2 * vol_weight * b1 * b1 * sigma_b_sq
        SSV_22 = calculator.Omega_s**2 * vol_weight * b2 * b2 * sigma_b_sq
        Cov_11 = N1 + SSV_11
        Cov_22 = N2 + SSV_22

        # Off-diagonal
        SSV_12 = calculator.Omega_s**2 * vol_weight * b1 * b2 * sigma_b_sq
        Cov_12 = SSV_12

        # Build mini matrix
        C = np.array([[Cov_11, Cov_12],
                      [Cov_12, Cov_22]])

        print(f"2x2 Covariance matrix:")
        print(f"  C_11 = {Cov_11:.3e} (shot: {N1:.3e}, SSV: {SSV_11:.3e})")
        print(f"  C_12 = {Cov_12:.3e}")
        print(f"  C_22 = {Cov_22:.3e} (shot: {N2:.3e}, SSV: {SSV_22:.3e})")

        # Tests
        # 1. Positive definite (all eigenvalues > 0)
        eigenvalues = np.linalg.eigvalsh(C)
        is_positive_definite = all(eig > 0 for eig in eigenvalues)

        # 2. Diagonal dominated by shot noise
        shot_dominates = (SSV_11 / N1 < 0.1) and (SSV_22 / N2 < 0.1)

        # 3. Correlations are small but non-zero
        corr_12 = Cov_12 / np.sqrt(Cov_11 * Cov_22)
        corr_small = 1e-6 < abs(corr_12) < 0.1

        print(f"\nMatrix properties:")
        print(f"  Positive definite: {'✅ YES' if is_positive_definite else '❌ NO'}")
        print(f"  Eigenvalues: {eigenvalues}")
        print(f"  Shot noise dominates diagonal: {'✅ YES' if shot_dominates else '❌ NO'}")
        print(f"  SSV/Shot ratios: {SSV_11/N1:.2e}, {SSV_22/N2:.2e}")
        print(f"  Correlation ρ_12 = {corr_12:.6f}")
        print(f"  Correlation in range [10^-6, 0.1]: {'✅ YES' if corr_small else '❌ NO'}")

        if is_positive_definite and shot_dominates:
            print(f"\n✅ PASS: Covariance matrix properties are correct")
            self.results['covariance'] = 'PASS'
        else:
            print(f"\n❌ FAIL: Covariance matrix has issues")
            self.results['covariance'] = 'FAIL'

        return is_positive_definite and shot_dominates

    def print_final_report(self):
        """Print final validation report"""
        print("\n" + "="*70)
        print("FINAL VALIDATION REPORT")
        print("="*70)

        total_tests = len(self.results)
        passed = sum(1 for v in self.results.values() if v == 'PASS')
        warnings = sum(1 for v in self.results.values() if v == 'WARNING')
        failed = sum(1 for v in self.results.values() if v == 'FAIL')

        print(f"\nTest Results:")
        for test_name, result in self.results.items():
            symbol = {
                'PASS': '✅',
                'FAIL': '❌',
                'WARNING': '⚠️'
            }.get(result, '?')
            print(f"  {symbol} {test_name}: {result}")

        print(f"\nSummary:")
        print(f"  Total tests: {total_tests}")
        print(f"  Passed: {passed}")
        print(f"  Warnings: {warnings}")
        print(f"  Failed: {failed}")

        if failed == 0 and warnings == 0:
            print(f"\n🎉 ALL TESTS PASSED - Implementation is CosmoLike-compatible!")
        elif failed == 0:
            print(f"\n✅ All critical tests passed (some warnings)")
        else:
            print(f"\n⚠️  Some tests failed - review implementation")


def main():
    """Run validation tests"""
    print("="*70)
    print("COSMOLIKE COMPATIBILITY VALIDATION TESTS")
    print("="*70)
    print("\nTesting: cosmolike_compatible_cluster_covariance_Y1.py")
    print(f"Against: CosmoLike official specifications")
    print(f"Reference: Krause & Eifler (2016), CosmoLike GitHub repositories")

    validator = CosmoLikeValidator()

    # Test Y1 configuration
    print("\n" + "#"*70)
    print("# TESTING Y1 CONFIGURATION")
    print("#"*70)

    calc_y1 = CosmoLikeClusterCovariance(survey_config="Y1")

    validator.test_binning_scheme(calc_y1)
    validator.test_survey_configuration(calc_y1, "Y1")
    validator.test_cluster_counts_realism(calc_y1)
    validator.test_bias_values(calc_y1)
    validator.test_covariance_properties(calc_y1)

    # Test Y10 configuration
    print("\n" + "#"*70)
    print("# TESTING Y10 CONFIGURATION")
    print("#"*70)

    calc_y10 = CosmoLikeClusterCovariance(survey_config="Y10")

    validator.test_binning_scheme(calc_y10)
    validator.test_survey_configuration(calc_y10, "Y10")

    # Print final report
    validator.print_final_report()

    print("\n" + "="*70)
    print("VALIDATION COMPLETE")
    print("="*70)
    print("\nFor detailed analysis, see: COSMOLIKE_VERIFICATION.md")
    print("For equation verification, see: FINAL_CORRECTIONS.md")


if __name__ == "__main__":
    main()
