#!/usr/bin/env python3
"""
Mass Bin Covariance Calculator

Extends the richness-bin covariance implementation to mass bins.
Compares diagonal elements with error bars from ../hmf_err7.py

Based on:
- cov_mat_cosmolike.py (our corrected implementation)
- ../hmf_err7.py (mass bin cluster counts and errors)
- Krause & Eifler (2016) equations
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import j1
import sys
import os

# Add parent directory to path to import hmf_err7 values
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


class MassBinCovariance:
    """
    Calculate covariance matrix for cluster number counts in MASS bins
    (as opposed to richness bins)
    """

    def __init__(self, survey_area=11808.0):
        """
        Initialize with survey parameters matching hmf_err7.py

        Parameters:
        -----------
        survey_area : float
            Survey area in square degrees (11808 from hmf_err7.py)
        """
        self.survey_area = survey_area
        self.Omega_s = survey_area * (np.pi/180.0)**2  # steradians
        self.h = 0.6727  # Match CosmoLike

        # Cosmological parameters (Krause & Eifler Table 1)
        self.Omega_m = 0.3156
        self.sigma_8 = 0.831
        self.n_s = 0.9645
        self.w_0 = -1.0
        self.w_a = 0.0

        # Mass bins from hmf_err7.py
        self.masses = np.array([1e14, 3e14, 5e14, 7e14, 1e15, 3e15])
        self.n_mass_bins = len(self.masses)

        # Calculate mass bin edges and widths in log space
        self.lnM = np.log(self.masses)
        self.dlnM = np.empty_like(self.lnM)
        self.dlnM[0] = self.lnM[1] - self.lnM[0]
        self.dlnM[-1] = self.lnM[-1] - self.lnM[-2]
        self.dlnM[1:-1] = 0.5 * (self.lnM[2:] - self.lnM[:-2])

        # Redshift bins from hmf_err7.py
        self.z_bins = [(0.2, 0.4), (0.4, 0.6), (0.6, 0.8)]
        self.n_z_bins = len(self.z_bins)

        # Mass-richness scatter (CosmoLike baseline)
        self.sigma_ln_M_lambda = 0.25

        print(f"Initialized mass bin covariance calculator:")
        print(f"  Survey area: {survey_area:.1f} deg² ({self.Omega_s:.4f} sr)")
        print(f"  Mass bins: {self.n_mass_bins}")
        print(f"  Mass range: {self.masses[0]:.1e} - {self.masses[-1]:.1e} M☉/h")
        print(f"  Redshift bins: {self.n_z_bins}")
        print(f"  Total bins: {self.n_mass_bins * self.n_z_bins}")

    def hubble_distance_function(self, z):
        """Comoving distance element dχ/dz in Mpc/h"""
        E_z = np.sqrt(self.Omega_m * (1 + z)**3 +
                     (1 - self.Omega_m) * (1 + z)**(3 * (1 + self.w_0 + self.w_a)) *
                     np.exp(-3 * self.w_a * z / (1 + z)))
        c_over_H0 = 2997.92458  # c/H0 in Mpc/h
        return c_over_H0 / E_z

    def comoving_distance(self, z):
        """Comoving distance χ(z) in Mpc/h"""
        z_array = np.linspace(0, z, 1000)
        integrand = self.hubble_distance_function(z_array)
        return np.trapz(integrand, z_array)

    def comoving_volume(self, z_min, z_max):
        """
        Comoving volume between z_min and z_max in (Mpc/h)³

        Matches the calculation in hmf_err7.py:
        V = Ω_s * (χ_max³ - χ_min³) / 3
        """
        chi_min = self.comoving_distance(z_min)
        chi_max = self.comoving_distance(z_max)

        volume = self.Omega_s * (chi_max**3 - chi_min**3) / 3.0
        return volume

    def comoving_volume_element(self, z):
        """Comoving volume element d²V/(dz dΩ) in (Mpc/h)³/sr"""
        chi = self.comoving_distance(z)
        dchi_dz = self.hubble_distance_function(z)
        return chi**2 * dchi_dz

    def tinker_mass_function(self, M, z):
        """
        Tinker et al. (2008) mass function dn/dlnM in (Mpc/h)^-3

        Returns number density per unit ln(M), matching hmf_err7.py format
        """
        # Tinker 2008 parameters for Δ = 200ρ_mean
        A = 0.186 * (1 + z)**(-0.14)
        a = 1.47 * (1 + z)**(-0.06)
        b = 2.57 * (1 + z)**(-0.17)
        c = 1.19

        # Matter density
        rho_crit0 = 2.775e11 * self.h**2  # (M_solar/h)/(Mpc/h)³
        rho_m0 = rho_crit0 * self.Omega_m
        rho_m_z = rho_m0 * (1 + z)**3

        # sigma(M) - CORRECTED version (not the broken Y1 version!)
        M_pivot = 3.0e14  # Cluster mass scale
        sigma_pivot = 0.8  # Calibrated to CAMB/CLASS
        sigma_M = sigma_pivot * (M / M_pivot)**(-0.5)

        # Growth factor
        D_z = 1.0 / (1 + z)
        sigma_M *= D_z
        sigma_M = np.maximum(sigma_M, 0.05)

        # Peak height
        delta_c = 1.686
        nu = delta_c / sigma_M

        # Tinker mass function
        f_nu = A * ((sigma_M/b)**(-a) + 1) * np.exp(-c/sigma_M**2)

        # dn/dlnM (not dn/dM!)
        dlnsigma_dlnM = -0.5
        dn_dlnM = (rho_m_z / M) * abs(dlnsigma_dlnM) * f_nu

        return dn_dlnM

    def tinker_halo_bias(self, M, z):
        """Tinker et al. (2010) halo bias"""
        # Use same sigma(M) as mass function
        M_pivot = 3.0e14
        sigma_pivot = 0.8
        sigma_M = sigma_pivot * (M / M_pivot)**(-0.5)

        D_z = 1.0 / (1 + z)
        sigma_M *= D_z
        sigma_M = np.maximum(sigma_M, 0.05)

        delta_c = 1.686
        nu = delta_c / sigma_M

        # Tinker 2010 bias parameters
        y = np.log10(200.0)
        A = 1.0 + 0.24 * y * np.exp(-(4.0/y)**4)
        a = 0.44 * y - 0.88
        B = 0.183
        b = 1.5
        C = 0.019 + 0.107 * y + 0.19 * np.exp(-(4.0/y)**4)
        c = 2.4

        bias = 1.0 + (A * nu**a + B * nu**b + C * nu**c) / (delta_c * D_z)
        return bias

    def cluster_counts_in_mass_bin(self, z_min, z_max, M_center, dlnM):
        """
        Expected cluster count in mass bin

        N = ∫ (dn/dlnM) dlnM dV
          = (dn/dlnM)|_{M_center} × dlnM × V

        Parameters:
        -----------
        z_min, z_max : float
            Redshift bin boundaries
        M_center : float
            Central mass of bin (M_solar/h)
        dlnM : float
            Width of mass bin in ln(M)

        Returns:
        --------
        N : float
            Expected cluster count in this mass-redshift bin
        """
        # Use bin center for evaluation
        z_center = 0.5 * (z_min + z_max)

        # Get dn/dlnM at this mass and redshift
        dn_dlnM = self.tinker_mass_function(M_center, z_center)

        # Get comoving volume
        volume = self.comoving_volume(z_min, z_max)

        # Cluster count: (dn/dlnM) × dlnM × V
        N = dn_dlnM * dlnM * volume

        return N

    def linear_power_spectrum(self, k, z):
        """Linear matter power spectrum (simplified)"""
        D_z = 1.0 / (1 + z)
        k_h = k * self.h
        P_k = (k_h)**self.n_s * np.exp(-(k_h / 0.1)**2)
        P_k *= (self.sigma_8 * D_z)**2
        return P_k

    def super_sample_variance(self, z):
        """Super-sample variance σ_b²(Ω_s; z)"""
        theta_s = np.sqrt(self.Omega_s / np.pi)
        chi_z = self.comoving_distance(z)

        def integrand(k_perp):
            x = k_perp * chi_z * theta_s
            if x < 1e-6:
                window_sq = 1.0
            else:
                window_sq = (2 * j1(x) / x)**2

            P_lin = self.linear_power_spectrum(k_perp, z)
            # CORRECTED: proper 2D integration measure
            return k_perp * P_lin * window_sq / (2 * np.pi)

        k_max = 10.0  # h/Mpc
        k_array = np.logspace(-4, np.log10(k_max), 1000)
        integrand_values = [integrand(k) for k in k_array]

        sigma_b_sq = np.trapz(integrand_values, k_array)
        return sigma_b_sq

    def calculate_covariance_matrix(self):
        """
        Calculate covariance matrix for cluster counts in mass bins

        Cov[N_i, N_j] = δ_ij N_i + SSV contribution

        where i,j index (mass bin, redshift bin) combinations
        """
        print("\nCalculating cluster counts in mass bins...")

        # Calculate all cluster counts
        N_counts = {}
        bins = []

        for i_z, (z_min, z_max) in enumerate(self.z_bins):
            z_label = f"z={z_min:.1f}-{z_max:.1f}"
            print(f"\n{z_label}:")

            for i_M, (M, dlnM) in enumerate(zip(self.masses, self.dlnM)):
                bin_id = (i_z, i_M)
                bins.append(bin_id)

                N = self.cluster_counts_in_mass_bin(z_min, z_max, M, dlnM)
                N_counts[bin_id] = N

                print(f"  M={M:.1e}: N={N:.2f}")

        n_total_bins = len(bins)
        print(f"\nTotal bins: {n_total_bins}")

        # Calculate covariance matrix
        print("\nCalculating covariance matrix...")

        covariance_matrix = np.zeros((n_total_bins, n_total_bins))

        for idx1, bin1 in enumerate(bins):
            i_z1, i_M1 = bin1
            z_min1, z_max1 = self.z_bins[i_z1]
            M1 = self.masses[i_M1]

            for idx2, bin2 in enumerate(bins):
                i_z2, i_M2 = bin2
                z_min2, z_max2 = self.z_bins[i_z2]
                M2 = self.masses[i_M2]

                # Shot noise (Poisson)
                if bin1 == bin2:
                    shot_noise = N_counts[bin1]
                else:
                    shot_noise = 0.0

                # Super-sample variance (only within same z-bin)
                if i_z1 == i_z2:
                    z_center = 0.5 * (z_min1 + z_max1)

                    # Halo biases
                    b1 = self.tinker_halo_bias(M1, z_center)
                    b2 = self.tinker_halo_bias(M2, z_center)

                    # Super-sample variance
                    sigma_b_sq = self.super_sample_variance(z_center)

                    # Volume weight
                    vol_weight = self.comoving_volume_element(z_center) * (z_max1 - z_min1)

                    # SSV contribution
                    super_sample_var = (self.Omega_s**2 * vol_weight *
                                      b1 * b2 * sigma_b_sq)
                else:
                    super_sample_var = 0.0

                # Total covariance
                covariance_matrix[idx1, idx2] = shot_noise + super_sample_var

        # Store results
        self.bins = bins
        self.N_counts = N_counts
        self.covariance_matrix = covariance_matrix

        return covariance_matrix

    def get_correlation_matrix(self):
        """Convert covariance to correlation matrix"""
        if not hasattr(self, 'covariance_matrix'):
            raise ValueError("Must calculate covariance matrix first!")

        diag = np.diag(self.covariance_matrix)
        sqrt_diag = np.sqrt(np.maximum(diag, 1e-20))
        correlation_matrix = self.covariance_matrix / np.outer(sqrt_diag, sqrt_diag)

        return correlation_matrix

    def compare_with_hmf_err7(self):
        """
        Compare diagonal covariance elements with error bars from hmf_err7.py

        The diagonal elements Cov[N_i, N_i] = σ²(N_i) should match
        the error bars squared from hmf_err7.py
        """
        print("\n" + "="*70)
        print("COMPARISON WITH HMF_ERR7.PY ERROR BARS")
        print("="*70)

        # Import values from hmf_err7.py
        try:
            # Mass bins and HMF values from hmf_err7.py
            masses_hmf = np.array([1e14, 3e14, 5e14, 7e14, 1e15, 3e15])

            # HMF values (dn/dlnM in (Mpc/h)^-3)
            hmf_z1 = np.array([4e-6, 2e-6, 1.2e-6, 7e-7, 4e-7, 5e-8])
            hmf_z2 = np.array([3e-6, 1.5e-6, 8e-7, 5e-7, 2e-7, 2e-8])
            hmf_z3 = np.array([2e-6, 1e-6, 6e-7, 3e-7, 1e-7, 1e-8])

            # Get volumes (need to recalculate from hmf_err7.py)
            vol_z1 = self.comoving_volume(0.2, 0.4)
            vol_z2 = self.comoving_volume(0.4, 0.6)
            vol_z3 = self.comoving_volume(0.6, 0.8)

            # Calculate counts from HMF: N = (dn/dlnM) × dlnM × V
            counts_hmf_z1 = hmf_z1 * self.dlnM * vol_z1
            counts_hmf_z2 = hmf_z2 * self.dlnM * vol_z2
            counts_hmf_z3 = hmf_z3 * self.dlnM * vol_z3

            # Poisson errors on counts
            err_counts_z1 = np.sqrt(counts_hmf_z1)
            err_counts_z2 = np.sqrt(counts_hmf_z2)
            err_counts_z3 = np.sqrt(counts_hmf_z3)

            # Store for plotting
            self.hmf_values = {
                'masses': masses_hmf,
                'z1': {'hmf': hmf_z1, 'counts': counts_hmf_z1, 'err_counts': err_counts_z1},
                'z2': {'hmf': hmf_z2, 'counts': counts_hmf_z2, 'err_counts': err_counts_z2},
                'z3': {'hmf': hmf_z3, 'counts': counts_hmf_z3, 'err_counts': err_counts_z3}
            }

            # Extract diagonal elements from covariance matrix
            diag_cov = np.diag(self.covariance_matrix)

            # Reshape by redshift bins
            sigma_from_cov = np.sqrt(diag_cov)

            print("\nComparison: sqrt(Cov_ii) vs Poisson errors")
            print("-"*70)

            # Compare bin by bin
            for i_z, z_label in enumerate(['z=0.2-0.4', 'z=0.4-0.6', 'z=0.6-0.8']):
                print(f"\n{z_label}:")
                print(f"  {'M [M☉/h]':<15} {'N_cov':<12} {'σ_cov':<12} {'σ_Pois':<12} {'Ratio':<10}")

                z_key = f'z{i_z+1}'
                err_poisson = self.hmf_values[z_key]['err_counts']
                counts_ref = self.hmf_values[z_key]['counts']

                for i_M in range(self.n_mass_bins):
                    idx = i_z * self.n_mass_bins + i_M

                    N_from_cov = self.N_counts[(i_z, i_M)]
                    sigma_cov = sigma_from_cov[idx]
                    sigma_pois = err_poisson[i_M]

                    ratio = sigma_cov / sigma_pois if sigma_pois > 0 else np.inf

                    print(f"  {self.masses[i_M]:<15.1e} {N_from_cov:<12.2f} "
                          f"{sigma_cov:<12.2f} {sigma_pois:<12.2f} {ratio:<10.3f}")

            # Store for plotting
            self.sigma_from_cov = sigma_from_cov.reshape(self.n_z_bins, self.n_mass_bins)
            self.err_poisson_hmf = np.array([
                err_counts_z1, err_counts_z2, err_counts_z3
            ])

            print("\n" + "="*70)
            return True

        except Exception as e:
            print(f"Error comparing with hmf_err7.py: {e}")
            return False

    def plot_comparison(self, save_path=None):
        """Create comparison plots"""
        if not hasattr(self, 'sigma_from_cov'):
            print("Run compare_with_hmf_err7() first!")
            return

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        # Panel A: Cluster counts
        ax = axes[0, 0]
        for i_z, z_label in enumerate(['z=0.2-0.4', 'z=0.4-0.6', 'z=0.6-0.8']):
            counts = [self.N_counts[(i_z, i_M)] for i_M in range(self.n_mass_bins)]
            ax.plot(self.masses, counts, 'o-', label=z_label, markersize=6)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$M\ [h^{-1} M_\odot]$', fontsize=12)
        ax.set_ylabel('Expected Cluster Count $N$', fontsize=12)
        ax.set_title('A. Cluster Counts in Mass Bins', fontsize=13, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Panel B: Error comparison
        ax = axes[0, 1]
        colors = ['C0', 'C1', 'C2']
        for i_z, (z_label, color) in enumerate(zip(
            ['z=0.2-0.4', 'z=0.4-0.6', 'z=0.6-0.8'], colors)):

            sigma_cov = self.sigma_from_cov[i_z, :]
            sigma_pois = self.err_poisson_hmf[i_z, :]

            ax.plot(self.masses, sigma_cov, 'o-', color=color,
                   label=f'{z_label} (Covariance)', markersize=6)
            ax.plot(self.masses, sigma_pois, 's--', color=color, alpha=0.6,
                   label=f'{z_label} (Poisson)', markersize=5)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$M\ [h^{-1} M_\odot]$', fontsize=12)
        ax.set_ylabel(r'Error $\sigma_N$', fontsize=12)
        ax.set_title('B. Error Comparison: √Cov vs Poisson', fontsize=13, fontweight='bold')
        ax.legend(fontsize=9, ncol=2)
        ax.grid(True, alpha=0.3)

        # Panel C: Ratio
        ax = axes[1, 0]
        for i_z, (z_label, color) in enumerate(zip(
            ['z=0.2-0.4', 'z=0.4-0.6', 'z=0.6-0.8'], colors)):

            sigma_cov = self.sigma_from_cov[i_z, :]
            sigma_pois = self.err_poisson_hmf[i_z, :]
            ratio = sigma_cov / sigma_pois

            ax.plot(self.masses, ratio, 'o-', color=color, label=z_label, markersize=6)

        ax.axhline(1.0, color='k', linestyle='--', alpha=0.5, label='Perfect match')
        ax.set_xscale('log')
        ax.set_xlabel(r'$M\ [h^{-1} M_\odot]$', fontsize=12)
        ax.set_ylabel(r'$\sigma_{\mathrm{Cov}} / \sigma_{\mathrm{Poisson}}$', fontsize=12)
        ax.set_title('C. Ratio: Covariance Error / Poisson Error', fontsize=13, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0.5, 2.0)

        # Panel D: Correlation matrix
        ax = axes[1, 1]
        corr_matrix = self.get_correlation_matrix()

        from matplotlib.colors import LogNorm
        im = ax.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        plt.colorbar(im, ax=ax, label='Correlation')

        # Add labels
        n_bins = len(self.bins)
        tick_positions = np.arange(0, n_bins, 3)
        ax.set_xticks(tick_positions)
        ax.set_yticks(tick_positions)
        ax.set_xlabel('Bin Index', fontsize=12)
        ax.set_ylabel('Bin Index', fontsize=12)
        ax.set_title('D. Correlation Matrix (Mass Bins)', fontsize=13, fontweight='bold')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"\nPlot saved to: {save_path}")

        plt.show()

        return fig

    def print_summary(self):
        """Print summary statistics"""
        if not hasattr(self, 'covariance_matrix'):
            print("No covariance matrix calculated yet!")
            return

        print("\n" + "="*70)
        print("MASS BIN COVARIANCE SUMMARY")
        print("="*70)

        print(f"\nSurvey parameters:")
        print(f"  Area: {self.survey_area:.1f} deg²")
        print(f"  Mass bins: {self.n_mass_bins}")
        print(f"  Redshift bins: {self.n_z_bins}")
        print(f"  Total bins: {len(self.bins)}")

        print(f"\nTotal expected clusters:")
        total = sum(self.N_counts.values())
        print(f"  {total:.1f} clusters")

        print(f"\nCovariance matrix properties:")
        diag = np.diag(self.covariance_matrix)
        print(f"  Diagonal range: {np.min(diag):.2e} to {np.max(diag):.2e}")

        corr_matrix = self.get_correlation_matrix()
        off_diag = corr_matrix[np.triu_indices_from(corr_matrix, k=1)]
        print(f"  Correlation range: {np.min(off_diag):.3f} to {np.max(off_diag):.3f}")


def main():
    """Main execution"""
    print("="*70)
    print("MASS BIN COVARIANCE CALCULATOR")
    print("="*70)
    print("\nExtending richness-bin covariance to mass bins")
    print("Comparing with error bars from ../hmf_err7.py\n")

    # Initialize calculator (use same survey area as hmf_err7.py)
    calc = MassBinCovariance(survey_area=11808.0)

    # Calculate covariance matrix
    cov_matrix = calc.calculate_covariance_matrix()

    # Print summary
    calc.print_summary()

    # Compare with hmf_err7.py
    success = calc.compare_with_hmf_err7()

    if success:
        # Create comparison plots
        save_path = '/Users/nesar/Projects/HEP/EDE/Codes/CovMat/mass_bin_comparison.png'
        calc.plot_comparison(save_path=save_path)

        # Save numerical results
        np.savetxt('/Users/nesar/Projects/HEP/EDE/Codes/CovMat/mass_bin_covariance.txt',
                   cov_matrix,
                   header='Mass bin covariance matrix (matching hmf_err7.py)')

        corr_matrix = calc.get_correlation_matrix()
        np.savetxt('/Users/nesar/Projects/HEP/EDE/Codes/CovMat/mass_bin_correlation.txt',
                   corr_matrix,
                   header='Mass bin correlation matrix')

        print("\nFiles created:")
        print("  - mass_bin_comparison.png")
        print("  - mass_bin_covariance.txt")
        print("  - mass_bin_correlation.txt")

    print("\n" + "="*70)
    print("DONE")
    print("="*70)


if __name__ == "__main__":
    main()
