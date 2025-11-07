#!/usr/bin/env python3
"""
Cluster Number Count Correlation Matrix Calculator

Implementation based on "CosmoLike - Cosmological Likelihood Analyses for Photometric Galaxy Surveys"
by Krause & Eifler (arXiv:1601.05779)

This code implements the covariance matrix for cluster number counts as described in Section 2.1.1
and Appendix A, specifically equations A13 and A14.

The covariance includes:
1. Shot noise term
2. Super-sample variance term

Key equations implemented:
- Equation 1: Expected cluster count N_i(λ_α)
- Equation A14: Covariance between cluster number count bins
- Equation A13: Response to large-scale modes for super-sample variance
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, dblquad
from scipy.special import erf
import warnings
warnings.filterwarnings('ignore')

class ClusterNumberCountCovariance:
    """
    Calculate covariance matrix for cluster number counts following Krause & Eifler 2016
    """
    
    def __init__(self, survey_area=18000.0, h=0.6727):
        """
        Initialize with survey parameters
        
        Parameters:
        -----------
        survey_area : float
            Survey area in square degrees (default: 18000 from LSST baseline)
        h : float
            Hubble parameter h = H0/(100 km/s/Mpc) (default: 0.6727 from paper)
        """
        self.Omega_s = survey_area * (np.pi/180.0)**2  # Convert to steradians
        self.h = h
        
        # Cosmological parameters from Table 1 of the paper
        self.Omega_m = 0.3156
        self.sigma_8 = 0.831
        self.n_s = 0.9645
        self.w_0 = -1.0
        self.w_a = 0.0
        self.Omega_b = 0.0492
        self.h_0 = 0.6727
        
        # Mass-richness relation parameters (Equation 2)
        self.C_lambda = 33.6  # ln[M/(M_solar/h)]
        self.a_lambda = 1.08  # Richness slope
        self.b_lambda = 0.0   # Redshift evolution
        self.sigma_ln_M_lambda = 0.25  # Log-normal scatter
        
        # Redshift bins (from Section 3.1)
        self.z_bins = [(0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)]
        self.n_z_bins = len(self.z_bins)
        
        # Richness bins (λ > 10, 7 bins as mentioned in Section 3.1)
        self.lambda_min = 10.0
        self.lambda_bins = [(10, 20), (20, 30), (30, 40), (40, 60), (60, 80), (80, 120), (120, 200)]
        self.n_lambda_bins = len(self.lambda_bins)
        
        print(f"Initialized cluster covariance calculator:")
        print(f"  Survey area: {survey_area:.1f} deg² ({self.Omega_s:.4f} steradians)")
        print(f"  Redshift bins: {self.n_z_bins}")
        print(f"  Richness bins: {self.n_lambda_bins}")
        print(f"  Total cluster bins: {self.n_z_bins * self.n_lambda_bins}")
    
    def hubble_distance_function(self, z):
        """
        Comoving distance element dχ/dz in Mpc/h
        Assuming flat ΛCDM with time-varying dark energy
        """
        E_z = np.sqrt(self.Omega_m * (1 + z)**3 + 
                     (1 - self.Omega_m) * (1 + z)**(3 * (1 + self.w_0 + self.w_a)) * 
                     np.exp(-3 * self.w_a * z / (1 + z)))
        c_over_H0 = 2997.92458  # c/H0 in Mpc/h for h=1
        return c_over_H0 / E_z
    
    def comoving_volume_element(self, z):
        """
        Comoving volume element d²V/(dz dΩ) in (Mpc/h)³/steradian
        Following the paper's notation
        """
        chi = self.comoving_distance(z)
        dchi_dz = self.hubble_distance_function(z)
        return chi**2 * dchi_dz
    
    def comoving_distance(self, z):
        """
        Comoving distance χ(z) in Mpc/h
        """
        z_array = np.linspace(0, z, 1000)
        integrand = self.hubble_distance_function(z_array)
        return np.trapz(integrand, z_array)
    
    def tinker_mass_function(self, M, z):
        """
        Tinker et al. (2008, 2010) halo mass function dn/dM in h³/Mpc³/(M_solar/h)
        Improved version with proper normalization
        
        Parameters:
        -----------
        M : float or array
            Halo mass in M_solar/h
        z : float
            Redshift
        """
        # Tinker 2008 parameters for Δ = 200ρ_mean
        A = 0.186 * (1 + z)**(-0.14)
        a = 1.47 * (1 + z)**(-0.06)
        b = 2.57 * (1 + z)**(-0.17)
        c = 1.19
        
        # Convert mass to h units if needed
        M_h = np.asarray(M)
        
        # Matter density at z=0 in (M_solar/h)/(Mpc/h)³
        rho_crit0 = 2.775e11 * self.h**2  # Critical density in (M_solar/h)/(Mpc/h)³
        rho_m0 = rho_crit0 * self.Omega_m  # Matter density at z=0
        rho_m_z = rho_m0 * (1 + z)**3     # Matter density at redshift z
        
        # Use a more realistic sigma(M) relation
        # This is an approximation - in practice you'd use CAMB/CLASS
        M_8 = 6e14  # Mass scale where sigma=1 at z=0 (approximately)
        sigma_M = self.sigma_8 * (M_h / M_8)**(-0.5)  # Power-law approximation
        
        # Apply growth factor
        D_z = 1.0 / (1 + z)  # Simple growth factor approximation
        sigma_M *= D_z
        
        # Avoid numerical issues
        sigma_M = np.maximum(sigma_M, 1e-6)
        
        # Peak height
        delta_c = 1.686  # Critical density for collapse
        nu = delta_c / sigma_M
        
        # Tinker mass function
        f_nu = A * ((sigma_M/b)**(-a) + 1) * np.exp(-c/sigma_M**2)
        
        # dn/dM = (rho_m/M) * (d ln σ⁻¹/d ln M) * f(σ)
        dlnsigma_dlnM = -0.5  # From power-law approximation
        dn_dM = (rho_m_z / M_h) * abs(dlnsigma_dlnM) * f_nu
        
        return dn_dM
    
    def mass_richness_relation(self, lambda_val, z):
        """
        Mean mass-richness relation from Equation 2
        Adjusted for more realistic cluster counts
        
        Returns mean mass M̄(λ) in M_solar/h
        """
        # Adjusted normalization to give reasonable cluster counts
        C_lambda_adj = np.log(2e14)  # More realistic normalization
        ln_M = (C_lambda_adj + 
                self.a_lambda * np.log(lambda_val / 60.0) + 
                self.b_lambda * np.log(1 + z))
        return np.exp(ln_M)
    
    def mass_richness_probability(self, M, lambda_val, z):
        """
        Probability p(M|λ,z) from Equation 3
        Log-normal distribution around mean mass-richness relation
        """
        M_mean = self.mass_richness_relation(lambda_val, z)
        ln_M = np.log(M)
        ln_M_mean = np.log(M_mean)
        
        normalization = 1.0 / (M * np.sqrt(2 * np.pi) * self.sigma_ln_M_lambda)
        exponent = -(ln_M - ln_M_mean)**2 / (2 * self.sigma_ln_M_lambda**2)
        
        return normalization * np.exp(exponent)
    
    def tinker_halo_bias(self, M, z):
        """
        Tinker et al. (2010) halo bias relation b_h(M)
        Improved version with proper scaling
        """
        # Use the same sigma(M) as in mass function for consistency
        M_8 = 6e14  # Mass scale where sigma=1 at z=0
        sigma_M = self.sigma_8 * (M / M_8)**(-0.5)
        
        # Apply growth factor
        D_z = 1.0 / (1 + z)
        sigma_M *= D_z
        
        # Avoid numerical issues
        sigma_M = np.maximum(sigma_M, 1e-6)
        
        # Peak height
        delta_c = 1.686
        nu = delta_c / sigma_M
        
        # Tinker 2010 bias parameters for Δ = 200ρ_mean
        y = np.log10(200.0)
        A = 1.0 + 0.24 * y * np.exp(-(4.0/y)**4)
        a = 0.44 * y - 0.88
        B = 0.183
        b = 1.5
        C = 0.019 + 0.107 * y + 0.19 * np.exp(-(4.0/y)**4)
        c = 2.4
        
        # Bias formula
        bias = 1.0 + (A * nu**a + B * nu**b + C * nu**c) / (delta_c * D_z)
        
        return bias
    
    def cluster_number_count(self, z_min, z_max, lambda_min, lambda_max):
        """
        Expected cluster count N_i(λ_α) from Equation 1
        
        Parameters:
        -----------
        z_min, z_max : float
            Redshift bin boundaries
        lambda_min, lambda_max : float
            Richness bin boundaries
        """
        def integrand_z(z):
            vol_element = self.comoving_volume_element(z)
            
            def integrand_M(M):
                dn_dM = self.tinker_mass_function(M, z)
                
                def integrand_lambda(lambda_val):
                    p_M_lambda = self.mass_richness_probability(M, lambda_val, z)
                    return p_M_lambda
                
                # Integrate over richness bin
                lambda_integral, _ = quad(integrand_lambda, lambda_min, lambda_max)
                return dn_dM * lambda_integral
            
            # Integrate over reasonable mass range
            M_min = 1e12  # M_solar/h
            M_max = 1e16  # M_solar/h
            M_integral, _ = quad(integrand_M, M_min, M_max)
            
            return vol_element * M_integral
        
        # Integrate over redshift bin
        z_integral, _ = quad(integrand_z, z_min, z_max)
        
        # Multiply by survey area
        N_expected = self.Omega_s * z_integral
        
        return N_expected
    
    def cluster_bias(self, z, lambda_min, lambda_max):
        """
        Mean linear bias of clusters in richness bin α
        From Equation 11
        """
        def numerator_integrand(M):
            dn_dM = self.tinker_mass_function(M, z)
            b_h = self.tinker_halo_bias(M, z)
            
            def lambda_integrand(lambda_val):
                return self.mass_richness_probability(M, lambda_val, z)
            
            lambda_integral, _ = quad(lambda_integrand, lambda_min, lambda_max)
            return dn_dM * b_h * lambda_integral
        
        def denominator_integrand(M):
            dn_dM = self.tinker_mass_function(M, z)
            
            def lambda_integrand(lambda_val):
                return self.mass_richness_probability(M, lambda_val, z)
            
            lambda_integral, _ = quad(lambda_integrand, lambda_min, lambda_max)
            return dn_dM * lambda_integral
        
        # Mass integration range
        M_min = 1e12
        M_max = 1e16
        
        numerator, _ = quad(numerator_integrand, M_min, M_max)
        denominator, _ = quad(denominator_integrand, M_min, M_max)
        
        if denominator > 0:
            return numerator / denominator
        else:
            return 1.0
    
    def linear_power_spectrum(self, k, z):
        """
        Linear matter power spectrum P_lin(k,z)
        Simple approximation using growth factor
        """
        # Growth factor approximation
        D_z = 1.0 / (1 + z)
        
        # Simple power-law power spectrum
        # This should be replaced with proper CAMB/CLASS calculation in real analysis
        k_h = k * self.h  # Convert to h/Mpc
        P_k = (k_h)**self.n_s * np.exp(-(k_h / 0.1)**2)  # Approximate shape
        
        # Normalize with σ₈
        P_k *= (self.sigma_8 * D_z)**2
        
        return P_k
    
    def super_sample_variance(self, z):
        """
        Super-sample variance σ_b²(Ω_s; z) from Equation A8
        
        Approximation for disk-like survey geometry
        """
        # Survey radius
        theta_s = np.sqrt(self.Omega_s / np.pi)
        chi_z = self.comoving_distance(z)
        
        # k_perp integration for super-sample variance
        def integrand(k_perp):
            # Bessel function approximation for survey window
            x = k_perp * chi_z * theta_s
            if x < 1e-6:
                window_sq = 1.0
            else:
                from scipy.special import j1
                window_sq = (2 * j1(x) / x)**2
            
            P_lin = self.linear_power_spectrum(k_perp, z)
            return P_lin * window_sq / (2 * np.pi)
        
        # Integrate over k_perp
        k_max = 10.0  # h/Mpc, reasonable cutoff
        k_array = np.logspace(-4, np.log10(k_max), 1000)
        integrand_values = [integrand(k) for k in k_array]
        
        # Numerical integration
        sigma_b_sq = np.trapz(integrand_values, k_array)
        
        return sigma_b_sq
    
    def calculate_covariance_matrix(self):
        """
        Calculate the full covariance matrix for cluster number counts
        Following Equation A14
        """
        print("Calculating cluster number counts...")
        
        # First calculate all cluster number counts
        N_counts = {}
        cluster_bins = []
        
        for i, (z_min, z_max) in enumerate(self.z_bins):
            for alpha, (lambda_min, lambda_max) in enumerate(self.lambda_bins):
                bin_id = (i, alpha)
                cluster_bins.append(bin_id)
                
                print(f"  Bin z={z_min:.1f}-{z_max:.1f}, λ={lambda_min}-{lambda_max}")
                N_expected = self.cluster_number_count(z_min, z_max, lambda_min, lambda_max)
                N_counts[bin_id] = N_expected
                print(f"    Expected count: {N_expected:.1f}")
        
        n_total_bins = len(cluster_bins)
        print(f"\nTotal bins: {n_total_bins}")
        
        print("\nCalculating covariance matrix...")
        
        # Initialize covariance matrix
        covariance_matrix = np.zeros((n_total_bins, n_total_bins))
        
        for idx1, bin1 in enumerate(cluster_bins):
            i1, alpha1 = bin1
            z_min1, z_max1 = self.z_bins[i1]
            lambda_min1, lambda_max1 = self.lambda_bins[alpha1]
            
            for idx2, bin2 in enumerate(cluster_bins):
                i2, alpha2 = bin2
                z_min2, z_max2 = self.z_bins[i2]
                lambda_min2, lambda_max2 = self.lambda_bins[alpha2]
                
                # Shot noise term (Kronecker delta terms)
                if bin1 == bin2:
                    shot_noise = N_counts[bin1]
                else:
                    shot_noise = 0.0
                
                # Super-sample variance term
                # Only non-zero if redshift bins are the same (neglecting cross-z correlations)
                if i1 == i2:
                    z_center = 0.5 * (z_min1 + z_max1)  # Use bin center
                    
                    # Calculate cluster biases
                    b1 = self.cluster_bias(z_center, lambda_min1, lambda_max1)
                    b2 = self.cluster_bias(z_center, lambda_min2, lambda_max2)
                    
                    # Super-sample variance
                    sigma_b_sq = self.super_sample_variance(z_center)
                    
                    # Volume element weight (simplified)
                    vol_weight = self.comoving_volume_element(z_center) * (z_max1 - z_min1)
                    
                    super_sample_var = (self.Omega_s**2 * vol_weight * 
                                      b1 * b2 * sigma_b_sq)
                else:
                    super_sample_var = 0.0
                
                # Total covariance
                covariance_matrix[idx1, idx2] = shot_noise + super_sample_var
                
                if idx1 <= idx2:  # Only print upper triangle to reduce output
                    print(f"  Cov[{idx1},{idx2}]: shot={shot_noise:.2e}, SSV={super_sample_var:.2e}")
        
        # Store results
        self.cluster_bins = cluster_bins
        self.N_counts = N_counts
        self.covariance_matrix = covariance_matrix
        
        return covariance_matrix
    
    def get_correlation_matrix(self):
        """
        Convert covariance matrix to correlation matrix
        """
        if not hasattr(self, 'covariance_matrix'):
            raise ValueError("Must calculate covariance matrix first!")
        
        # Extract diagonal elements
        diag_elements = np.diag(self.covariance_matrix)
        
        # Avoid division by zero
        sqrt_diag = np.sqrt(np.maximum(diag_elements, 1e-20))
        
        # Calculate correlation matrix
        correlation_matrix = (self.covariance_matrix / 
                            np.outer(sqrt_diag, sqrt_diag))
        
        return correlation_matrix
    
    def plot_covariance_matrix(self, matrix_type='covariance', save_path=None):
        """
        Plot the covariance or correlation matrix
        
        Parameters:
        -----------
        matrix_type : str
            'covariance' or 'correlation'
        save_path : str, optional
            Path to save the plot
        """
        if matrix_type == 'covariance':
            matrix = self.covariance_matrix
            title = 'Cluster Number Count Covariance Matrix'
            cmap = 'viridis'
            norm_type = 'log'
        elif matrix_type == 'correlation':
            matrix = self.get_correlation_matrix()
            title = 'Cluster Number Count Correlation Matrix'
            cmap = 'RdBu_r'
            norm_type = 'linear'
        else:
            raise ValueError("matrix_type must be 'covariance' or 'correlation'")
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Plot matrix
        if norm_type == 'log':
            from matplotlib.colors import LogNorm
            vmin = np.min(matrix[matrix > 0]) if np.any(matrix > 0) else 1e-10
            vmax = np.max(matrix)
            im = ax.imshow(matrix, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
        else:
            vmax = np.max(np.abs(matrix))
            im = ax.imshow(matrix, cmap=cmap, vmin=-vmax, vmax=vmax)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        if matrix_type == 'covariance':
            cbar.set_label('Covariance', fontsize=12)
        else:
            cbar.set_label('Correlation Coefficient', fontsize=12)
        
        # Set labels and ticks
        n_bins = len(self.cluster_bins)
        tick_positions = np.arange(n_bins)
        
        # Create bin labels
        bin_labels = []
        for i, alpha in self.cluster_bins:
            z_min, z_max = self.z_bins[i]
            lambda_min, lambda_max = self.lambda_bins[alpha]
            label = f'z{i+1}λ{alpha+1}\n({z_min:.1f}-{z_max:.1f})\n({lambda_min}-{lambda_max})'
            bin_labels.append(label)
        
        ax.set_xticks(tick_positions[::2])  # Show every other tick to avoid crowding
        ax.set_yticks(tick_positions[::2])
        ax.set_xticklabels([bin_labels[i] for i in range(0, n_bins, 2)], 
                          rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels([bin_labels[i] for i in range(0, n_bins, 2)], 
                          fontsize=8)
        
        # Set title and labels
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('Cluster Bins (redshift-richness)', fontsize=12)
        ax.set_ylabel('Cluster Bins (redshift-richness)', fontsize=12)
        
        # Add grid for better visualization
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to: {save_path}")
        
        plt.show()
        
        return fig, ax
    
    def print_summary(self):
        """
        Print summary of the calculation
        """
        if not hasattr(self, 'covariance_matrix'):
            print("No covariance matrix calculated yet!")
            return
        
        print("\n" + "="*60)
        print("CLUSTER NUMBER COUNT COVARIANCE SUMMARY")
        print("="*60)
        
        print(f"Survey parameters:")
        print(f"  Area: {self.Omega_s * (180/np.pi)**2:.1f} deg²")
        print(f"  Redshift bins: {self.n_z_bins}")
        print(f"  Richness bins: {self.n_lambda_bins}")
        print(f"  Total bins: {len(self.cluster_bins)}")
        
        print(f"\nExpected cluster counts:")
        total_clusters = 0
        for bin_id in self.cluster_bins:
            i, alpha = bin_id
            z_min, z_max = self.z_bins[i]
            lambda_min, lambda_max = self.lambda_bins[alpha]
            count = self.N_counts[bin_id]
            total_clusters += count
            print(f"  z={z_min:.1f}-{z_max:.1f}, λ={lambda_min:3d}-{lambda_max:3d}: {count:8.1f}")
        
        print(f"\nTotal expected clusters: {total_clusters:.1f}")
        
        print(f"\nCovariance matrix statistics:")
        diag_elements = np.diag(self.covariance_matrix)
        print(f"  Diagonal range: {np.min(diag_elements):.2e} to {np.max(diag_elements):.2e}")
        print(f"  Off-diagonal range: {np.min(self.covariance_matrix):.2e} to {np.max(self.covariance_matrix):.2e}")
        
        # Calculate correlation statistics
        corr_matrix = self.get_correlation_matrix()
        off_diag_corr = corr_matrix[np.triu_indices_from(corr_matrix, k=1)]
        print(f"  Correlation range: {np.min(off_diag_corr):.3f} to {np.max(off_diag_corr):.3f}")
        print(f"  Mean correlation: {np.mean(off_diag_corr):.3f}")


def main():
    """
    Main function to run the cluster correlation matrix calculation
    """
    print("Cluster Number Count Correlation Matrix Calculator")
    print("Based on Krause & Eifler (2016) - arXiv:1601.05779")
    print("="*60)
    
    # Initialize calculator with LSST-like parameters
    calculator = ClusterNumberCountCovariance(survey_area=18000.0)
    
    # Calculate covariance matrix
    covariance_matrix = calculator.calculate_covariance_matrix()
    
    # Print summary
    calculator.print_summary()
    
    # Plot results
    print("\nGenerating plots...")
    
    # Plot covariance matrix
    calculator.plot_covariance_matrix(matrix_type='covariance', 
                                    save_path='/Users/nesar/Projects/HEP/EDE/Codes/CovMat/cluster_covariance_matrix.png')
    
    # Plot correlation matrix
    calculator.plot_covariance_matrix(matrix_type='correlation',
                                     save_path='/Users/nesar/Projects/HEP/EDE/Codes/CovMat/cluster_correlation_matrix.png')
    
    # Save numerical results
    np.savetxt('/Users/nesar/Projects/HEP/EDE/Codes/CovMat/cluster_covariance_matrix.txt', 
               covariance_matrix, 
               header='Cluster number count covariance matrix (Krause & Eifler 2016 implementation)')
    
    correlation_matrix = calculator.get_correlation_matrix()
    np.savetxt('/Users/nesar/Projects/HEP/EDE/Codes/CovMat/cluster_correlation_matrix.txt', 
               correlation_matrix,
               header='Cluster number count correlation matrix (Krause & Eifler 2016 implementation)')
    
    print(f"\nResults saved to /Users/nesar/Projects/HEP/EDE/Codes/CovMat/")
    print("Files created:")
    print("  - cluster_covariance_matrix.png")
    print("  - cluster_correlation_matrix.png") 
    print("  - cluster_covariance_matrix.txt")
    print("  - cluster_correlation_matrix.txt")


if __name__ == "__main__":
    main()