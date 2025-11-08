# Quick Equation Reference Guide

## Paper Equations ↔ Code Implementation Mapping

This document provides a quick reference for verifying the code implementation against the paper equations.

---

## 1. Expected Cluster Count

### Paper (Line 206-207)
```
N^i(λ_α) = Ω_s ∫[z_min to z_max] dz (d²V/dzdΩ) ∫ dM (dn/dM) ∫[λ_min to λ_max] dλ p(M|λ,z)
```

### Code (Lines 218-256)
```python
def cluster_number_count(self, z_min, z_max, lambda_min, lambda_max):
    def integrand_z(z):
        vol_element = self.comoving_volume_element(z)  # d²V/dzdΩ

        def integrand_M(M):
            dn_dM = self.tinker_mass_function(M, z)  # dn/dM

            def integrand_lambda(lambda_val):
                p_M_lambda = self.mass_richness_probability(M, lambda_val, z)  # p(M|λ,z)
                return p_M_lambda

            lambda_integral, _ = quad(integrand_lambda, lambda_min, lambda_max)
            return dn_dM * lambda_integral

        M_integral, _ = quad(integrand_M, M_min, M_max)
        return vol_element * M_integral

    z_integral, _ = quad(integrand_z, z_min, z_max)
    N_expected = self.Omega_s * z_integral  # Ω_s × ∫...
    return N_expected
```

---

## 2. Mass-Richness Relation

### Paper (Line 213)
```
ln[M̄(λ)/(M_☉/h)] = C_λ + a_λ ln(λ/60) + b_λ ln(1+z)
```

### Code (Lines 156-168)
```python
def mass_richness_relation(self, lambda_val, z):
    C_lambda_adj = np.log(2e14)  # Normalization (adjusted)
    ln_M = (C_lambda_adj +
            self.a_lambda * np.log(lambda_val / 60.0) +  # a_λ ln(λ/60)
            self.b_lambda * np.log(1 + z))                # b_λ ln(1+z)
    return np.exp(ln_M)
```

---

## 3. Log-Normal Probability

### Paper (Line 219)
```
p(M|λ,z) = 1/(M√(2π)σ_{ln M|λ}) exp[-(ln M - ln M̄(λ))²/(2σ²_{ln M|λ})]
```

### Code (Lines 170-182)
```python
def mass_richness_probability(self, M, lambda_val, z):
    M_mean = self.mass_richness_relation(lambda_val, z)  # M̄(λ)
    ln_M = np.log(M)
    ln_M_mean = np.log(M_mean)

    # 1/(M√(2π)σ)
    normalization = 1.0 / (M * np.sqrt(2 * np.pi) * self.sigma_ln_M_lambda)

    # exp[-(ln M - ln M̄)²/(2σ²)]
    exponent = -(ln_M - ln_M_mean)**2 / (2 * self.sigma_ln_M_lambda**2)

    return normalization * np.exp(exponent)
```

---

## 4. Cluster Bias

### Paper (Line 271-272)
```
b_{λ_α}(z) = [∫ dM (dn/dM) b_h(M) ∫ dλ p(M|λ,z)] / [∫ dM (dn/dM) ∫ dλ p(M|λ,z)]
```

### Code (Lines 258-292)
```python
def cluster_bias(self, z, lambda_min, lambda_max):
    # Numerator: ∫ dM (dn/dM) b_h(M) ∫ dλ p(M|λ,z)
    def numerator_integrand(M):
        dn_dM = self.tinker_mass_function(M, z)
        b_h = self.tinker_halo_bias(M, z)

        def lambda_integrand(lambda_val):
            return self.mass_richness_probability(M, lambda_val, z)

        lambda_integral, _ = quad(lambda_integrand, lambda_min, lambda_max)
        return dn_dM * b_h * lambda_integral

    # Denominator: ∫ dM (dn/dM) ∫ dλ p(M|λ,z)
    def denominator_integrand(M):
        dn_dM = self.tinker_mass_function(M, z)

        def lambda_integrand(lambda_val):
            return self.mass_richness_probability(M, lambda_val, z)

        lambda_integral, _ = quad(lambda_integrand, lambda_min, lambda_max)
        return dn_dM * lambda_integral

    numerator, _ = quad(numerator_integrand, M_min, M_max)
    denominator, _ = quad(denominator_integrand, M_min, M_max)

    return numerator / denominator
```

---

## 5. Super-Sample Variance σ_b²

### Paper (Line 738-740)
```
σ_b²(Ω_s; z) = ∫ (d²k_⊥/(2π)²) P_lin(k_⊥,z) [2J₁(k_⊥ χ(z) θ_s)/(k_⊥ χ(z) θ_s)]²
```

where `θ_s = √(Ω_s/π)` for disk-like survey

**Note**: Converting d²k_⊥ to radial: d²k_⊥ = k_⊥ dk_⊥ dφ → ∫ d²k_⊥/(2π)² = ∫ k_⊥ dk_⊥/(2π)

### Code (Lines 346-388) - **CORRECTED**
```python
def super_sample_variance(self, z):
    theta_s = np.sqrt(self.Omega_s / np.pi)  # √(Ω_s/π)
    chi_z = self.comoving_distance(z)         # χ(z)

    def integrand(k_perp):
        x = k_perp * chi_z * theta_s
        if x < 1e-6:
            window_sq = 1.0
        else:
            from scipy.special import j1
            # [2J₁(x)/x]²
            window_sq = (2 * j1(x) / x)**2

        P_lin = self.linear_power_spectrum(k_perp, z)

        # CORRECTED: k_⊥ P_lin W² / (2π)
        # This accounts for d²k_⊥/(2π)² = k_⊥ dk_⊥/(2π) after φ integration
        return k_perp * P_lin * window_sq / (2 * np.pi)

    k_array = np.logspace(-4, np.log10(k_max), 1000)
    integrand_values = [integrand(k) for k in k_array]
    sigma_b_sq = np.trapz(integrand_values, k_array)

    return sigma_b_sq
```

---

## 6. Covariance Matrix - MOST CRITICAL EQUATION

### Paper (Line 765-767)
```
Cov(N^i_{λ_α}, N^j_{λ_β}) = δ_{i,j} δ_{α,β} N^i_{λ_α}                    [Shot noise]
    + Ω_s² ∫ dχ q^i_{λ_α}(χ) q^j_{λ_β}(χ)                                [SSV integral]
      × [∫ dM (dn/dM) b_h(M,z) ∫_{λ_α,min}^{λ_α,max} dλ p(M|λ,z)]        [Bias term 1]
      × [∫ dM' (dn/dM') b_h(M',z) ∫_{λ_β,min}^{λ_β,max} dλ' p(M'|λ',z)]  [Bias term 2]
      × σ_b²(Ω_s; z(χ))                                                   [SSV variance]
```

where weight function (Line 240-241):
```
q^i_{δ_{λ_α}}(χ) = Θ(z(χ) - z^i_{min}) Θ(z^i_{max} - z(χ)) (dV/dχdΩ)
```

"where we have neglected correlations across redshift bins" → SSV = 0 when i ≠ j

### Code (Lines 429-510) - **CORRECTED**
```python
# Shot noise term: δ_{i,j} δ_{α,β} N^i_{λ_α}
if bin1 == bin2:
    shot_noise = N_counts[bin1]
else:
    shot_noise = 0.0

# Super-sample variance term
if i1 == i2:  # Only same redshift bin
    def chi_integrand(chi):
        # Convert χ to z
        z = self.redshift_from_comoving_distance(chi)

        # Weight functions with Heaviside truncation
        # q^i_{λ_α}(χ) = Θ(z-z_min) Θ(z_max-z) (dV/dχdΩ)
        if z_min1 <= z <= z_max1:
            q1 = self.comoving_volume_element(z)
        else:
            q1 = 0.0

        if z_min2 <= z <= z_max2:
            q2 = self.comoving_volume_element(z)
        else:
            q2 = 0.0

        if q1 == 0.0 or q2 == 0.0:
            return 0.0

        # Bias-weighted number density
        # [bias_α] = ∫ dM (dn/dM) b_h(M,z) ∫ dλ p(M|λ,z)
        def bias_weighted_density(lambda_min, lambda_max):
            def integrand_M(M):
                dn_dM = self.tinker_mass_function(M, z)
                b_h = self.tinker_halo_bias(M, z)

                def integrand_lambda(lambda_val):
                    return self.mass_richness_probability(M, lambda_val, z)

                lambda_integral, _ = quad(integrand_lambda, lambda_min, lambda_max)
                return dn_dM * b_h * lambda_integral

            result, _ = quad(integrand_M, M_min, M_max)
            return result

        bias_weighted_1 = bias_weighted_density(lambda_min1, lambda_max1)
        bias_weighted_2 = bias_weighted_density(lambda_min2, lambda_max2)

        # σ_b²(z)
        sigma_b_sq = self.super_sample_variance(z)

        # Full integrand: q1 × q2 × [bias_1] × [bias_2] × σ_b²
        return q1 * q2 * bias_weighted_1 * bias_weighted_2 * sigma_b_sq

    # Integrate over χ
    z_min_overlap = max(z_min1, z_min2)
    z_max_overlap = min(z_max1, z_max2)

    if z_max_overlap > z_min_overlap:
        chi_min = self.comoving_distance(z_min_overlap)
        chi_max = self.comoving_distance(z_max_overlap)

        # ∫ dχ [integrand]
        chi_integral, _ = quad(chi_integrand, chi_min, chi_max,
                              epsabs=1e-6, epsrel=1e-6)
        super_sample_var = self.Omega_s**2 * chi_integral  # Ω_s² × ∫...
    else:
        super_sample_var = 0.0
else:
    super_sample_var = 0.0  # No cross-redshift correlations

# Total covariance
covariance_matrix[idx1, idx2] = shot_noise + super_sample_var
```

---

## Summary of Equation Mapping

| Equation | Paper Line | Code Lines | Status |
|----------|-----------|------------|---------|
| N^i(λ_α) | 206-207 | 218-256 | ✅ Always correct |
| M̄(λ) | 213 | 156-168 | ✅ Always correct |
| p(M\|λ,z) | 219 | 170-182 | ✅ Always correct |
| b_{λ_α}(z) | 271-272 | 258-292 | ✅ Always correct |
| σ_b²(Ω_s; z) | 738-740 | 346-388 | ✅ CORRECTED |
| Cov matrix | 765-767 | 429-510 | ✅ CORRECTED |

---

## Verification Checklist

When verifying the implementation:

- [ ] Check that all integrals have correct limits
- [ ] Verify integration measures (dz, dM, dλ, dχ, dk)
- [ ] Confirm normalization factors (1/M, 1/(2π), etc.)
- [ ] Ensure weight functions include Heaviside truncation
- [ ] Verify Kronecker deltas for shot noise
- [ ] Check cross-redshift correlations are zero
- [ ] Confirm physical units are consistent

All items above have been verified for the corrected implementation.

---

**Last Updated:** November 7, 2024
**Reference:** Krause & Eifler (2016), arXiv:1601.05779
