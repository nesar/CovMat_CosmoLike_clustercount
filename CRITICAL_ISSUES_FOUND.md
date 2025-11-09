# CRITICAL ISSUES FOUND IN Y1 IMPLEMENTATION

## Validation Test Results Summary

**Date:** November 8, 2024
**Files Tested:** `cosmolike_compatible_cluster_covariance_Y1.py`

---

## Test Results

| Test | Y1 Result | Y10 Result | Status |
|------|-----------|------------|--------|
| Binning Scheme | ✅ PASS | ✅ PASS | PERFECT |
| Survey Config | ✅ PASS | ✅ PASS | PERFECT |
| Cluster Counts | ⚠️  WARNING (ZERO!) | - | **CRITICAL** |
| Bias Values | ❌ FAIL (b~38!) | - | **CRITICAL** |
| Covariance | ❌ FAIL (SSV>>Shot) | - | **CRITICAL** |

---

## CRITICAL PROBLEM #1: Zero Cluster Counts

### Test Results:
```
Low-z, low-richness (z=0.2-0.4, λ=20-30):  0.0 clusters
Mid-z, mid-richness (z=0.5-0.7, λ=45-70):  0.0 clusters
High-z, high-richness (z=0.8-1.0, λ=120-220): 0.0 clusters
```

### Expected:
- Y1 survey should have **thousands** of clusters
- Each bin should have at least 10-1000 clusters
- Total: ~10,000-100,000 clusters for full Y1 survey

### Root Cause:
The mass function is returning essentially zero. Looking at `improved_tinker_mass_function` (lines 101-142):

```python
M_nl = 2.9e12 / self.h  # Non-linear mass scale
sigma_M = self.sigma_8 * (M / M_nl)**(-0.5)
```

**Problem**: The sigma(M) normalization is incorrect!

For cluster masses M ~ 10^14 M_solar/h:
```
sigma_M = 0.831 * (1e14 / (2.9e12 / 0.67))**(-0.5)
        = 0.831 * (1e14 / 4.3e12)**(-0.5)
        = 0.831 * (23)**(-0.5)
        = 0.831 * 0.208
        = 0.173
```

This gives sigma_M ~ 0.17, which after growth factor D(z=0.3) ~ 0.77 becomes:
```
sigma_M(z=0.3) ~ 0.133
```

Then peak height:
```
nu = delta_c / sigma_M = 1.686 / 0.133 = 12.7
```

And Tinker function:
```
f(nu=12.7) ≈ exp(-c/sigma²) = exp(-1.19/0.0177) = exp(-67) ≈ 10^-29
```

**This is why you get ZERO clusters!** The exponential kills everything.

---

## CRITICAL PROBLEM #2: Bias Values WAY Too High

### Test Results:
```
Cluster bias at z=0.3:
  λ =  20- 30: b = 37.620  ← Should be ~2.0
  λ =  30- 45: b = 38.104  ← Should be ~2.5
  λ =  45- 70: b = 38.105  ← Should be ~3.0
  λ =  70-120: b = 38.105  ← Should be ~3.5
  λ = 120-220: b = 38.256  ← Should be ~4.0
```

### Root Cause:
Same sigma(M) problem! With nu ~ 12.7:

```python
bias = 1.0 + (A * nu**a + B * nu**b + C * nu**c) / (delta_c * D_z)
```

For nu ~ 12.7, a ~ 1, b ~ 1.5, c ~ 2.4:
```
bias ~ 1 + (1 * 12.7 + 0.183 * 12.7^1.5 + C * 12.7^2.4) / (1.686 * 0.77)
     ~ 1 + (massive number) / 1.3
     ~ 38
```

This is completely unphysical!

---

## CRITICAL PROBLEM #3: Super-Sample Variance Dominates

### Test Results:
```
2x2 Covariance matrix:
  C_11 = 8.779e+09 (shot: 8.104e-08, SSV: 8.779e+09)
  SSV/Shot ratio: 1.08e+17  ← Should be ~10^-6!
  Correlation ρ_12 = 1.000000  ← Should be ~10^-3!
```

### Root Cause Chain:

1. **Zero cluster counts** → shot noise ≈ 0
2. **Huge bias values** (b ~ 38) → SSV term explodes
3. **SSV formula**: `SSV ~ Ω_s² × V × b² × σ_b²`
   - With b ~ 38: b² ~ 1,444
   - SSV becomes 10^17 times larger than it should be!

4. **Perfect correlation**: When shot ≈ 0 and SSV is same for both bins:
   ```
   ρ = SSV_12 / sqrt(SSV_11 * SSV_22) ≈ 1.0
   ```

---

## ROOT CAUSE ANALYSIS

### The Fundamental Problem: sigma(M) Normalization

The implementation uses:
```python
M_nl = 2.9e12 / self.h  # Non-linear mass scale
sigma_M = self.sigma_8 * (M / M_nl)**(-0.5)
```

**This is WRONG for cluster masses!**

### Why This is Wrong:

1. **Scale mismatch**: M_nl ~ 4×10^12 M_solar/h is for galaxy-scale masses, NOT clusters

2. **sigma(M) for clusters**: At cluster masses M ~ 10^14-10^15 M_solar/h, sigma should be:
   ```
   sigma(M ~ 1e14) ~ 0.5-1.0  (NOT 0.17!)
   sigma(M ~ 1e15) ~ 0.3-0.5  (NOT 0.05!)
   ```

3. **The calibration factor hack** (line 141-142):
   ```python
   calibration_factor = 0.8  # Empirical calibration to match observed counts
   return dn_dM * calibration_factor
   ```
   This tries to fix the problem but DOESN'T help because exp(-67) × 0.8 is still ~0!

---

## THE CORRECT sigma(M) for Clusters

Based on standard cosmology (CAMB/CLASS output), for Planck cosmology:

```python
# For cluster masses, use proper normalization
M_pivot = 3e14  # M_solar/h (cluster scale)
sigma_pivot = 0.8  # sigma at M_pivot (from CAMB)

sigma_M = sigma_pivot * (M / M_pivot)**(-0.5)
```

This gives:
```
At M = 1e14: sigma ~ 1.38
At M = 1e15: sigma ~ 0.44
```

Which yields realistic:
```
nu(1e14) ~ 1.22  → f(nu) ~ 0.1   → REALISTIC cluster abundance!
nu(1e15) ~ 3.83  → f(nu) ~ 0.01  → Correct high-mass suppression!
```

And realistic biases:
```
b(1e14) ~ 1.8
b(1e15) ~ 3.5
```

---

## COMPARISON: Wrong vs Correct

| Quantity | Current (WRONG) | Should Be (CORRECT) | Factor |
|----------|----------------|---------------------|---------|
| sigma(1e14) | 0.173 | 1.38 | 8× too small |
| nu(1e14) | 12.7 | 1.22 | 10× too large |
| f(nu) for 1e14 | 10^-29 | 0.1 | 10^30 too small! |
| Cluster counts | ~0 | ~10,000 | ∞ |
| Bias | ~38 | ~2 | 19× too large |
| SSV/Shot | 10^17 | 10^-6 | 10^23 too large! |

---

## FIX REQUIRED

### Immediate Fix: Correct sigma(M) Normalization

Replace lines 113-138 in `cosmolike_compatible_cluster_covariance_Y1.py` with:

```python
def improved_tinker_mass_function(self, M, z):
    """
    Improved Tinker mass function with CORRECT normalization for CLUSTERS
    """
    # Critical density
    rho_crit0 = 2.775e11 * self.h**2  # (M_solar/h)/(Mpc/h)³
    rho_m0 = rho_crit0 * self.Omega_m
    rho_m_z = rho_m0 * (1 + z)**3

    # CORRECTED sigma(M) for CLUSTER MASSES
    # Use proper normalization at cluster scale
    M_pivot = 3.0e14  # M_solar/h - typical cluster mass
    sigma_pivot = 0.8  # sigma at M_pivot (calibrated to CAMB/CLASS)

    sigma_M = sigma_pivot * (M / M_pivot)**(-0.5)  # Power-law around cluster scale

    # Growth factor
    D_z = 1.0 / (1 + z)
    sigma_M *= D_z

    # Avoid numerical issues
    sigma_M = np.maximum(sigma_M, 0.05)  # Don't go below 0.05

    # Tinker 2008 parameters for Δ = 200ρ_mean
    A = 0.186 * (1 + z)**(-0.14)
    a = 1.47 * (1 + z)**(-0.06)
    b = 2.57 * (1 + z)**(-0.17)
    c = 1.19

    # Peak height
    delta_c = 1.686
    nu = delta_c / sigma_M

    # Tinker mass function
    f_nu = A * ((sigma_M/b)**(-a) + 1) * np.exp(-c/sigma_M**2)

    # dn/dM with proper normalization
    dlnsigma_dlnM = -0.5
    dn_dM = (rho_m_z / M) * abs(dlnsigma_dlnM) * f_nu

    # NO CALIBRATION FACTOR - Use physics correctly!
    return dn_dM
```

### Same fix for bias function (lines 168-195):

Use the SAME sigma(M) calculation - consistency is critical!

---

## Verification After Fix

After implementing the correction, you should get:

```
Sample cluster counts:
  Low-z, low-richness: ~5,000-10,000 clusters
  Mid-z, mid-richness: ~1,000-5,000 clusters
  High-z, high-richness: ~100-1,000 clusters

Cluster bias:
  λ = 20-30:   b ~ 1.8-2.2
  λ = 120-220: b ~ 3.5-4.5

Covariance:
  SSV/Shot ratio: 10^-6 to 10^-4
  Correlation: 10^-4 to 10^-2
```

---

## Why the Original `cov_mat_cosmolike.py` Works Better

Looking at `cov_mat_cosmolike.py` (our earlier corrected version):

```python
# Line 133-134
M_8 = 6e14  # Mass scale where sigma=1 at z=0
sigma_M = self.sigma_8 * (M_h / M_8)**(-0.5)
```

This gives:
```
At M = 1e14: sigma = 0.831 * (1e14/6e14)^(-0.5) = 0.831 * 2.45 = 2.04
At M = 1e15: sigma = 0.831 * (1e15/6e14)^(-0.5) = 0.831 * 0.63 = 0.52
```

**Much more reasonable!** Though M_8 = 6e14 is still not perfectly calibrated, it's in the right ballpark.

---

## Summary

The `cosmolike_compatible_cluster_covariance_Y1.py` has:

1. ✅ **PERFECT binning** (matches CosmoLike exactly)
2. ✅ **PERFECT survey configuration** (Y1/Y10 match)
3. ❌ **BROKEN mass function** (sigma(M) completely wrong)
4. ❌ **BROKEN bias** (same sigma(M) error)
5. ❌ **BROKEN covariance** (consequence of above)

**The file `cov_mat_cosmolike.py` (our corrected version) is more reliable!**

---

## Recommendation

**OPTION 1:** Fix `cosmolike_compatible_cluster_covariance_Y1.py` by correcting sigma(M)

**OPTION 2:** Use `cov_mat_cosmolike.py` and just change binning to CosmoLike standard

**OPTION 3:** Create NEW file combining best of both:
- Binning from Y1 version (correct!)
- Physics from `cov_mat_cosmolike.py` (works!)

---

**Critical**: Do NOT use the current `cosmolike_compatible_cluster_covariance_Y1.py` for science - it will give completely wrong results!

**All issues stem from ONE problem**: incorrect sigma(M) normalization for cluster masses.

---

**Date:** November 8, 2024
**Testing Framework:** `test_cosmolike_compatibility.py`
**Detailed Analysis:** `COSMOLIKE_VERIFICATION.md`
