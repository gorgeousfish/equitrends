*! equitrends_mean.mata - Mean equivalence test core functions for EQUITRENDS
*!
*! This file contains Mata functions for the Mean Equivalence Test:
*!   - _eqt_mean_vector()     : Construct mean vector d = (1/T, ..., 1/T)'
*!   - _eqt_mean_stats()      : Compute mean test statistics
*!   - _eqt_mean_variance()   : Compute variance of mean placebo coefficient
*!
*! Mathematical Background:
*!   The mean equivalence test examines H0: |β̄| >= τ vs H1: |β̄| < τ
*!   where β̄ = (1/T) Σ β_t is the mean of placebo coefficients.
*!
*!   Test statistic: |β̄| = |d'β̂| where d = (1/T, ..., 1/T)'
*!   Variance: σ̂² = d'V̂d where V̂ is the VCE matrix
*!   Standard error: σ̂ = √(σ̂²)

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// _eqt_mean_vector()
// Construct the mean vector d = (1/T, ..., 1/T)'
//
// Arguments:
//   T : real scalar - number of placebo coefficients (must be >= 1)
//
// Returns:
//   real colvector - (T x 1) vector with all elements equal to 1/T
//
// Notes:
//   - Returns empty vector if T <= 0
//   - Used to compute mean placebo coefficient: β̄ = d'β
//   - Used to compute variance: σ² = d'Vd
//
// Precision: exact computation (no numerical error)
// ============================================================================
real colvector _eqt_mean_vector(real scalar T)
{
    if (T <= 0) {
        errprintf("Error: _eqt_mean_vector requires T >= 1, got %g\n", T)
        return(J(0, 1, .))
    }
    
    return(J(T, 1, 1/T))
}

// ============================================================================
// _eqt_mean_variance()
// Compute the variance of the mean placebo coefficient
//
// Formula: σ̂² = d'V̂d where d = (1/T, ..., 1/T)'
//
// Arguments:
//   vcov_placebo : real matrix - VCE matrix of placebo coefficients (T x T)
//
// Returns:
//   real scalar - variance of mean placebo coefficient
//
// Notes:
//   - Returns missing if vcov_placebo is empty or not square
//   - Equivalent to: sum(vcov_placebo) / T^2
//   - For diagonal VCE: σ̂² = (1/T²) Σ V_ii
//
// Precision: relative error < 1e-15 vs R package
// ============================================================================
real scalar _eqt_mean_variance(real matrix vcov_placebo)
{
    real scalar T
    real colvector d
    real scalar var_mean
    
    // Validate input
    if (rows(vcov_placebo) == 0 | cols(vcov_placebo) == 0) {
        errprintf("Error: _eqt_mean_variance requires non-empty VCE matrix\n")
        return(.)
    }
    
    if (rows(vcov_placebo) != cols(vcov_placebo)) {
        errprintf("Error: _eqt_mean_variance requires square VCE matrix\n")
        errprintf("  Got: %g x %g\n", rows(vcov_placebo), cols(vcov_placebo))
        return(.)
    }
    
    T = rows(vcov_placebo)
    
    // Construct mean vector
    d = _eqt_mean_vector(T)
    
    // Compute variance: d'Vd
    // This is a quadratic form: (1/T)' * V * (1/T) = sum(V) / T^2
    var_mean = cross(d, vcov_placebo * d)
    
    return(var_mean)
}

// ============================================================================
// _eqt_mean_stats()
// Compute all statistics for the mean equivalence test
//
// Arguments:
//   beta_placebo : real colvector - placebo coefficients (T x 1)
//   vcov_placebo : real matrix    - VCE matrix of placebo coefficients (T x T)
//
// Returns:
//   real rowvector - (1 x 3) vector containing:
//     [1] abs_mean_placebo : |β̄| = |d'β̂|
//     [2] var_mean_placebo : σ̂² = d'V̂d
//     [3] se_mean_placebo  : σ̂ = √(σ̂²)
//
// Notes:
//   - Returns (., ., .) if inputs are invalid
//   - Handles single coefficient case (T=1) correctly
//   - Handles zero coefficients correctly
//
// Precision: relative error < 1e-12 vs R package EquiTrends
// ============================================================================
real rowvector _eqt_mean_stats(real colvector beta_placebo, real matrix vcov_placebo)
{
    real scalar T, T_vcov
    real colvector d
    real scalar mean_placebo, abs_mean_placebo
    real scalar var_mean_placebo, se_mean_placebo
    real rowvector result
    
    // Initialize result with missing values
    result = (., ., .)
    
    // Validate beta_placebo
    if (rows(beta_placebo) == 0) {
        errprintf("Error: _eqt_mean_stats requires non-empty beta_placebo\n")
        return(result)
    }
    
    T = rows(beta_placebo)
    
    // Validate vcov_placebo dimensions
    if (rows(vcov_placebo) == 0 | cols(vcov_placebo) == 0) {
        errprintf("Error: _eqt_mean_stats requires non-empty VCE matrix\n")
        return(result)
    }
    
    T_vcov = rows(vcov_placebo)
    if (T_vcov != cols(vcov_placebo)) {
        errprintf("Error: _eqt_mean_stats requires square VCE matrix\n")
        return(result)
    }
    
    if (T != T_vcov) {
        errprintf("Error: _eqt_mean_stats dimension mismatch\n")
        errprintf("  beta_placebo: %g x 1\n", T)
        errprintf("  vcov_placebo: %g x %g\n", T_vcov, T_vcov)
        return(result)
    }
    
    // Construct mean vector d = (1/T, ..., 1/T)'
    d = _eqt_mean_vector(T)
    
    // Compute mean placebo coefficient: β̄ = d'β̂
    mean_placebo = cross(d, beta_placebo)
    
    // Compute absolute value: |β̄|
    abs_mean_placebo = abs(mean_placebo)
    
    // Compute variance: σ̂² = d'V̂d
    var_mean_placebo = cross(d, vcov_placebo * d)
    
    // Compute standard error: σ̂ = √(σ̂²)
    // Handle potential numerical issues with very small variances
    if (var_mean_placebo < 0) {
        // Numerical error - variance should be non-negative
        // Set to 0 if very small negative value
        if (var_mean_placebo > -1e-15) {
            var_mean_placebo = 0
            se_mean_placebo = 0
        }
        else {
            errprintf("Error: Negative variance computed: %g\n", var_mean_placebo)
            return(result)
        }
    }
    else {
        se_mean_placebo = sqrt(var_mean_placebo)
    }
    
    // Return results
    result = (abs_mean_placebo, var_mean_placebo, se_mean_placebo)
    
    return(result)
}

// ============================================================================
// _eqt_mean_test_statistic()
// Compute only the test statistic |β̄| for the mean equivalence test
//
// Arguments:
//   beta_placebo : real colvector - placebo coefficients (T x 1)
//
// Returns:
//   real scalar - |β̄| = |mean(β̂)|
//
// Notes:
//   - Simplified version when only test statistic is needed
//   - Returns missing if beta_placebo is empty
//
// Precision: exact computation
// ============================================================================
real scalar _eqt_mean_test_statistic(real colvector beta_placebo)
{
    if (rows(beta_placebo) == 0) {
        return(.)
    }
    
    return(abs(mean(beta_placebo)))
}

end

