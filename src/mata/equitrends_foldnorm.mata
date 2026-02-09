*! equitrends_foldnorm.mata - Folded Normal Distribution and IU Test Functions
*!
*! Implements the Intersection-Union (IU) equivalence test for pre-trend
*! analysis in difference-in-differences estimation.
*!
*! Core Components:
*!   - _eqt_pfoldnorm()               : Folded normal CDF
*!   - _eqt_qfoldnorm()               : Folded normal quantile function
*!   - _eqt_iu_critical_values()      : Critical values for IU test
*!   - _eqt_iu_test()                 : IU test decision rule
*!   - _eqt_iu_min_threshold_single() : Minimum threshold for single coefficient
*!   - _eqt_iu_min_threshold()        : Overall minimum equivalence threshold
*!   - _eqt_mean_min_threshold()      : Minimum threshold for Mean method
*!
*! Mathematical Background:
*!   If X ~ N(μ, σ²), then |X| follows a Folded Normal Distribution N_F(μ, σ²).
*!   
*!   CDF: F_{N_F}(x; μ, σ) = Φ((x-μ)/σ) + Φ((x+μ)/σ) - 1,  for x ≥ 0
*!   
*!   IU Test: Reject H₀: ||β||_∞ ≥ δ iff |β̂_t| < Q_{N_F(δ,se_t²)}(α) ∀t

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// INPUT VALIDATION FUNCTIONS
// ============================================================================

// ============================================================================
// _eqt_validate_foldnorm_params()
// Validate parameters for folded normal distribution functions
//
// Arguments:
//   sd              : real scalar - standard deviation (must be > 0)
//   p               : real scalar - probability for qfoldnorm (must be in (0,1))
//   alpha           : real scalar - significance level (must be in (0,1))
//   equiv_threshold : real scalar - equivalence threshold (must be > 0, per paper Section 3.1)
//   check_p         : real scalar - whether to check p parameter (1=yes, 0=no)
//   check_alpha     : real scalar - whether to check alpha parameter (1=yes, 0=no)
//   check_threshold : real scalar - whether to check equiv_threshold (1=yes, 0=no)
//
// Returns:
//   real scalar - 1 if all parameters valid, 0 otherwise (with error message)
// ============================================================================
real scalar _eqt_validate_foldnorm_params(real scalar sd,
                                          | real scalar p,
                                          real scalar alpha,
                                          real scalar equiv_threshold,
                                          real scalar check_p,
                                          real scalar check_alpha,
                                          real scalar check_threshold)
{
    // Set defaults for optional parameters
    if (args() < 5) check_p = 0
    if (args() < 6) check_alpha = 0
    if (args() < 7) check_threshold = 0
    
    // Validate sd > 0
    if (sd <= 0) {
        errprintf("Error: sd must be positive (got %g)\n", sd)
        return(0)
    }
    
    // Validate p in (0, 1) if requested
    if (check_p) {
        if (p <= 0 | p >= 1) {
            errprintf("Error: p must be between 0 and 1 (got %g)\n", p)
            return(0)
        }
    }
    
    // Validate alpha in (0, 1) if requested
    if (check_alpha) {
        if (alpha <= 0 | alpha >= 1) {
            errprintf("Error: alpha must be between 0 and 1 (got %g)\n", alpha)
            return(0)
        }
    }
    
    // Validate equiv_threshold > 0 if requested
    // δ > 0 is required; if δ = 0, H₀: ||β||_∞ ≥ 0 is always true
    if (check_threshold) {
        if (equiv_threshold <= 0) {
            errprintf("Error: equiv_threshold must be strictly positive (got %g)\n", equiv_threshold)
            return(0)
        }
    }
    
    return(1)
}

// ============================================================================
// FOLDED NORMAL DISTRIBUTION FUNCTIONS
// ============================================================================

// ============================================================================
// _eqt_pfoldnorm()
// Cumulative distribution function (CDF) of the Folded Normal Distribution
//
// F_{N_F}(x; μ, σ) = Φ((x-|μ|)/σ) + Φ((x+|μ|)/σ) - 1,  for x ≥ 0
// F_{N_F}(x; μ, σ) = 0,                                  for x < 0
//
// Special case (Half-Normal): F_{N_F}(x; 0, σ) = 2Φ(x/σ) - 1
//
// Arguments:
//   x    : real scalar - quantile value (support: x ≥ 0)
//   mean : real scalar - location parameter μ
//   sd   : real scalar - scale parameter σ (must be > 0)
//
// Returns:
//   real scalar - CDF value F(x; μ, σ) in [0, 1]
// ============================================================================
real scalar _eqt_pfoldnorm(real scalar x, real scalar mean, real scalar sd)
{
    real scalar z1, z2, cdf_val, mu_abs
    
    // Folded normal has support [0, ∞)
    if (x < 0) {
        return(0)
    }
    
    // Folded normal uses |μ| in the CDF formula
    mu_abs = abs(mean)
    
    // Edge case: when sd → 0, distribution degenerates to point mass at |μ|
    if (sd < 1e-15) {
        return(x >= mu_abs ? 1 : 0)
    }
    
    // Compute standardized values
    z1 = (x - mu_abs) / sd
    z2 = (x + mu_abs) / sd
    
    // CDF formula: Φ((x-|μ|)/σ) + Φ((x+|μ|)/σ) - 1
    cdf_val = normal(z1) + normal(z2) - 1
    
    return(cdf_val)
}

// ============================================================================
// _eqt_qfoldnorm()
// Quantile function (inverse CDF) of the Folded Normal Distribution
//
// Finds x such that F_{N_F}(x; μ, σ) = p via bisection method
//
// Arguments:
//   p    : real scalar - probability value in (0, 1)
//   mean : real scalar - location parameter μ
//   sd   : real scalar - scale parameter σ (must be > 0)
//
// Returns:
//   real scalar - quantile value x such that F(x; μ, σ) = p
// ============================================================================
real scalar _eqt_qfoldnorm(real scalar p, real scalar mean, real scalar sd)
{
    real scalar lower, upper, mid, f_lower, f_upper, f_mid
    real scalar tol, max_iter, iter
    real scalar c, mu_abs, result
    
    // Validate inputs
    if (p <= 0 | p >= 1) {
        errprintf("Error: p must be between 0 and 1 (got %g)\n", p)
        return(.)
    }
    if (sd <= 0) {
        errprintf("Error: sd must be positive (got %g)\n", sd)
        return(.)
    }
    
    mu_abs = abs(mean)
    
    // Edge case: when sd → 0, distribution degenerates to point mass at |μ|
    if (sd < 1e-15) {
        return(mu_abs)
    }
    
    // Edge case: when |μ|/σ → ∞, folded normal approaches N(|μ|, σ²)
    c = mu_abs / sd
    if (c > 1e10) {
        result = mu_abs + invnormal(p) * sd
        return(max((0, result)))
    }
    
    // Set search interval [0, mean + 10*sd]
    // Upper bound should be large enough to cover high quantiles
    lower = 0
    upper = mu_abs + 10 * sd
    
    // Extend upper bound if needed for very high probabilities
    while (_eqt_pfoldnorm(upper, mean, sd) < p) {
        upper = upper * 2
        if (upper > 1e10) {
            errprintf("Warning: Search interval may be insufficient\n")
            break
        }
    }
    
    // Bisection method for root finding
    // Find x such that pfoldnorm(x, mean, sd) - p = 0
    // Use tight tolerance for high precision
    tol = 1e-14
    max_iter = 100
    
    f_lower = _eqt_pfoldnorm(lower, mean, sd) - p
    f_upper = _eqt_pfoldnorm(upper, mean, sd) - p
    
    // Check if root is bracketed
    if (f_lower * f_upper > 0) {
        errprintf("Error: Root not bracketed in qfoldnorm\n")
        return(.)
    }
    
    for (iter = 1; iter <= max_iter; iter++) {
        mid = (lower + upper) / 2
        f_mid = _eqt_pfoldnorm(mid, mean, sd) - p
        
        // Check convergence based on function value
        if (abs(f_mid) < tol) {
            return(mid)
        }
        
        // Also check if interval is small enough relative to mid
        // This ensures relative precision
        if ((upper - lower) < 1e-14 * max((1, abs(mid)))) {
            return(mid)
        }
        
        // Update interval
        if (f_mid * f_lower < 0) {
            upper = mid
            f_upper = f_mid
        } else {
            lower = mid
            f_lower = f_mid
        }
    }
    
    // Return best estimate if max iterations reached
    return(mid)
}

// ============================================================================
// IU TEST CORE FUNCTIONS
// ============================================================================

// ============================================================================
// _eqt_iu_critical_values()
// Calculate critical values for the IU (Intersection-Union) test
//
// Formula: c_t(α) = Q_{N_F}(α; δ, se_t) for each t
//
// Arguments:
//   alpha           : real scalar - significance level
//   equiv_threshold : real scalar - equivalence threshold δ
//   beta_se         : real colvector - standard errors of placebo coefficients
//
// Returns:
//   real colvector - critical values (same length as beta_se)
// ============================================================================
real colvector _eqt_iu_critical_values(real scalar alpha, 
                                        real scalar equiv_threshold,
                                        real colvector beta_se)
{
    real colvector crit_values
    real scalar T, t
    
    T = rows(beta_se)
    crit_values = J(T, 1, .)
    
    // Calculate critical value for each coefficient
    for (t = 1; t <= T; t++) {
        crit_values[t] = _eqt_qfoldnorm(alpha, equiv_threshold, beta_se[t])
    }
    
    return(crit_values)
}

// ============================================================================
// _eqt_iu_test()
// IU (Intersection-Union) test decision rule
//
// Decision rule: Reject H₀ iff |β̂_t| < c_t(α) for ALL t ∈ {1,...,T}
//
// Uses strict inequality (<) per the equivalence testing framework
//
// Arguments:
//   abs_beta    : real colvector - absolute values of placebo coefficients
//   crit_values : real colvector - critical values (same length as abs_beta)
//
// Returns:
//   real scalar - 1 if reject H₀ (equivalence concluded), 0 otherwise
// ============================================================================
real scalar _eqt_iu_test(real colvector abs_beta, real colvector crit_values)
{
    real colvector reject_vec
    real scalar reject_H0
    
    // Strict inequality: |β̂_t| < c_t(α)
    reject_vec = abs_beta :< crit_values
    
    // IU principle: reject H₀ only if ALL coefficients satisfy the condition
    reject_H0 = (sum(reject_vec) == rows(reject_vec))
    
    return(reject_H0)
}

// ============================================================================
// _eqt_iu_min_threshold_single()
// Find minimum equivalence threshold for a single coefficient
//
// Finds δ* such that F_{N_F}(|β̂|; δ*, se) = α
// This is the smallest threshold at which equivalence can be concluded
//
// Arguments:
//   abs_coef : real scalar - absolute value of coefficient |β̂|
//   sd       : real scalar - standard error of coefficient
//   alpha    : real scalar - significance level
//
// Returns:
//   real scalar - minimum equivalence threshold δ*
// ============================================================================
real scalar _eqt_iu_min_threshold_single(real scalar abs_coef, 
                                          real scalar sd, 
                                          real scalar alpha)
{
    real scalar lower, upper, mid, f_lower, f_upper, f_mid
    real scalar tol, max_iter, iter
    real scalar target_p
    
    // Validate inputs
    // Note: sd=0 can occur legitimately when a coefficient is normalized to zero
    // (e.g., the reference period in DiD). In this case, silently return missing
    // so that the overall min_threshold calculation can proceed by ignoring this coefficient.
    if (sd <= 0) {
        return(.)
    }
    if (alpha <= 0 | alpha >= 1) {
        errprintf("Error: alpha must be between 0 and 1 (got %g)\n", alpha)
        return(.)
    }
    
    // Degenerate case: when |β̂| = 0, no finite δ satisfies the equation
    // Return upper boundary of search interval
    if (abs_coef == 0) {
        return(10 * sd)
    }
    
    // Search interval extends to accommodate large |β̂|/sd ratios
    lower = max((0, abs_coef - 10 * sd))
    upper = max((10 * sd, abs_coef + 10 * sd))
    
    // Solve f(δ) = pfoldnorm(abs_coef, δ, sd) - alpha = 0
    // Note: pfoldnorm is decreasing in the mean parameter δ
    
    target_p = alpha
    tol = 1e-14
    max_iter = 200
    
    f_lower = _eqt_pfoldnorm(abs_coef, lower, sd) - target_p
    f_upper = _eqt_pfoldnorm(abs_coef, upper, sd) - target_p
    
    // If root is not bracketed, return boundary value closest to target
    if (f_lower * f_upper > 0) {
        if (abs(f_lower) < abs(f_upper)) {
            return(lower)
        } else {
            return(upper)
        }
    }
    
    // Bisection method
    for (iter = 1; iter <= max_iter; iter++) {
        mid = (lower + upper) / 2
        f_mid = _eqt_pfoldnorm(abs_coef, mid, sd) - target_p
        
        // Check convergence based on function value
        if (abs(f_mid) < tol) {
            return(mid)
        }
        
        // Also check if interval is too small to continue
        if (upper - lower < 1e-14 * max((1, abs(mid)))) {
            return(mid)
        }
        
        // Update interval (note: f is decreasing in δ)
        if (f_mid * f_lower < 0) {
            upper = mid
            f_upper = f_mid
        } else {
            lower = mid
            f_lower = f_mid
        }
    }
    
    // Return best estimate if max iterations reached
    return(mid)
}

// ============================================================================
// _eqt_iu_min_threshold()
// Calculate overall minimum equivalence threshold for IU test
//
// Formula: δ* = max_{l ∈ {1,...,T}} δ_l*
// The overall threshold is the maximum of individual coefficient thresholds
//
// Arguments:
//   abs_beta       : real colvector - absolute values of placebo coefficients
//   beta_se        : real colvector - standard errors (same length as abs_beta)
//   alpha          : real scalar - significance level
//   min_thresholds : real colvector (optional output) - individual thresholds
//
// Returns:
//   real scalar - overall minimum equivalence threshold δ*
// ============================================================================
real scalar _eqt_iu_min_threshold(real colvector abs_beta,
                                   real colvector beta_se,
                                   real scalar alpha,
                                   | real colvector min_thresholds)
{
    real scalar T, l, min_threshold
    real colvector valid_thresholds
    real colvector valid_indices
    
    // Validate inputs
    T = rows(abs_beta)
    if (T != rows(beta_se)) {
        errprintf("Error: abs_beta and beta_se must have the same length\n")
        return(.)
    }
    if (T == 0) {
        errprintf("Error: beta_se cannot be empty\n")
        return(.)
    }
    
    // Initialize output vector for individual thresholds
    min_thresholds = J(T, 1, .)
    
    // Calculate minimum threshold for each coefficient
    for (l = 1; l <= T; l++) {
        min_thresholds[l] = _eqt_iu_min_threshold_single(abs_beta[l], beta_se[l], alpha)
    }
    
    // Filter out missing values (coefficients with sd <= 0 are excluded)
    valid_indices = selectindex(min_thresholds :< .)
    if (length(valid_indices) == 0) {
        errprintf("Error: All placebo coefficients have sd<=0, cannot calculate min_threshold\n")
        return(.)
    }
    valid_thresholds = min_thresholds[valid_indices]
    
    // Overall minimum threshold is the maximum of valid individual thresholds
    min_threshold = max(valid_thresholds)
    
    return(min_threshold)
}

// ============================================================================
// MEAN TEST MINIMUM THRESHOLD FUNCTION
// ============================================================================

// ============================================================================
// _eqt_mean_min_threshold()
// Find minimum equivalence threshold for Mean method
//
// Finds τ* such that F_{N_F}(|β̄|; τ*, sd) = α
// This is the smallest threshold at which equivalence of the average
// placebo coefficient can be concluded
//
// Arguments:
//   abs_mean_coef : real scalar - absolute value of mean coefficient |β̄|
//   sd            : real scalar - standard error of the mean coefficient
//   alpha         : real scalar - significance level
//
// Returns:
//   real scalar - minimum equivalence threshold τ*
// ============================================================================
real scalar _eqt_mean_min_threshold(real scalar abs_mean_coef, 
                                     real scalar sd, 
                                     real scalar alpha)
{
    real scalar lower, upper, mid, f_lower, f_upper, f_mid
    real scalar tol, max_iter, iter
    real scalar target_p
    
    // Validate inputs
    if (sd <= 0) {
        errprintf("Error: sd must be positive (got %g)\n", sd)
        return(.)
    }
    if (alpha <= 0 | alpha >= 1) {
        errprintf("Error: alpha must be between 0 and 1 (got %g)\n", alpha)
        return(.)
    }
    if (abs_mean_coef < 0) {
        errprintf("Error: abs_mean_coef must be non-negative (got %g)\n", abs_mean_coef)
        return(.)
    }
    
    // Degenerate case: when |β̄| = 0, return upper boundary of search interval
    if (abs_mean_coef == 0) {
        return(4 * sd)
    }
    
    // Search interval for bisection
    lower = max((0, abs_mean_coef - 4 * sd))
    upper = abs_mean_coef + 4 * sd
    
    // Solve f(τ) = pfoldnorm(abs_mean_coef, τ, sd) - alpha = 0
    // Note: pfoldnorm is decreasing in the mean parameter τ
    
    target_p = alpha
    tol = 1e-14
    max_iter = 200
    
    f_lower = _eqt_pfoldnorm(abs_mean_coef, lower, sd) - target_p
    f_upper = _eqt_pfoldnorm(abs_mean_coef, upper, sd) - target_p
    
    // If root is not bracketed, return boundary value closest to target
    if (f_lower * f_upper > 0) {
        if (abs(f_lower) < abs(f_upper)) {
            return(lower)
        } else {
            return(upper)
        }
    }
    
    // Bisection method
    for (iter = 1; iter <= max_iter; iter++) {
        mid = (lower + upper) / 2
        f_mid = _eqt_pfoldnorm(abs_mean_coef, mid, sd) - target_p
        
        // Check convergence based on function value
        if (abs(f_mid) < tol) {
            return(mid)
        }
        
        // Also check if interval is too small to continue
        if (upper - lower < 1e-14 * max((1, abs(mid)))) {
            return(mid)
        }
        
        // Update interval (note: f is decreasing in τ)
        if (f_mid * f_lower < 0) {
            upper = mid
            f_upper = f_mid
        } else {
            lower = mid
            f_lower = f_mid
        }
    }
    
    // Return best estimate if max iterations reached
    return(mid)
}

// ============================================================================
// ADO INTERFACE HELPER FUNCTIONS
// ============================================================================

// ============================================================================
// _eqt_ado_iu_with_threshold()
// Helper function for ado interface - IU test with specified threshold
//
// Arguments:
//   abs_beta_name : string scalar - name of Stata matrix with |β|
//   se_name       : string scalar - name of Stata matrix with SE
//   alpha         : real scalar - significance level
//   threshold     : real scalar - equivalence threshold δ
//
// Returns: (via Stata r() scalars and matrices)
//   r(reject)          - 1 if reject H0, 0 otherwise
//   r(equiv_threshold) - threshold used
//   r(alpha)           - alpha used
//   r(n_coefs)         - number of coefficients
//   r(critical_values) - matrix of critical values
// ============================================================================
void _eqt_ado_iu_with_threshold(string scalar abs_beta_name, 
                                 string scalar se_name,
                                 real scalar alpha,
                                 real scalar threshold)
{
    real colvector abs_beta, se, crit_values
    real scalar reject, n_coefs
    
    // Get matrices from Stata
    abs_beta = st_matrix(abs_beta_name)
    se = st_matrix(se_name)
    
    // Ensure column vectors
    if (cols(abs_beta) > rows(abs_beta)) abs_beta = abs_beta'
    if (cols(se) > rows(se)) se = se'
    
    n_coefs = rows(abs_beta)
    
    // Calculate critical values
    crit_values = _eqt_iu_critical_values(alpha, threshold, se)
    
    // Perform IU test
    reject = _eqt_iu_test(abs_beta, crit_values)
    
    // Return results to Stata
    st_numscalar("r(reject)", reject)
    st_numscalar("r(equiv_threshold)", threshold)
    st_numscalar("r(alpha)", alpha)
    st_numscalar("r(n_coefs)", n_coefs)
    st_matrix("r(critical_values)", crit_values)
    st_matrixrowstripe("r(critical_values)", (J(n_coefs, 1, ""), strofreal(1::n_coefs)))
    st_matrixcolstripe("r(critical_values)", ("", "crit_value"))
}

// ============================================================================
// _eqt_ado_iu_min_threshold()
// Helper function for ado interface - calculate minimum threshold
//
// Arguments:
//   abs_beta_name : string scalar - name of Stata matrix with |β|
//   se_name       : string scalar - name of Stata matrix with SE
//   alpha         : real scalar - significance level
//
// Returns: (via Stata r() scalars and matrices)
//   r(min_equiv_threshold) - overall minimum threshold
//   r(alpha)               - alpha used
//   r(n_coefs)             - number of coefficients
//   r(min_thresholds)      - matrix of individual minimum thresholds
// ============================================================================
void _eqt_ado_iu_min_threshold(string scalar abs_beta_name,
                                string scalar se_name,
                                real scalar alpha)
{
    real colvector abs_beta, se, min_thresholds
    real scalar min_threshold, n_coefs
    
    // Get matrices from Stata
    abs_beta = st_matrix(abs_beta_name)
    se = st_matrix(se_name)
    
    // Ensure column vectors
    if (cols(abs_beta) > rows(abs_beta)) abs_beta = abs_beta'
    if (cols(se) > rows(se)) se = se'
    
    n_coefs = rows(abs_beta)
    
    // Calculate minimum threshold
    min_thresholds = J(n_coefs, 1, .)
    min_threshold = _eqt_iu_min_threshold(abs_beta, se, alpha, min_thresholds)
    
    // Return results to Stata
    st_numscalar("r(min_equiv_threshold)", min_threshold)
    st_numscalar("r(alpha)", alpha)
    st_numscalar("r(n_coefs)", n_coefs)
    st_matrix("r(min_thresholds)", min_thresholds)
    st_matrixrowstripe("r(min_thresholds)", (J(n_coefs, 1, ""), strofreal(1::n_coefs)))
    st_matrixcolstripe("r(min_thresholds)", ("", "min_threshold"))
}

end
