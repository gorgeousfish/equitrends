*! equitrends_rms.mata - RMS Equivalence Test functions for EQUITRENDS
*!
*! This file contains Mata functions for the RMS Equivalence Test:
*!   - _eqt_sample_without_replacement() : Random sampling without replacement
*!   - _eqt_twfe_subsample()             : TWFE estimation on subsample
*!   - _eqt_rms_test()                   : Main RMS test calculation
*!   - _eqt_rms_ci()                     : Confidence interval for RMS
*!
*! Mathematical Framework:
*!   H0: β_RMS >= ζ  vs  H1: β_RMS < ζ
*!   where β_RMS = sqrt((1/T) Σ β_l²) is the root mean square of placebo coefficients
*!
*!   Self-normalized pivot statistic:
*!   M̂_n = (β̂_RMS²(1) - β_RMS²) / V̂_n
*!
*!   Rejection rule: β̂_RMS² < ζ² + Q_W(α) * V̂_n

version 16.0

mata:
mata set matastrict on

// ============================================================================
// _eqt_sample_without_replacement()
// Random sampling without replacement using Fisher-Yates partial shuffle
//
// Arguments:
//   N     : real scalar - population size (1 to N)
//   sub_N : real scalar - sample size to draw
//
// Returns:
//   real colvector - indices of sampled elements (1-indexed, length sub_N)
//
// Algorithm:
//   Fisher-Yates partial shuffle - only shuffle first sub_N elements
//   This is O(sub_N) instead of O(N) for full shuffle
//
// Notes:
//   - Uses Stata's runiform() for random number generation
//   - Sampling is by individuals, not observations
//   - No duplicate indices in output
// ============================================================================
real colvector _eqt_sample_without_replacement(real scalar N, real scalar sub_N)
{
    real colvector population, sample_idx
    real scalar i, j, temp
    
    // Validate inputs
    if (sub_N > N) {
        errprintf("Error: sample size (%g) cannot exceed population size (%g)\n", sub_N, N)
        _error(3498)
    }
    if (sub_N <= 0) {
        errprintf("Error: sample size must be positive, got %g\n", sub_N)
        _error(3498)
    }
    
    // Create population indices 1, 2, ..., N
    population = (1::N)
    
    // Fisher-Yates partial shuffle (only first sub_N elements)
    // For i = 1 to sub_N:
    //   Pick random j from [i, N]
    //   Swap population[i] and population[j]
    for (i = 1; i <= sub_N; i++) {
        // Random index from i to N (inclusive)
        // runiform() returns [0, 1), so floor(runiform() * (N - i + 1)) gives [0, N-i]
        j = floor(runiform(1, 1) * (N - i + 1)) + i
        
        // Swap population[i] and population[j]
        if (i != j) {
            temp = population[i]
            population[i] = population[j]
            population[j] = temp
        }
    }
    
    // Return first sub_N elements (the sampled indices)
    sample_idx = population[1::sub_N]
    
    return(sample_idx)
}


// ============================================================================
// _eqt_twfe_subsample()
// Estimate TWFE model on a subsample of individuals and return placebo coefficients
//
// Arguments:
//   Y              : real colvector - dependent variable (N_obs x 1)
//   X              : real matrix    - control variables (N_obs x p), can be empty
//   ID             : real colvector - individual identifiers (N_obs x 1)
//   G              : real colvector - treatment group indicator (N_obs x 1)
//   period         : real colvector - time period (N_obs x 1)
//   subset_indiv   : real colvector - indices of individuals to include
//   placebo_periods: real colvector - placebo period values
//   baseperiod     : real scalar    - base period value
//
// Returns:
//   real colvector - placebo coefficients (n_placebo x 1)
// ============================================================================
real colvector _eqt_twfe_subsample(real colvector Y, real matrix X,
                                    real colvector ID, real colvector G,
                                    real colvector period,
                                    real colvector subset_indiv,
                                    real colvector placebo_periods,
                                    real scalar baseperiod)
{
    real colvector keep_rows, Y_sub, ID_sub, G_sub, period_sub
    real matrix X_sub, placebo_vars, X_full, X_full_dm
    real colvector Y_dm, coefs, placebo_coefs
    real scalar n_obs, n_placebo, n_controls, i, n_sub, n_t_sub
    
    // ========================================================================
    // Step 1: Identify rows to keep (observations for individuals in subset)
    // ========================================================================
    
    n_obs = rows(Y)
    keep_rows = J(n_obs, 1, 0)
    
    for (i = 1; i <= n_obs; i++) {
        if (anyof(subset_indiv, ID[i])) {
            keep_rows[i] = 1
        }
    }
    
    // Convert to selection index
    keep_rows = selectindex(keep_rows)
    
    if (rows(keep_rows) == 0) {
        errprintf("Error: No observations found for subset individuals\n")
        _error(3498)
    }
    
    // ========================================================================
    // Step 2: Extract subsample data
    // ========================================================================
    
    Y_sub = Y[keep_rows]
    ID_sub = ID[keep_rows]
    G_sub = G[keep_rows]
    period_sub = period[keep_rows]
    
    n_controls = cols(X)
    if (n_controls > 0) {
        X_sub = X[keep_rows, .]
    }
    else {
        X_sub = J(rows(Y_sub), 0, .)
    }
    
    // ========================================================================
    // Step 3: Construct placebo dummy variables
    // placebo_t = (period == t) * G for each placebo period t
    // ========================================================================
    
    n_placebo = rows(placebo_periods)
    placebo_vars = J(rows(Y_sub), n_placebo, 0)
    
    for (i = 1; i <= n_placebo; i++) {
        // placebo_t = 1 if (period == placebo_periods[i]) AND (G == 1)
        placebo_vars[., i] = (period_sub :== placebo_periods[i]) :* G_sub
    }
    
    // Build full X matrix: [placebo_vars, X_controls]
    if (n_controls > 0) {
        X_full = placebo_vars, X_sub
    }
    else {
        X_full = placebo_vars
    }
    
    // ========================================================================
    // Step 4: Double demean (two-way fixed effects transformation)
    // Standard double demeaning is applied to remove individual and time
    // fixed effects: x_dm = x - mean_i - mean_t + mean
    // ========================================================================
    
    // Standard double demean: x_dm = x - mean_i - mean_t + mean
    Y_dm = _eqt_standard_double_demean_vec(Y_sub, ID_sub, period_sub)
    X_full_dm = _eqt_standard_double_demean(X_full, ID_sub, period_sub)
    
    // ========================================================================
    // Step 5: OLS estimation
    // ========================================================================
    
    // Solve (X'X)^{-1} X'Y using Cholesky decomposition
    coefs = _eqt_ols_cholesky(cross(X_full_dm, X_full_dm), cross(X_full_dm, Y_dm))
    
    // ========================================================================
    // Step 6: Extract placebo coefficients (first n_placebo elements)
    // ========================================================================
    
    placebo_coefs = coefs[1::n_placebo]
    
    return(placebo_coefs)
}


// ============================================================================
// _eqt_rms_test()
// Main RMS equivalence test calculation
//
// Arguments:
//   Y, X, ID, G, period, placebo_periods, baseperiod : Data inputs
//   threshold_specified : 1 if threshold is specified, 0 otherwise
//   threshold           : Equivalence threshold (if specified)
//   alpha               : Significance level for hypothesis test
//   no_lambda           : Number of subsamples for self-normalization
//   level               : Confidence level for CI (e.g., 95 for 95% CI)
//
// Returns (via st_numscalar and st_matrix):
//   r(MS_placebo), r(RMS_placebo), r(V_n), r(Q_W)
//   r(RMS_ci_lower), r(RMS_ci_upper), r(MS_ci_lower), r(MS_ci_upper), r(level)
//   r(MS_critical), r(RMS_critical), r(reject) [if threshold specified]
//   r(min_threshold) [if threshold not specified]
// ============================================================================
void _eqt_rms_test(real colvector Y, real matrix X,
                   real colvector ID, real colvector G,
                   real colvector period,
                   real colvector placebo_periods, real scalar baseperiod,
                   real scalar threshold_specified, real scalar threshold,
                   real scalar alpha, real scalar no_lambda, real scalar level)
{
    real colvector individuals, subset_index, subset_indiv
    real colvector placebo_coefs, diff_vec
    real matrix MS_lambda
    real scalar N, k, lambda_k, sub_N
    real scalar MS_full, RMS_full, V_n, Q_W
    real scalar MS_critical, RMS_critical, reject, min_threshold
    real scalar ci_alpha, q_lower, q_upper
    real scalar MS_ci_lower, MS_ci_upper, RMS_ci_lower, RMS_ci_upper
    
    // ========================================================================
    // Input validation: no_lambda must be at least 2
    // The V_n variance estimator requires computing differences between
    // subsample estimates, which necessitates at least two subsamples.
    // ========================================================================
    
    if (no_lambda < 2) {
        errprintf("Error: no_lambda must be at least 2 (required for V_n calculation)\n")
        errprintf("The V_n formula requires at least 2 subsamples to compute variance.\n")
        _error(3498)
    }
    
    // ========================================================================
    // Extract unique individual identifiers
    // Preserves first-appearance order to ensure consistent subsample selection
    // ========================================================================
    
    individuals = _equitrends_uniq_preserve_order(ID)
    N = rows(individuals)
    
    // ========================================================================
    // Data validation: ensure at least one individual is present
    // ========================================================================
    
    if (N == 0) {
        errprintf("Error: No valid individuals found in data.\n")
        errprintf("All ID values may be missing, or data is empty.\n")
        _error(3498)
    }
    
    // ========================================================================
    // Initialize MS_lambda storage
    // ========================================================================
    
    MS_lambda = J(no_lambda, 1, .)
    
    // ========================================================================
    // Subsample loop: estimate β̂(λ) for λ ∈ {1/K, 2/K, ..., 1}
    // ========================================================================
    
    for (k = 1; k <= no_lambda; k++) {
        lambda_k = k / no_lambda
        sub_N = floor(lambda_k * N)
        
        // Ensure subsample contains at least one individual
        if (sub_N < 1) {
            errprintf("Error: sub_N=0 at lambda_k=%g with N=%g individuals.\n", lambda_k, N)
            errprintf("N must be >= no_lambda. Please reduce no_lambda or use more individuals.\n")
            _error(3498)
        }
        
        // Random sample of individuals (without replacement)
        subset_index = _eqt_sample_without_replacement(N, sub_N)
        subset_indiv = individuals[subset_index]
        
        // Estimate TWFE on subsample and obtain placebo coefficients
        placebo_coefs = _eqt_twfe_subsample(Y, X, ID, G, period, 
                                            subset_indiv, placebo_periods, baseperiod)
        
        // Compute MS(λ) = (1/T) Σ β̂_l²(λ)
        MS_lambda[k] = mean(placebo_coefs :^ 2)
    }
    
    // The final iteration (k = no_lambda) corresponds to λ = 1 (full sample)
    // placebo_coefs now contains β̂(1), the full-sample placebo coefficients
    
    // ========================================================================
    // Calculate MS_full and RMS_full
    // ========================================================================
    
    MS_full = MS_lambda[no_lambda]
    RMS_full = sqrt(MS_full)
    
    // ========================================================================
    // Calculate self-normalized variance estimator V̂_n
    // V̂_n = sqrt( (1/(K-1)) Σ_{k=1}^{K-1} (MS(λ_k) - MS(1))² )
    // This is used for the pivotal test statistic M̂_n
    // ========================================================================
    
    diff_vec = MS_lambda[1::(no_lambda-1)] :- MS_full
    V_n = sqrt(mean(diff_vec :^ 2))
    
    // ========================================================================
    // Get W distribution critical value
    // ========================================================================
    
    Q_W = _equitrends_W_critical_value(alpha)
    
    // ========================================================================
    // Calculate test results
    // ========================================================================
    
    if (threshold_specified == 1) {
        // Critical value: ζ² + Q_W(α) * V̂_n
        MS_critical = threshold^2 + Q_W * V_n
        
        // Negative critical value indicates threshold is too small
        if (MS_critical < 0) {
            errprintf("The critical value is negative. Please enter a higher equivalence threshold.\n")
            _error(3498)
        }
        
        RMS_critical = sqrt(MS_critical)
        
        // Rejection rule: reject H0 if β̂_RMS² < ζ² + Q_W(α) * V̂_n
        reject = (MS_full < MS_critical) ? 1 : 0
        
        // Return results
        st_numscalar("r(MS_critical)", MS_critical)
        st_numscalar("r(RMS_critical)", RMS_critical)
        st_numscalar("r(reject)", reject)
    }
    else {
        // Minimum equivalence threshold: ζ* = sqrt(β̂_RMS² - Q_W(α) * V̂_n)
        // Since Q_W(α) < 0 for standard significance levels,
        // β̂_RMS² - Q_W(α) * V̂_n > 0 is guaranteed
        min_threshold = sqrt(MS_full - Q_W * V_n)
        
        st_numscalar("r(min_threshold)", min_threshold)
    }
    
    // ========================================================================
    // Return common results
    // ========================================================================
    
    st_numscalar("r(MS_placebo)", MS_full)
    st_numscalar("r(RMS_placebo)", RMS_full)
    st_numscalar("r(V_n)", V_n)
    st_numscalar("r(Q_W)", Q_W)
    
    // Return placebo coefficients as row vector
    st_matrix("r(placebo_coefs)", placebo_coefs')
    
    // Return MS_lambda as row vector
    st_matrix("r(MS_lambda)", MS_lambda')
    
    // ========================================================================
    // Confidence Interval Calculation
    // Based on Remark 4.2(b) of Dette & Schumann (2024):
    //   CI_α(β²_RMS) = [ β̂²_RMS + Q_W(α/2)·V̂_n,  β̂²_RMS + Q_W(1-α/2)·V̂_n ]
    //   CI_α(β_RMS)  = [ √max(0, lower_MS),  √upper_MS ]
    // ========================================================================
    
    // Convert confidence level to significance level
    ci_alpha = (100 - level) / 100
    
    // Get W distribution quantiles for CI bounds
    // Q_W(α/2) for lower bound, Q_W(1-α/2) for upper bound
    q_lower = _equitrends_W_critical_value(ci_alpha / 2)
    q_upper = _equitrends_W_critical_value(1 - ci_alpha / 2)
    
    // CI for β²_RMS (mean squared placebo coefficient)
    MS_ci_lower = MS_full + q_lower * V_n
    MS_ci_upper = MS_full + q_upper * V_n
    
    // CI for β_RMS (root mean squared placebo coefficient)
    // Ensure lower bound is non-negative
    RMS_ci_lower = sqrt(max((0, MS_ci_lower)))
    RMS_ci_upper = sqrt(MS_ci_upper)
    
    // Return CI results
    st_numscalar("r(RMS_ci_lower)", RMS_ci_lower)
    st_numscalar("r(RMS_ci_upper)", RMS_ci_upper)
    st_numscalar("r(MS_ci_lower)", MS_ci_lower)
    st_numscalar("r(MS_ci_upper)", MS_ci_upper)
    st_numscalar("r(level)", level)
}

end
