*! equitrends_meantest_core.mata - Core Mata functions for meanequivtest
*!
*! This file contains the core Mata functions for the meanequivtest command:
*!   - _meanequivtest_core_mata()  : Mean equivalence test core calculation
*!
*! These functions are called by _meanequivtest_core.ado
*!
*! Mathematical Framework:
*!   H0: |β̄| >= τ  vs  H1: |β̄| < τ
*!   where β̄ = (1/T) Σ β_t is the mean of placebo coefficients
*!
*!   Test statistic: |β̄| = |d'β̂| where d = (1/T, ..., 1/T)'
*!   Variance: σ̂² = d'V̂d
*!   Rejection rule: |β̄| < Q_{N_F(τ, σ̂²)}(α)

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// _meanequivtest_core_mata()
// Mean equivalence test core calculation
//
// Arguments:
//   y_name              : string scalar - outcome variable name
//   x_names             : string scalar - control variable names (space-separated)
//   id_name             : string scalar - individual ID variable name
//   g_name              : string scalar - treatment group variable name
//   period_name         : string scalar - time period variable name
//   no_placebos         : real scalar   - number of placebo coefficients
//   threshold_specified : real scalar   - 1 if threshold specified, 0 otherwise
//   threshold           : real scalar   - equivalence threshold
//   alpha               : real scalar   - significance level
//   vce_type            : string scalar - VCE type (ols/hc1/hc3/cr0)
//   cluster_name        : string scalar - cluster variable name
//   pretreat_periods    : real colvector - pretreatment periods
//   baseperiod          : real scalar   - base period
//
// Returns (via st_numscalar/st_matrix):
//   r(placebo_coefs)    - placebo coefficient matrix
//   r(vcov_placebo)     - placebo VCE matrix
//   r(abs_mean_placebo) - |β̄|
//   r(var_mean_placebo) - σ̂²
//   r(se_mean_placebo)  - σ̂
//   r(critical_value)   - critical value (when threshold specified)
//   r(p_value)          - p-value (when threshold specified)
//   r(reject)           - rejection decision (when threshold specified)
//   r(min_threshold)    - minimum threshold (when threshold not specified)
// ============================================================================
void _meanequivtest_core_mata(string scalar y_name, string scalar x_names,
                               string scalar id_name, string scalar g_name,
                               string scalar period_name,
                               real scalar no_placebos, real scalar threshold_specified,
                               real scalar threshold, real scalar alpha,
                               string scalar vce_type, string scalar cluster_name,
                               real colvector pretreat_periods, real scalar baseperiod)
{
    real colvector Y, ID, G, period, cluster_var
    real matrix X, X_dm, vcov_mat, vcov_placebo, placebo_vars
    real colvector Y_dm, unconstrained_coefs, placebo_coefs
    real rowvector mean_stats
    real scalar abs_mean_placebo, var_mean_placebo, se_mean_placebo
    real scalar critical_value, p_value, reject, min_threshold
    real scalar i, n_controls, n, n_t
    string rowvector x_varnames
    real matrix X_full, X_full_dm
    
    // Variables for multicollinearity detection
    real rowvector problematic_x, removed_placebo_idx, removed_control_idx
    real scalar effective_no_placebos, idx
    string scalar period_list
    
    // Variables for period label preservation
    real colvector effective_placebo_periods
    string colvector placebo_rownames
    
    // Variables for missing value handling
    real colvector valid_idx
    real scalar na_omitted
    
    // =========================================================================
    // Load data
    // =========================================================================
    
    Y = st_data(., y_name)
    ID = st_data(., id_name)
    G = st_data(., g_name)
    period = st_data(., period_name)
    
    // Load cluster variable if provided
    if (cluster_name != "") {
        cluster_var = st_data(., cluster_name)
    }
    else {
        cluster_var = ID  // Default to individual ID for clustering
    }
    
    // Load X variables if provided
    if (x_names != "") {
        x_varnames = tokens(x_names)
        X = st_data(., x_varnames)
        n_controls = cols(X)
    }
    else {
        X = J(rows(Y), 0, .)
        n_controls = 0
    }
    
    // =========================================================================
    // Handle missing values
    // =========================================================================
    valid_idx = _eqt_omit_na(Y, ID, G, period, X, cluster_var)
    na_omitted = rows(Y) - rows(valid_idx)
    if (na_omitted > 0) {
        printf("{txt}Warning: %g rows of pre-treatment data omitted due to NAs.\n", na_omitted)
        Y = Y[valid_idx]
        ID = ID[valid_idx]
        G = G[valid_idx]
        period = period[valid_idx]
        cluster_var = cluster_var[valid_idx]
        if (n_controls > 0) {
            X = X[valid_idx, .]
        }
    }
    
    // =========================================================================
    // Filter to pretreatment periods only
    // =========================================================================
    
    real colvector pretreat_idx, all_periods
    
    // Check if pretreat_periods is a sentinel value (1x1 matrix with -999999)
    real scalar has_pretreat
    has_pretreat = !(rows(pretreat_periods) == 1 && pretreat_periods[1] == -999999)
    
    if (has_pretreat) {
        // User specified pretreatment periods
        // Filter data to only include these periods
        pretreat_idx = J(rows(Y), 1, 0)
        for (i = 1; i <= rows(pretreat_periods); i++) {
            pretreat_idx = pretreat_idx :| (period :== pretreat_periods[i])
        }
        pretreat_idx = selectindex(pretreat_idx)
        
        Y = Y[pretreat_idx]
        ID = ID[pretreat_idx]
        G = G[pretreat_idx]
        period = period[pretreat_idx]
        cluster_var = cluster_var[pretreat_idx]
        if (n_controls > 0) {
            X = X[pretreat_idx, .]
        }
    }
    
    // =========================================================================
    // Validate non-empty data after filtering
    // =========================================================================
    if (rows(Y) == 0) {
        errprintf("Error: No observations remain after filtering to pretreatment periods.\n")
        errprintf("Possible causes:\n")
        errprintf("  - The specified pretreatment periods do not exist in the data\n")
        errprintf("  - All matching observations were removed due to missing values\n")
        _error(3498)
    }
    
    // =========================================================================
    // Determine placebo periods (pretreatment periods excluding base period)
    // =========================================================================
    
    real colvector placebo_periods
    real scalar base_period_val
    
    if (has_pretreat && baseperiod != -999999) {
        // User specified pretreatment periods and base period
        placebo_periods = select(pretreat_periods, pretreat_periods :!= baseperiod)
        base_period_val = baseperiod
    }
    else if (has_pretreat) {
        // User specified pretreatment periods but no base period
        base_period_val = max(pretreat_periods)
        placebo_periods = select(pretreat_periods, pretreat_periods :!= base_period_val)
    }
    else {
        // No pretreatment specified - use all periods
        real colvector unique_periods
        unique_periods = uniqrows(period)
        if (baseperiod != -999999) {
            base_period_val = baseperiod
        }
        else {
            base_period_val = max(unique_periods)
        }
        placebo_periods = select(unique_periods, unique_periods :!= base_period_val)
    }
    
    // =========================================================================
    // Create placebo indicators: placebo_t = (period == t) * G
    // =========================================================================
    
    placebo_vars = J(rows(Y), rows(placebo_periods), 0)
    for (i = 1; i <= rows(placebo_periods); i++) {
        // placebo_t = 1 if (period == placebo_periods[i]) AND (G == 1)
        placebo_vars[., i] = (period :== placebo_periods[i]) :* G
    }
    
    // Build full X matrix: [placebo_vars, X_controls]
    if (n_controls > 0) {
        X_full = placebo_vars, X
    }
    else {
        X_full = placebo_vars
    }
    
    // =========================================================================
    // Double demean (two-way fixed effects transformation)
    // =========================================================================
    
    Y_dm = _eqt_standard_double_demean_vec(Y, ID, period)
    X_full_dm = _eqt_standard_double_demean(X_full, ID, period)
    
    // Get panel dimensions
    n = rows(uniqrows(ID))
    n_t = rows(uniqrows(period))
    
    // =========================================================================
    // Detect and remove multicollinear variables
    // =========================================================================
    effective_no_placebos = no_placebos
    problematic_x = _eqt_detect_multicol_qr(X_full_dm)
    
    if (cols(problematic_x) > 0) {
        // Separate placebo and control removals
        removed_placebo_idx = J(1, 0, .)
        removed_control_idx = J(1, 0, .)
        
        for (i = 1; i <= cols(problematic_x); i++) {
            idx = problematic_x[i]
            if (idx <= no_placebos) {
                removed_placebo_idx = removed_placebo_idx, idx
            } else {
                removed_control_idx = removed_control_idx, idx
            }
        }
        
        // Display warning for removed placebos
        if (cols(removed_placebo_idx) > 0) {
            period_list = strofreal(placebo_periods[removed_placebo_idx[1]])
            for (i = 2; i <= cols(removed_placebo_idx); i++) {
                period_list = period_list + ", " + strofreal(placebo_periods[removed_placebo_idx[i]])
            }
            printf("{txt}Warning: The placebo corresponding to period(s) %s removed due to multicolinearity.\n", period_list)
        }
        
        // Display warning for removed controls
        if (cols(removed_control_idx) > 0) {
            real scalar ctrl_idx
            string scalar var_name_list
            ctrl_idx = removed_control_idx[1] - no_placebos
            if (ctrl_idx >= 1 & ctrl_idx <= cols(x_varnames)) {
                var_name_list = x_varnames[ctrl_idx]
            } else {
                var_name_list = strofreal(ctrl_idx)
            }
            for (i = 2; i <= cols(removed_control_idx); i++) {
                ctrl_idx = removed_control_idx[i] - no_placebos
                if (ctrl_idx >= 1 & ctrl_idx <= cols(x_varnames)) {
                    var_name_list = var_name_list + ", " + x_varnames[ctrl_idx]
                } else {
                    var_name_list = var_name_list + ", " + strofreal(ctrl_idx)
                }
            }
            printf("{txt}Warning: The following control variables were removed due to multicolinearity: %s\n", var_name_list)
        }
        
        // Remove problematic columns from demeaned X
        X_full_dm = _eqt_remove_multicol(X_full_dm, problematic_x)
        
        // Update effective number of placebos
        effective_no_placebos = no_placebos - cols(removed_placebo_idx)
        
        // Validate remaining placebo coefficients
        if (effective_no_placebos <= 0) {
            errprintf("Error: All placebo coefficients were removed due to multicollinearity.\n")
            errprintf("The Mean equivalence test cannot be performed without placebo coefficients.\n")
            errprintf("This indicates severe multicollinearity in your data.\n")
            _error(459)
        }
        
        // Update effective placebo periods
        if (cols(removed_placebo_idx) > 0) {
            real colvector keep_idx
            keep_idx = J(rows(placebo_periods), 1, 1)
            for (i = 1; i <= cols(removed_placebo_idx); i++) {
                keep_idx[removed_placebo_idx[i]] = 0
            }
            effective_placebo_periods = select(placebo_periods, keep_idx)
        }
        else {
            effective_placebo_periods = placebo_periods
        }
    }
    else {
        // No multicollinearity - use original placebo_periods
        effective_placebo_periods = placebo_periods
    }
    
    // =========================================================================
    // Compute OLS coefficients
    // =========================================================================
    
    unconstrained_coefs = _eqt_ols_cholesky(cross(X_full_dm, X_full_dm), cross(X_full_dm, Y_dm))
    
    // Extract placebo coefficients (using effective_no_placebos)
    placebo_coefs = unconstrained_coefs[1::effective_no_placebos]
    
    // =========================================================================
    // Compute variance-covariance matrix
    // 
    // Custom Mata VCE types:
    //   ols/homoskedastic     - Homoskedastic OLS VCE
    //   hc1/heteroskedastic   - HC1 White heteroskedasticity-robust VCE
    //   hac                   - HC3 Arellano panel-robust VCE
    //   cluster               - CR0 cluster-robust VCE (R package compatible)
    //   hc1_cluster           - HC1 cluster-robust VCE
    //
    // Native Stata VCE types (via regress command):
    //   hc2                   - HC2 leverage-adjusted VCE
    //   hc3                   - HC3 more conservative leverage-adjusted VCE
    //   cr0                   - CR0 cluster-robust without adjustment
    //   cr1                   - CR1 cluster-robust with DF adjustment (Stata default)
    // =========================================================================
    
    if (vce_type == "ols" || vce_type == "homoskedastic") {
        vcov_mat = _eqt_vcov_ols(X_full_dm, Y_dm, n, n_t)
    }
    else if (vce_type == "hc1_cluster") {
        // HC1 with clustering by individual (custom Mata)
        vcov_mat = _eqt_vcov_hc1_cluster(X_full_dm, Y_dm, unconstrained_coefs, cluster_var)
    }
    else if (vce_type == "hc1" || vce_type == "heteroskedastic" || vce_type == "hc" || vce_type == "robust") {
        // HC1 heteroskedasticity-robust standard errors (custom Mata)
        vcov_mat = _eqt_vcov_hc1(X_full_dm, Y_dm, unconstrained_coefs)
    }
    else if (vce_type == "hac") {
        // HC3 Arellano panel-robust VCE (custom Mata, for backward compatibility)
        vcov_mat = _eqt_vcov_hc3(X_full_dm, Y_dm, unconstrained_coefs, ID)
    }
    else if (vce_type == "cluster") {
        // CR0 cluster-robust VCE (custom Mata, R package compatible)
        vcov_mat = _eqt_vcov_cr0(X_full_dm, Y_dm, unconstrained_coefs, cluster_var)
    }
    // Native Stata VCE types
    else if (vce_type == "hc2") {
        // HC2 leverage-adjusted VCE (native Stata)
        vcov_mat = _eqt_vcov_hc2_native(X_full_dm, Y_dm)
    }
    else if (vce_type == "hc3") {
        // HC3 more conservative leverage-adjusted VCE (native Stata)
        vcov_mat = _eqt_vcov_hc3_native(X_full_dm, Y_dm)
    }
    else if (vce_type == "cr0") {
        // CR0 cluster-robust VCE without adjustment (native Stata)
        vcov_mat = _eqt_vcov_cr0_native(X_full_dm, Y_dm, cluster_var)
    }
    else if (vce_type == "cr1") {
        // CR1 cluster-robust VCE with DF adjustment (native Stata, Stata default)
        vcov_mat = _eqt_vcov_cr1_native(X_full_dm, Y_dm, cluster_var)
    }
    else {
        // Default to OLS
        vcov_mat = _eqt_vcov_ols(X_full_dm, Y_dm, n, n_t)
    }
    
    // Extract placebo VCE submatrix (using effective_no_placebos)
    vcov_placebo = _eqt_extract_vcov_placebo(vcov_mat, effective_no_placebos)
    
    // =========================================================================
    // Compute mean test statistics
    // =========================================================================
    
    mean_stats = _eqt_mean_stats(placebo_coefs, vcov_placebo)
    abs_mean_placebo = mean_stats[1]
    var_mean_placebo = mean_stats[2]
    se_mean_placebo = mean_stats[3]
    
    // =========================================================================
    // Perform hypothesis test or find minimum threshold
    // =========================================================================
    
    // Initialize return values
    critical_value = .
    p_value = .
    reject = .
    min_threshold = .
    
    if (threshold_specified == 1) {
        // ====================================================================
        // Specified threshold: compute critical value, p-value, and decision
        // ====================================================================
        
        // Critical value: Q_{N_F(τ, σ̂)}(α)
        critical_value = _eqt_qfoldnorm(alpha, threshold, se_mean_placebo)
        
        // P-value: P(|X| <= |β̄|) where X ~ N(τ, σ̂²)
        p_value = _eqt_pfoldnorm(abs_mean_placebo, threshold, se_mean_placebo)
        
        // Reject H0 if |β̄| < critical_value
        reject = (abs_mean_placebo < critical_value) ? 1 : 0
    }
    else {
        // ====================================================================
        // No threshold specified: find minimum threshold
        // ====================================================================
        
        // Find τ* such that pfoldnorm(|β̄|, τ*, σ̂) = α
        min_threshold = _eqt_mean_min_threshold(abs_mean_placebo, se_mean_placebo, alpha)
    }
    
    // =========================================================================
    // Return results to Stata
    // =========================================================================
    
    // Build row names using period labels
    placebo_rownames = J(effective_no_placebos, 1, "")
    for (i = 1; i <= effective_no_placebos; i++) {
        placebo_rownames[i] = "placebo_" + strofreal(effective_placebo_periods[i])
    }
    
    // Return placebo coefficients matrix
    st_matrix("r(placebo_coefs)", placebo_coefs)
    st_matrixrowstripe("r(placebo_coefs)", (J(effective_no_placebos, 1, ""), placebo_rownames))
    st_matrixcolstripe("r(placebo_coefs)", ("", "coef"))
    
    // Return VCE matrix
    st_matrix("r(vcov_placebo)", vcov_placebo)
    st_matrixrowstripe("r(vcov_placebo)", (J(effective_no_placebos, 1, ""), placebo_rownames))
    st_matrixcolstripe("r(vcov_placebo)", (J(effective_no_placebos, 1, ""), placebo_rownames))
    
    // Return placebo periods for display
    st_matrix("r(placebo_periods)", effective_placebo_periods)
    
    // Return scalars
    st_numscalar("r(abs_mean_placebo)", abs_mean_placebo)
    st_numscalar("r(var_mean_placebo)", var_mean_placebo)
    st_numscalar("r(se_mean_placebo)", se_mean_placebo)
    
    if (threshold_specified == 1) {
        st_numscalar("r(critical_value)", critical_value)
        st_numscalar("r(p_value)", p_value)
        st_numscalar("r(reject)", reject)
    }
    else {
        st_numscalar("r(min_threshold)", min_threshold)
    }
    
    // Return effective number of placebos (after multicollinearity removal)
    st_numscalar("r(no_placebos)", effective_no_placebos)
}

end
