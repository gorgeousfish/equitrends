*! equitrends_maxtest_core.mata - Core Mata functions for max equivalence test
*!
*! Implements the intersection-union (IU) and bootstrap tests for the hypothesis:
*!   H0: ||beta||_inf >= delta vs H1: ||beta||_inf < delta
*! where beta is the vector of placebo coefficients from a two-way fixed effects model.
*!
*! Main functions:
*!   - _maxequivtest_core_iu_mata()      : IU method core calculation
*!   - _maxequivtest_core_bootstrap()    : Bootstrap method core calculation

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// _maxequivtest_core_iu_mata()
// IU method core calculation for maxequivtest
//
// Arguments:
//   y_name           : string scalar - outcome variable name
//   x_names          : string scalar - control variable names (space-separated)
//   id_name          : string scalar - individual ID variable name
//   g_name           : string scalar - treatment group variable name
//   period_name      : string scalar - time period variable name
//   no_placebos      : real scalar   - number of placebo coefficients
//   delta_specified  : real scalar   - 1 if threshold specified, 0 otherwise
//   delta            : real scalar   - equivalence threshold
//   alpha            : real scalar   - significance level
//   vce_type         : string scalar - VCE type
//   cluster_name     : string scalar - cluster variable name
//   pretreat_periods : real colvector - pretreatment periods
//   baseperiod       : real scalar   - base period
//
// Returns (via st_numscalar/st_matrix):
//   r(placebo_coefs)   - placebo coefficient matrix
//   r(max_abs_coef)    - maximum absolute placebo coefficient
//   r(critical_value)  - critical value (when delta specified)
//   r(reject)          - rejection decision (when delta specified)
//   r(min_delta)       - minimum threshold (when delta not specified)
//   r(sigma_hathat_c)  - constrained variance estimate
// ============================================================================
void _maxequivtest_core_iu_mata(string scalar y_name, string scalar x_names,
                                 string scalar id_name, string scalar g_name,
                                 string scalar period_name,
                                 real scalar no_placebos, real scalar delta_specified,
                                 real scalar delta, real scalar alpha,
                                 string scalar vce_type, string scalar cluster_name,
                                 real colvector pretreat_periods, real scalar baseperiod)
{
    real colvector Y, ID, G, period, cluster_var
    real matrix X, X_dm, XtX_inv, vcov_mat, placebo_vars
    real colvector Y_dm, unconstrained_coefs, placebo_coefs, placebo_se
    real colvector abs_beta, crit_values
    real scalar max_abs_coef, sigma_hathat_c, critical_value, reject, min_delta
    real scalar i, n_controls, has_pretreat, n, n_t
    string rowvector x_varnames
    real matrix X_full, X_full_dm
    real colvector pretreat_idx, placebo_periods
    real scalar base_period_val
    
    // Variables for multicollinearity detection
    real rowvector problematic_x, removed_placebo_idx, removed_control_idx
    real scalar effective_no_placebos, idx
    string scalar period_list
    
    // Variables for preserving real period labels
    real colvector effective_placebo_periods
    string colvector placebo_rownames
    
    // Variables for NA handling
    real colvector valid_idx
    real scalar na_omitted
    
    // Load FULL data (both treatment and control groups)
    Y = st_data(., y_name)
    ID = st_data(., id_name)
    G = st_data(., g_name)
    period = st_data(., period_name)
    
    // Load cluster variable if provided
    // When cluster variable not specified, default to individual ID
    if (cluster_name != "") {
        cluster_var = st_data(., cluster_name)
    }
    else {
        // Default to ID for cluster VCE
        cluster_var = ID
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
    
    // ========================================================================
    // Handle missing values
    // Remove all rows with NA in Y, ID, G, period, X, or cluster before processing
    // ========================================================================
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
    
    // ========================================================================
    // Filter to pretreatment periods
    // Data is filtered prior to estimation to ensure only pretreatment
    // observations are used when computing placebo coefficients
    // ========================================================================
    
    // Check if pretreat_periods is a sentinel value (1x1 matrix with -999999)
    has_pretreat = !(rows(pretreat_periods) == 1 && pretreat_periods[1] == -999999)
    
    if (has_pretreat) {
        // Filter data to user-specified pretreatment periods
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
    
    // ========================================================================
    // Check for empty data after pretreatment filtering
    // ========================================================================
    if (rows(Y) == 0) {
        errprintf("Error: No observations remain after filtering to pretreatment periods.\n")
        errprintf("Possible causes:\n")
        errprintf("  - The specified pretreatment periods do not exist in the data\n")
        errprintf("  - All matching observations were removed due to missing values\n")
        _error(3498)
    }
    
    // Determine placebo periods (pretreatment periods excluding base period)
    if (has_pretreat && baseperiod != -999999) {
        placebo_periods = select(pretreat_periods, pretreat_periods :!= baseperiod)
        base_period_val = baseperiod
    }
    else if (has_pretreat) {
        base_period_val = max(pretreat_periods)
        placebo_periods = select(pretreat_periods, pretreat_periods :!= base_period_val)
    }
    else {
        // No pretreatment specified - use all periods
        real colvector unique_periods
        unique_periods = uniqrows(period)
        // Check if user specified baseperiod; if not, use maximum period
        if (baseperiod != -999999) {
            base_period_val = baseperiod
        }
        else {
            base_period_val = max(unique_periods)
        }
        placebo_periods = select(unique_periods, unique_periods :!= base_period_val)
    }
    
    // Create placebo indicators: placebo_t = G * I(period == t)
    // This is the interaction between treatment group and time period
    placebo_vars = J(rows(Y), rows(placebo_periods), 0)
    for (i = 1; i <= rows(placebo_periods); i++) {
        placebo_vars[., i] = G :* (period :== placebo_periods[i])
    }
    
    // Build full X matrix: [placebo_vars, X_controls]
    if (n_controls > 0) {
        X_full = placebo_vars, X
    }
    else {
        X_full = placebo_vars
    }
    
    // ========================================================================
    // Standard two-way fixed effects demeaning
    // Apply within-transformation: x_dm = x - mean_i - mean_t + mean
    // ========================================================================
    Y_dm = _eqt_standard_double_demean_vec(Y, ID, period)
    X_full_dm = _eqt_standard_double_demean(X_full, ID, period)
    
    // Get panel dimensions
    n = rows(uniqrows(ID))
    n_t = rows(uniqrows(period))
    
    // ========================================================================
    // Check and remove multicollinearity before OLS and VCE calculation
    // This ensures VCE calculation with invsym(X'X) does not fail
    // ========================================================================
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
        
        // Update effective_placebo_periods by removing multicollinear periods
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
    
    // Check if all placebos were removed due to multicollinearity
    if (effective_no_placebos <= 0) {
        errprintf("Error: All placebo coefficients were removed due to multicollinearity.\n")
        errprintf("The IU equivalence test cannot be performed without placebo coefficients.\n")
        errprintf("This indicates severe multicollinearity in your data.\n")
        _error(459)
    }
    
    // Compute unconstrained OLS (using reduced X if multicollinearity was found)
    unconstrained_coefs = _eqt_ols_cholesky(cross(X_full_dm, X_full_dm), cross(X_full_dm, Y_dm))
    
    // Extract placebo coefficients (using effective_no_placebos)
    placebo_coefs = unconstrained_coefs[1::effective_no_placebos]
    
    // Compute max absolute coefficient
    max_abs_coef = max(abs(placebo_coefs))
    
    // ========================================================================
    // Compute VCE matrix based on vce_type
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
    // ========================================================================
    
    if (vce_type == "ols" || vce_type == "homoskedastic") {
        vcov_mat = _eqt_vcov_ols(X_full_dm, Y_dm, n, n_t)
    }
    else if (vce_type == "hc1_cluster") {
        // HC1 cluster-robust VCE (custom Mata)
        vcov_mat = _eqt_vcov_hc1_cluster(X_full_dm, Y_dm, unconstrained_coefs, cluster_var)
    }
    else if (vce_type == "hc1" || vce_type == "heteroskedastic" || vce_type == "hc" || vce_type == "robust") {
        // HC1 White heteroskedasticity-robust VCE (custom Mata)
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
    
    // Extract placebo standard errors from VCE matrix (using effective_no_placebos)
    placebo_se = sqrt(diagonal(vcov_mat[1::effective_no_placebos, 1::effective_no_placebos]))
    
    // Compute sigma_hathat_c (residual variance) - for display purposes
    sigma_hathat_c = _eqt_sigma_hathat_c(unconstrained_coefs, X_full_dm, Y_dm, ID, period)
    
    // Initialize return values
    critical_value = .
    reject = .
    min_delta = .
    
    if (delta_specified == 1) {
        // Compute critical values using IU method
        abs_beta = abs(placebo_coefs)
        
        // Use _eqt_iu_critical_values (plural) which takes all SEs at once
        crit_values = _eqt_iu_critical_values(alpha, delta, placebo_se)
        
        // Perform IU test
        reject = _eqt_iu_test(abs_beta, crit_values)
        
        // Return max critical value for display
        critical_value = max(crit_values)
    }
    else {
        // Compute minimum threshold using IU method (using effective_no_placebos)
        real colvector min_thresholds
        min_thresholds = J(effective_no_placebos, 1, .)
        min_delta = _eqt_iu_min_threshold(abs(placebo_coefs), placebo_se, alpha, min_thresholds)
    }
    
    // Build row names using real period labels (e.g., "placebo_1", "placebo_2", "placebo_5")
    placebo_rownames = J(effective_no_placebos, 1, "")
    for (i = 1; i <= effective_no_placebos; i++) {
        placebo_rownames[i] = "placebo_" + strofreal(effective_placebo_periods[i])
    }
    
    // Return results to Stata (using effective_no_placebos after multicollinearity removal)
    st_matrix("r(placebo_coefs)", placebo_coefs)
    st_matrixrowstripe("r(placebo_coefs)", (J(effective_no_placebos, 1, ""), placebo_rownames))
    st_matrixcolstripe("r(placebo_coefs)", ("", "coef"))
    
    // Return standard errors
    st_matrix("r(placebo_se)", placebo_se)
    st_matrixrowstripe("r(placebo_se)", (J(effective_no_placebos, 1, ""), placebo_rownames))
    st_matrixcolstripe("r(placebo_se)", ("", "se"))
    
    // Return variance-covariance matrix
    st_matrix("r(vcov_placebo)", vcov_mat[1::effective_no_placebos, 1::effective_no_placebos])
    
    st_numscalar("r(max_abs_coef)", max_abs_coef)
    st_numscalar("r(critical_value)", critical_value)
    st_numscalar("r(reject)", reject)
    st_numscalar("r(min_delta)", min_delta)
    st_numscalar("r(sigma_hathat_c)", sigma_hathat_c)
    
    // Return critical values or min thresholds for each coefficient
    if (delta_specified == 1) {
        st_matrix("r(crit_values)", crit_values)
        st_matrixrowstripe("r(crit_values)", (J(effective_no_placebos, 1, ""), placebo_rownames))
        st_matrixcolstripe("r(crit_values)", ("", "cv"))
    }
    else {
        st_matrix("r(min_thresholds)", min_thresholds)
        st_matrixrowstripe("r(min_thresholds)", (J(effective_no_placebos, 1, ""), placebo_rownames))
        st_matrixcolstripe("r(min_thresholds)", ("", "mt"))
    }
    
    // Return placebo_periods for display use
    st_matrix("r(placebo_periods)", effective_placebo_periods)
    
    // Return effective_no_placebos (actual number of placebos after multicollinearity removal)
    st_numscalar("r(no_placebos)", effective_no_placebos)
}


// ============================================================================
// _maxequivtest_core_bootstrap()
// Bootstrap method core calculation for maxequivtest
//
// Arguments:
//   y_name           : string scalar - outcome variable name
//   x_names          : string scalar - control variable names (space-separated)
//   id_name          : string scalar - individual ID variable name
//   g_name           : string scalar - treatment group variable name
//   period_name      : string scalar - time period variable name
//   no_placebos      : real scalar   - number of placebo coefficients
//   delta_specified  : real scalar   - 1 if threshold specified, 0 otherwise
//   delta            : real scalar   - equivalence threshold
//   alpha            : real scalar   - significance level
//   B                : real scalar   - number of bootstrap replications
//   type             : string scalar - "Boot" or "Wild"
//   seed             : real scalar   - random seed (0 = no seed)
//   pretreat_periods : real colvector - pretreatment periods
//   baseperiod       : real scalar   - base period
//
// Returns (via st_numscalar/st_matrix):
//   r(placebo_coefs)   - placebo coefficient matrix
//   r(max_abs_coef)    - maximum absolute placebo coefficient
//   r(critical_value)  - critical value (when delta specified)
//   r(reject)          - rejection decision (when delta specified)
//   r(min_delta)       - minimum threshold (when delta not specified)
//   r(sigma_hathat_c)  - constrained variance estimate
// ============================================================================
void _maxequivtest_core_bootstrap(string scalar y_name, string scalar x_names,
                                   string scalar id_name, string scalar g_name,
                                   string scalar period_name,
                                   real scalar no_placebos, real scalar delta_specified,
                                   real scalar delta, real scalar alpha, real scalar B,
                                   string scalar type, real scalar seed,
                                   real colvector pretreat_periods, real scalar baseperiod)
{
    real colvector Y, ID, G, period
    real matrix X, WD, X_dm, placebo_vars
    real colvector Y_dm, unconstrained_coefs, placebo_coefs
    real scalar max_abs_coef, sigma_hathat_c, critical_value, reject, min_delta
    real scalar i, n_controls, has_pretreat
    string rowvector x_varnames
    real matrix X_full, X_full_dm
    struct _eqt_bootstrap_result scalar boot_result
    struct _eqt_min_delta_result scalar min_delta_result
    real colvector pretreat_idx, placebo_periods
    real scalar base_period_val
    
    // Variables for NA handling
    real colvector valid_idx
    real scalar na_omitted
    
    // Load FULL data (both treatment and control groups)
    Y = st_data(., y_name)
    ID = st_data(., id_name)
    G = st_data(., g_name)
    period = st_data(., period_name)
    
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
    
    // ========================================================================
    // Handle missing values
    // Remove all rows with NA in Y, ID, G, period, or X before processing
    // ========================================================================
    valid_idx = _eqt_omit_na(Y, ID, G, period, X)
    na_omitted = rows(Y) - rows(valid_idx)
    if (na_omitted > 0) {
        printf("{txt}Warning: %g rows of pre-treatment data omitted due to NAs.\n", na_omitted)
        Y = Y[valid_idx]
        ID = ID[valid_idx]
        G = G[valid_idx]
        period = period[valid_idx]
        if (n_controls > 0) {
            X = X[valid_idx, .]
        }
    }
    
    // ========================================================================
    // Filter to pretreatment periods
    // Data is filtered prior to estimation to ensure only pretreatment
    // observations are used when computing placebo coefficients
    // ========================================================================
    
    // Check if pretreat_periods is a sentinel value (1x1 matrix with -999999)
    has_pretreat = !(rows(pretreat_periods) == 1 && pretreat_periods[1] == -999999)
    
    if (has_pretreat) {
        // Filter data to user-specified pretreatment periods
        pretreat_idx = J(rows(Y), 1, 0)
        for (i = 1; i <= rows(pretreat_periods); i++) {
            pretreat_idx = pretreat_idx :| (period :== pretreat_periods[i])
        }
        pretreat_idx = selectindex(pretreat_idx)
        
        Y = Y[pretreat_idx]
        ID = ID[pretreat_idx]
        G = G[pretreat_idx]
        period = period[pretreat_idx]
        if (n_controls > 0) {
            X = X[pretreat_idx, .]
        }
    }
    
    // ========================================================================
    // Check for empty data after pretreatment filtering
    // ========================================================================
    if (rows(Y) == 0) {
        errprintf("Error: No observations remain after filtering to pretreatment periods.\n")
        errprintf("Possible causes:\n")
        errprintf("  - The specified pretreatment periods do not exist in the data\n")
        errprintf("  - All matching observations were removed due to missing values\n")
        _error(3498)
    }
    
    // Determine placebo periods (pretreatment periods excluding base period)
    if (has_pretreat && baseperiod != -999999) {
        placebo_periods = select(pretreat_periods, pretreat_periods :!= baseperiod)
        base_period_val = baseperiod
    }
    else if (has_pretreat) {
        base_period_val = max(pretreat_periods)
        placebo_periods = select(pretreat_periods, pretreat_periods :!= base_period_val)
    }
    else {
        // No pretreatment specified - use all periods
        real colvector unique_periods
        unique_periods = uniqrows(period)
        // Check if user specified baseperiod; if not, use maximum period
        if (baseperiod != -999999) {
            base_period_val = baseperiod
        }
        else {
            base_period_val = max(unique_periods)
        }
        placebo_periods = select(unique_periods, unique_periods :!= base_period_val)
    }
    
    // Create placebo indicators: placebo_t = G * I(period == t)
    // This is the interaction between treatment group and time period
    placebo_vars = J(rows(Y), rows(placebo_periods), 0)
    for (i = 1; i <= rows(placebo_periods); i++) {
        placebo_vars[., i] = G :* (period :== placebo_periods[i])
    }
    
    // Build full X matrix
    if (n_controls > 0) {
        X_full = placebo_vars, X
    }
    else {
        X_full = placebo_vars
    }
    
    // Double demean
    WD = _eqt_construct_WD(period, ID)
    Y_dm = _eqt_double_demean(Y, ID, period, WD)
    X_full_dm = _eqt_double_demean(X_full, ID, period, WD)
    
    // Check and remove multicollinearity before computing placebo_coefs
    real rowvector problematic_wd, problematic_x
    real scalar effective_no_placebos, idx
    real rowvector removed_placebo_idx, removed_control_idx
    string scalar period_list
    
    // Variables for preserving real period labels
    real colvector effective_placebo_periods
    string colvector placebo_rownames
    
    // Check for multicollinearity in WD matrix
    problematic_wd = _eqt_detect_multicol_qr(WD)
    if (cols(problematic_wd) > 0) {
        WD = _eqt_remove_multicol(WD, problematic_wd)
        // Re-demean with updated WD
        Y_dm = _eqt_double_demean(Y, ID, period, WD)
        X_full_dm = _eqt_double_demean(X_full, ID, period, WD)
    }
    
    // Check for multicollinearity in X_full_dm
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
        
        // Also remove from raw X_full for passing to bootstrap test
        X_full = _eqt_remove_multicol(X_full, problematic_x)
        
        // Update effective number of placebos
        effective_no_placebos = no_placebos - cols(removed_placebo_idx)
        
        // Update effective_placebo_periods by removing multicollinear periods
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
    
    // Compute unconstrained OLS (using reduced X if multicollinearity was found)
    unconstrained_coefs = _eqt_ols_cholesky(cross(X_full_dm, X_full_dm), cross(X_full_dm, Y_dm))
    
    // Extract placebo coefficients (using effective_no_placebos)
    placebo_coefs = unconstrained_coefs[1::effective_no_placebos]
    
    // Compute max absolute coefficient
    max_abs_coef = max(abs(placebo_coefs))
    
    // Compute sigma_hathat_c
    sigma_hathat_c = _eqt_sigma_hathat_c(unconstrained_coefs, X_full_dm, Y_dm, ID, period)
    
    // Initialize return values
    critical_value = .
    reject = .
    min_delta = .
    
    if (delta_specified == 1) {
        // Run bootstrap test with specified delta
        // Pass raw data (Y, X_full), not demeaned data
        // _eqt_bootstrap_test will handle double-demeaning internally
        // Note: X_full may have been reduced by multicollinearity removal above
        boot_result = _eqt_bootstrap_test(Y, X_full, ID, period, effective_no_placebos,
                                           delta, alpha, B, type, seed, x_varnames)
        
        critical_value = boot_result.critical_value
        reject = boot_result.reject_H0
        // Use sigma_hathat_c from bootstrap result (computed with constrained coefficients)
        sigma_hathat_c = boot_result.sigma_hathat_c
    }
    else {
        // Compute minimum threshold using bootstrap method
        // Pass raw data - _eqt_min_delta_bootstrap handles demeaning internally
        // Note: X_full may have been reduced by multicollinearity removal above
        min_delta_result = _eqt_min_delta_bootstrap(Y, X_full, ID, period, effective_no_placebos,
                                                     alpha, B, type, seed)
        
        min_delta = min_delta_result.min_delta
    }
    
    // Build row names using real period labels (e.g., "placebo_1", "placebo_2", "placebo_5")
    placebo_rownames = J(effective_no_placebos, 1, "")
    for (i = 1; i <= effective_no_placebos; i++) {
        placebo_rownames[i] = "placebo_" + strofreal(effective_placebo_periods[i])
    }
    
    // Return results to Stata
    st_matrix("r(placebo_coefs)", placebo_coefs)
    st_matrixrowstripe("r(placebo_coefs)", (J(effective_no_placebos, 1, ""), placebo_rownames))
    st_matrixcolstripe("r(placebo_coefs)", ("", "coef"))
    
    st_numscalar("r(max_abs_coef)", max_abs_coef)
    st_numscalar("r(critical_value)", critical_value)
    st_numscalar("r(reject)", reject)
    st_numscalar("r(min_delta)", min_delta)
    st_numscalar("r(sigma_hathat_c)", sigma_hathat_c)
    
    // Return effective_no_placebos (actual number of placebos after multicollinearity removal)
    st_numscalar("r(no_placebos)", effective_no_placebos)
    
    // Return placebo_periods for display use
    st_matrix("r(placebo_periods)", effective_placebo_periods)
}

end
