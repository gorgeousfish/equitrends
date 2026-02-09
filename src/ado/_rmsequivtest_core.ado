*! _rmsequivtest_core.ado - Core computation for RMS equivalence test
*!
*! Implements the pivotal self-normalized test for the RMS hypothesis (3.3):
*!   H0: β_RMS ≥ ζ  vs  H1: β_RMS < ζ
*! where β_RMS = ||β||/√T is the root mean squared placebo coefficient.
*!
*! The test statistic is constructed using sequential subsampling:
*!
*!   M̂_n = (β̂²_RMS(1) - β²_RMS) / V̂_n
*!
*! where the self-normalized variance V̂_n is computed from subsample estimates:
*!
*!   V̂_n = [∫(β̂²_RMS(λ) - β̂²_RMS(1))² ν(dλ)]^{1/2}
*!
*! The null hypothesis is rejected when:
*!
*!   β̂²_RMS < ζ² + Q_W(α) · V̂_n
*!
*! where Q_W(α) is the α-quantile of the limiting W distribution.
*!
*! Arguments:
*!   y, id, g, period : varnames - Panel data identifiers
*!   no_placebos      : integer  - Number of placebo coefficients (T)
*!   threshold        : real     - Equivalence threshold ζ > 0
*!   alpha            : real     - Significance level
*!   nolambda         : integer  - Number of subsample fractions
*!   level            : real     - Confidence level for CI (default 95)
*!
*! Returns:
*!   r(MS_placebo)    : scalar - Mean squared placebo coefficient β̂²_RMS
*!   r(RMS_placebo)   : scalar - Root mean squared coefficient √(β̂²_RMS)
*!   r(V_n)           : scalar - Self-normalized variance estimate
*!   r(Q_W)           : scalar - Critical value from W distribution
*!   r(placebo_coefs) : matrix - Estimated placebo coefficients β̂
*!   r(MS_lambda)     : matrix - Subsample MS estimates
*!   r(reject)        : scalar - Test decision (1 = reject H0)
*!   r(min_threshold) : scalar - Minimum equivalence threshold ζ*
*!   r(RMS_ci_lower)  : scalar - Lower bound of RMS CI
*!   r(RMS_ci_upper)  : scalar - Upper bound of RMS CI
*!   r(MS_ci_lower)   : scalar - Lower bound of MS CI
*!   r(MS_ci_upper)   : scalar - Upper bound of MS CI
*!   r(level)         : scalar - Confidence level

version 16.0

program define _rmsequivtest_core, rclass
    
    syntax, Y(varname numeric) ID(varname numeric) G(varname numeric) Period(varname numeric) ///
        NOPlacebos(integer) ///
        THRESHOLDSpecified(integer) THRESHold(real) ///
        Alpha(real) NOLambda(integer) LEVel(real) ///
        [X(varlist numeric)] [PREperiod(numlist)] [BASEperiod(real -999999)] ///
        [SEED(integer -999999)]
    
    // =========================================================================
    // Random Seed Initialization
    // =========================================================================
    
    if `seed' != -999999 {
        set seed `seed'
    }
    
    // =========================================================================
    // Pre-treatment Period Matrix Construction
    // =========================================================================
    
    tempname pretreat_mat
    if "`preperiod'" != "" {
        local n_pretreat : word count `preperiod'
        matrix `pretreat_mat' = J(`n_pretreat', 1, .)
        local i = 1
        foreach val of local preperiod {
            matrix `pretreat_mat'[`i', 1] = `val'
            local i = `i' + 1
        }
    }
    else {
        // Sentinel value indicates no user-specified pre-treatment periods
        matrix `pretreat_mat' = J(1, 1, -999999)
    }
    
    // =========================================================================
    // Mata Core Computation
    // =========================================================================
    
    local x_vars "`x'"
    
    capture noisily mata: _rmsequivtest_core_mata("`y'", "`x_vars'", "`id'", "`g'", "`period'", ///
        `noplacebos', `thresholdspecified', `threshold', `alpha', ///
        `nolambda', `level', st_matrix("`pretreat_mat'"), `baseperiod')
    
    if _rc {
        display as error "rmsequivtest: Mata computation failed"
        display as error "  Error code: " _rc
        exit _rc
    }
    
    // =========================================================================
    // Return Value Validation
    // =========================================================================
    
    capture confirm matrix r(placebo_coefs)
    if _rc {
        display as error "rmsequivtest: Mata function did not return placebo_coefs matrix"
        exit 498
    }
    
    capture confirm matrix r(MS_lambda)
    if _rc {
        display as error "rmsequivtest: Mata function did not return MS_lambda matrix"
        exit 498
    }
    
    if missing(r(MS_placebo)) {
        display as error "rmsequivtest: Mata function returned invalid MS_placebo"
        display as error "  Possible causes:"
        display as error "    - Data contains missing values in key variables"
        display as error "    - Subsample estimation failed"
        display as error "    - Numerical computation error"
        exit 498
    }
    
    // =========================================================================
    // Return Results
    // =========================================================================
    
    tempname placebo_coefs_mat MS_lambda_mat
    matrix `placebo_coefs_mat' = r(placebo_coefs)
    matrix `MS_lambda_mat' = r(MS_lambda)
    
    local MS_placebo = r(MS_placebo)
    local RMS_placebo = r(RMS_placebo)
    local V_n = r(V_n)
    local Q_W = r(Q_W)
    
    // CI results
    local RMS_ci_lower = r(RMS_ci_lower)
    local RMS_ci_upper = r(RMS_ci_upper)
    local MS_ci_lower = r(MS_ci_lower)
    local MS_ci_upper = r(MS_ci_upper)
    local level_ret = r(level)
    
    if `thresholdspecified' == 1 {
        local MS_critical = r(MS_critical)
        local RMS_critical = r(RMS_critical)
        local reject = r(reject)
    }
    else {
        local min_threshold = r(min_threshold)
    }
    
    return matrix placebo_coefs = `placebo_coefs_mat'
    return matrix MS_lambda = `MS_lambda_mat'
    
    return scalar MS_placebo = `MS_placebo'
    return scalar RMS_placebo = `RMS_placebo'
    return scalar V_n = `V_n'
    return scalar Q_W = `Q_W'
    
    // Return CI results
    return scalar RMS_ci_lower = `RMS_ci_lower'
    return scalar RMS_ci_upper = `RMS_ci_upper'
    return scalar MS_ci_lower = `MS_ci_lower'
    return scalar MS_ci_upper = `MS_ci_upper'
    return scalar level = `level_ret'
    
    if `thresholdspecified' == 1 {
        return scalar MS_critical = `MS_critical'
        return scalar RMS_critical = `RMS_critical'
        return scalar reject = `reject'
    }
    else {
        return scalar min_threshold = `min_threshold'
    }
    
end


// ============================================================================
// Mata Implementation
// ============================================================================
mata:

// ----------------------------------------------------------------------------
// _rmsequivtest_core_mata()
// Core computation for RMS equivalence test using self-normalized subsampling.
//
// The function estimates placebo coefficients from the TWFE model and computes
// the self-normalized test statistic for the RMS hypothesis.
// ----------------------------------------------------------------------------
void _rmsequivtest_core_mata(string scalar y_name, string scalar x_names,
                              string scalar id_name, string scalar g_name,
                              string scalar period_name,
                              real scalar no_placebos, real scalar threshold_specified,
                              real scalar threshold, real scalar alpha,
                              real scalar no_lambda, real scalar level,
                              real colvector pretreat_periods, real scalar baseperiod)
{
    real colvector Y, ID, G, period
    real matrix X
    real colvector placebo_periods, unique_periods
    real scalar base_period_val, i
    string rowvector x_varnames
    real colvector pretreat_idx
    real colvector valid_idx
    real scalar na_omitted
    
    // =========================================================================
    // Data Loading
    // =========================================================================
    
    Y = st_data(., y_name)
    ID = st_data(., id_name)
    G = st_data(., g_name)
    period = st_data(., period_name)
    
    if (x_names != "") {
        x_varnames = tokens(x_names)
        X = st_data(., x_varnames)
    }
    else {
        X = J(rows(Y), 0, .)
    }
    
    // =========================================================================
    // Missing Value Handling
    // =========================================================================
    
    valid_idx = _eqt_omit_na(Y, ID, G, period, X)
    na_omitted = rows(Y) - rows(valid_idx)
    if (na_omitted > 0) {
        printf("{txt}Warning: %g rows of pre-treatment data omitted due to NAs.\n", na_omitted)
        Y = Y[valid_idx]
        ID = ID[valid_idx]
        G = G[valid_idx]
        period = period[valid_idx]
        if (cols(X) > 0) {
            X = X[valid_idx, .]
        }
    }
    
    // =========================================================================
    // Base Period and Placebo Period Determination
    // =========================================================================
    
    real scalar has_pretreat
    has_pretreat = !(rows(pretreat_periods) == 1 && pretreat_periods[1] == -999999)
    
    if (has_pretreat && baseperiod != -999999) {
        // Both pre-treatment periods and base period specified by user
        placebo_periods = select(pretreat_periods, pretreat_periods :!= baseperiod)
        base_period_val = baseperiod
    }
    else if (has_pretreat) {
        // Pre-treatment periods specified; base period defaults to maximum
        base_period_val = max(pretreat_periods)
        placebo_periods = select(pretreat_periods, pretreat_periods :!= base_period_val)
    }
    else {
        // No pre-treatment specified; use all available periods
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
    // Pre-treatment Period Filtering
    // Base period is included for TWFE double-demeaning transformation
    // =========================================================================
    
    if (has_pretreat) {
        pretreat_idx = J(rows(Y), 1, 0)
        for (i = 1; i <= rows(pretreat_periods); i++) {
            pretreat_idx = pretreat_idx :| (period :== pretreat_periods[i])
        }
        pretreat_idx = pretreat_idx :| (period :== base_period_val)
        pretreat_idx = selectindex(pretreat_idx)
        
        Y = Y[pretreat_idx]
        ID = ID[pretreat_idx]
        G = G[pretreat_idx]
        period = period[pretreat_idx]
        if (cols(X) > 0) {
            X = X[pretreat_idx, .]
        }
    }
    
    // =========================================================================
    // Data Availability Validation
    // =========================================================================
    
    if (rows(Y) == 0) {
        errprintf("Error: No observations remain after filtering to pretreatment periods.\n")
        errprintf("Possible causes:\n")
        errprintf("  - The specified pretreatment periods do not exist in the data\n")
        errprintf("  - All matching observations were removed due to missing values\n")
        _error(3498)
    }
    
    // Sort placebo periods for consistent output ordering
    // Test statistics (MS, RMS, critical values) are order-invariant
    placebo_periods = sort(placebo_periods, 1)
    
    // =========================================================================
    // RMS Test Execution
    // =========================================================================
    
    _eqt_rms_test(Y, X, ID, G, period, placebo_periods, base_period_val,
                  threshold_specified, threshold, alpha, no_lambda, level)
}

end
