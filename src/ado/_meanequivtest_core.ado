*! _meanequivtest_core.ado - Core computation for the mean equivalence test
*!
*! This module implements the folded normal test for hypothesis (3.2) from
*! Dette & Schumann (2024):
*!   H0: |beta_bar| >= tau  vs  H1: |beta_bar| < tau
*! where beta_bar = (1/T) * sum_{t=1}^{T} beta_t is the mean placebo effect.
*!
*! Test statistic: |beta_bar_hat| = |d'beta_hat| with d = (1/T, ..., 1/T)'
*! Variance estimator: sigma_hat^2 = d' * V_hat * d
*! Rejection rule: Reject H0 when |beta_bar_hat| < Q_{FN(tau, sigma_hat^2)}(alpha)
*! where Q_{FN} denotes the quantile function of the folded normal distribution.
*!
*! The minimum equivalence threshold tau* is computed when no threshold is
*! specified, representing the smallest tau for which H0 can be rejected.

program define _meanequivtest_core, rclass
    version 16.0
    
    syntax, Y(varname numeric) ID(varname numeric) G(varname numeric) Period(varname numeric) ///
        NOplacebos(integer) ///
        THRESHOLDspecified(integer) THRESHold(real) ///
        Alpha(real) VCE(string) ///
        [X(varlist numeric)] [Cluster(varname numeric)] ///
        [PREtreatment(numlist)] [BASEperiod(real -999999)]
    
    local no_placebos = `noplacebos'
    local threshold_specified = `thresholdspecified'
    
    // -------------------------------------------------------------------------
    // Pre-treatment Period Processing
    // -------------------------------------------------------------------------
    
    tempname pretreat_mat
    if "`pretreatment'" != "" {
        local n_pretreat : word count `pretreatment'
        matrix `pretreat_mat' = J(`n_pretreat', 1, .)
        local i = 1
        foreach val of local pretreatment {
            matrix `pretreat_mat'[`i', 1] = `val'
            local i = `i' + 1
        }
    }
    else {
        // Sentinel value indicating no explicit pretreatment specification
        matrix `pretreat_mat' = J(1, 1, -999999)
    }
    
    // -------------------------------------------------------------------------
    // Mata Core Computation
    // -------------------------------------------------------------------------
    
    local x_vars "`x'"
    
    capture noisily mata: _meanequivtest_core_mata("`y'", "`x_vars'", "`id'", "`g'", "`period'", ///
        `no_placebos', `threshold_specified', `threshold', `alpha', ///
        "`vce'", "`cluster'", st_matrix("`pretreat_mat'"), `baseperiod')
    
    if _rc {
        display as error "meanequivtest: Mata computation failed"
        display as error "  Error code: " _rc
        exit _rc
    }
    
    // -------------------------------------------------------------------------
    // Validate Return Values
    // -------------------------------------------------------------------------
    
    capture confirm matrix r(placebo_coefs)
    if _rc {
        display as error "meanequivtest: Mata function did not return placebo_coefs matrix"
        exit 498
    }
    
    capture confirm matrix r(vcov_placebo)
    if _rc {
        display as error "meanequivtest: Mata function did not return vcov_placebo matrix"
        exit 498
    }
    
    if missing(r(abs_mean_placebo)) {
        display as error "meanequivtest: Mata function returned invalid abs_mean_placebo"
        display as error "  Possible causes:"
        display as error "    - Data contains missing values in key variables"
        display as error "    - VCE matrix computation failed"
        display as error "    - Numerical computation error (e.g., singular matrix)"
        exit 498
    }
    
    // -------------------------------------------------------------------------
    // Process and Return Results
    // -------------------------------------------------------------------------
    
    tempname placebo_coefs_mat vcov_placebo_mat placebo_periods_mat
    matrix `placebo_coefs_mat' = r(placebo_coefs)
    matrix `vcov_placebo_mat' = r(vcov_placebo)
    matrix `placebo_periods_mat' = r(placebo_periods)
    
    local abs_mean_placebo = r(abs_mean_placebo)
    local var_mean_placebo = r(var_mean_placebo)
    local se_mean_placebo = r(se_mean_placebo)
    
    // Effective dimension after multicollinearity adjustment
    local effective_no_placebos = r(no_placebos)
    
    // Construct placebo coefficient labels
    local effective_placebo_names ""
    forvalues i = 1/`effective_no_placebos' {
        local period_val = `placebo_periods_mat'[`i', 1]
        local effective_placebo_names "`effective_placebo_names' placebo_`period_val'"
    }
    local effective_placebo_names = strtrim("`effective_placebo_names'")
    
    if `threshold_specified' == 1 {
        local critical_value = r(critical_value)
        local p_value = r(p_value)
        local reject = r(reject)
    }
    else {
        local min_threshold = r(min_threshold)
    }
    
    // Return matrices
    return matrix placebo_coefs = `placebo_coefs_mat'
    return matrix vcov_placebo = `vcov_placebo_mat'
    
    // Return scalars
    return scalar abs_mean_placebo = `abs_mean_placebo'
    return scalar var_mean_placebo = `var_mean_placebo'
    return scalar se_mean_placebo = `se_mean_placebo'
    return scalar no_placebos = `effective_no_placebos'
    return local placebo_names = "`effective_placebo_names'"
    
    if `threshold_specified' == 1 {
        return scalar critical_value = `critical_value'
        return scalar p_value = `p_value'
        return scalar reject = `reject'
    }
    else {
        return scalar min_threshold = `min_threshold'
    }
    
end
