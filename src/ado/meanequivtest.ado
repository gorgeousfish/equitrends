*! meanequivtest.ado - Mean equivalence test for pre-trends in DiD estimation
*!
*! Implements the mean equivalence test for assessing the parallel trends
*! assumption in Difference-in-Differences estimation. This test evaluates
*! whether the average placebo coefficient lies within an equivalence bound.
*!
*! Hypothesis:
*!   H0: |beta_bar| >= tau  vs  H1: |beta_bar| < tau
*!   where beta_bar = (1/T) * sum_{l=1}^{T} beta_l
*!
*! The test statistic |beta_bar| is compared against the alpha-quantile of
*! a folded normal distribution with mean tau and variance sigma^2,
*! where sigma^2 = 1' * Sigma * 1 / n.

program define meanequivtest, eclass
    version 16.0
    
    equitrends_init
    
    local cmdline "meanequivtest `0'"
    
    // -------------------------------------------------------------------------
    // Syntax Parsing
    // -------------------------------------------------------------------------
    
    // Parse command options
    // Note: Group() accepts g/gr/gro/grou/group; Time() accepts t/ti/tim/time
    // Period() is an alias for Time() for backward compatibility
    syntax varname(numeric), ID(varname numeric) ///
        [Group(varname numeric) G(varname numeric)] ///
        [Time(varname numeric) Period(varname numeric)] ///
        [X(varlist numeric)] [Cluster(varname numeric)] ///
        [PREtreatment(numlist)] [BASEperiod(real -999999)] ///
        [THRESHold(real -999999)] [Alpha(real 0.05)] ///
        [VCE(string)] [NODISplay]
    
    local y `varlist'
    
    // Handle Group/G alias
    if "`group'" != "" & "`g'" != "" {
        display as error "may not specify both group() and g()"
        exit 198
    }
    if "`group'" == "" & "`g'" == "" {
        display as error "option group() or g() required"
        exit 198
    }
    if "`group'" != "" {
        local g `group'
    }
    
    // Handle Time/Period alias
    if "`time'" != "" & "`period'" != "" {
        display as error "may not specify both time() and period()"
        exit 198
    }
    if "`time'" == "" & "`period'" == "" {
        display as error "option time() or period() required"
        exit 198
    }
    if "`time'" != "" {
        local period `time'
    }
    
    // -------------------------------------------------------------------------
    // Parameter Validation
    // -------------------------------------------------------------------------
    
    local parse_cmd "_meanequivtest_parse, y(`y') id(`id') g(`g') period(`period')"
    
    if "`x'" != "" {
        local parse_cmd "`parse_cmd' x(`x')"
    }
    if "`cluster'" != "" {
        local parse_cmd "`parse_cmd' cluster(`cluster')"
    }
    if "`pretreatment'" != "" {
        local parse_cmd "`parse_cmd' pretreatment(`pretreatment')"
    }
    if `baseperiod' != -999999 {
        local parse_cmd "`parse_cmd' baseperiod(`baseperiod')"
    }
    if `threshold' != -999999 {
        local parse_cmd "`parse_cmd' threshold(`threshold')"
    }
    local parse_cmd "`parse_cmd' alpha(`alpha')"
    if "`vce'" != "" {
        local parse_cmd "`parse_cmd' vce(`vce')"
    }
    
    `parse_cmd'
    
    local threshold_specified = r(threshold_specified)
    local threshold_val = r(threshold)
    local alpha_val = r(alpha)
    local vce_type = r(vce_type)
    
    // Extract cluster variable from vce() option if not specified separately
    local cluster_parsed "`r(cluster)'"
    if "`cluster_parsed'" != "" & "`cluster_parsed'" != "." & "`cluster'" == "" {
        local cluster "`cluster_parsed'"
    }
    
    // -------------------------------------------------------------------------
    // Data Processing
    // -------------------------------------------------------------------------
    
    local dataproc_cmd "_equitrends_dataproc, y(`y') id(`id') g(`g') period(`period')"
    
    if "`x'" != "" {
        local dataproc_cmd "`dataproc_cmd' x(`x')"
    }
    if "`cluster'" != "" {
        local dataproc_cmd "`dataproc_cmd' cluster(`cluster')"
    }
    if "`pretreatment'" != "" {
        local dataproc_cmd "`dataproc_cmd' pretreatment(`pretreatment')"
    }
    if `baseperiod' != -999999 {
        local dataproc_cmd "`dataproc_cmd' baseperiod(`baseperiod')"
    }
    
    quietly `dataproc_cmd'
    
    local baseperiod_val = r(baseperiod)
    local is_balanced = r(balanced_panel)
    local no_placebos = r(no_placebos)
    local n_individuals = r(n)
    local n_periods = r(n_t)
    local N_obs = r(N_obs)
    local dp_placebo_names "`r(placebo_names)'"
    local dp_no_periods "`r(no_periods)'"
    
    // Extract T_min and T_max for unbalanced panels
    local dp_t_min = -999999
    local dp_t_max = -999999
    if `is_balanced' == 0 {
        local dp_t_min : word 1 of `dp_no_periods'
        local dp_t_max : word 2 of `dp_no_periods'
    }
    
    // -------------------------------------------------------------------------
    // Core Computation
    // Mean placebo coefficient is computed and tested against the folded
    // normal distribution quantile.
    // -------------------------------------------------------------------------
    
    local core_cmd "_meanequivtest_core, y(`y') id(`id') g(`g') period(`period')"
    local core_cmd "`core_cmd' noplacebos(`no_placebos')"
    local core_cmd "`core_cmd' thresholdspecified(`threshold_specified') threshold(`threshold_val')"
    local core_cmd "`core_cmd' alpha(`alpha_val') vce(`vce_type')"
    
    if "`x'" != "" {
        local core_cmd "`core_cmd' x(`x')"
    }
    if "`cluster'" != "" {
        local core_cmd "`core_cmd' cluster(`cluster')"
    }
    if "`pretreatment'" != "" {
        local core_cmd "`core_cmd' pretreatment(`pretreatment')"
    }
    if `baseperiod_val' != -999999 {
        local core_cmd "`core_cmd' baseperiod(`baseperiod_val')"
    }
    
    `core_cmd'
    
    // Extract test statistics
    local abs_mean_placebo = r(abs_mean_placebo)
    local var_mean_placebo = r(var_mean_placebo)
    local se_mean_placebo = r(se_mean_placebo)
    
    // Adjust for collinear placebo variables removed during estimation
    local effective_no_placebos = r(no_placebos)
    local effective_placebo_names "`r(placebo_names)'"
    
    tempname b_placebo V_placebo
    matrix `b_placebo' = r(placebo_coefs)
    matrix `V_placebo' = r(vcov_placebo)
    
    if `threshold_specified' == 1 {
        local critical_value = r(critical_value)
        local p_value = r(p_value)
        local reject = r(reject)
    }
    else {
        local min_threshold = r(min_threshold)
    }
    
    // -------------------------------------------------------------------------
    // Results Display
    // -------------------------------------------------------------------------
    
    if "`nodisplay'" == "" {
        if `threshold_specified' == 1 {
            _meanequivtest_display, abs_mean_placebo(`abs_mean_placebo') ///
                se_mean_placebo(`se_mean_placebo') ///
                threshold_specified(`threshold_specified') ///
                threshold(`threshold_val') critical_value(`critical_value') ///
                p_value(`p_value') reject(`reject') ///
                alpha(`alpha_val') n(`N_obs') n_g(`n_individuals') n_t(`n_periods') ///
                no_placebos(`effective_no_placebos') baseperiod(`baseperiod_val') ///
                balanced(`is_balanced') vce(`vce_type') ///
                t_min(`dp_t_min') t_max(`dp_t_max')
        }
        else {
            _meanequivtest_display, abs_mean_placebo(`abs_mean_placebo') ///
                se_mean_placebo(`se_mean_placebo') ///
                threshold_specified(`threshold_specified') ///
                min_threshold(`min_threshold') ///
                alpha(`alpha_val') n(`N_obs') n_g(`n_individuals') n_t(`n_periods') ///
                no_placebos(`effective_no_placebos') baseperiod(`baseperiod_val') ///
                balanced(`is_balanced') vce(`vce_type') ///
                t_min(`dp_t_min') t_max(`dp_t_max')
        }
    }
    
    // -------------------------------------------------------------------------
    // Store Estimation Results
    // -------------------------------------------------------------------------
    
    ereturn clear
    
    // Coefficient vector and variance-covariance matrix
    ereturn matrix b_placebo = `b_placebo'
    ereturn matrix V_placebo = `V_placebo'
    
    // Sample and test information
    ereturn scalar N = `N_obs'
    ereturn scalar N_g = `n_individuals'
    ereturn scalar T = `n_periods'
    ereturn scalar no_placebos = `effective_no_placebos'
    ereturn scalar abs_mean_placebo = `abs_mean_placebo'
    ereturn scalar var_mean_placebo = `var_mean_placebo'
    ereturn scalar se_mean_placebo = `se_mean_placebo'
    ereturn scalar alpha = `alpha_val'
    ereturn scalar base_period = `baseperiod_val'
    ereturn scalar is_balanced = `is_balanced'
    
    // Store T_min and T_max for unbalanced panels
    if `is_balanced' == 0 {
        ereturn scalar T_min = `dp_t_min'
        ereturn scalar T_max = `dp_t_max'
    }
    
    if `threshold_specified' == 1 {
        ereturn scalar threshold = `threshold_val'
        ereturn scalar critical_value = `critical_value'
        ereturn scalar p_value = `p_value'
        ereturn scalar reject = `reject'
    }
    else {
        ereturn scalar min_threshold = `min_threshold'
    }
    
    // Command and variable information
    ereturn local cmd "meanequivtest"
    ereturn local cmdline "`cmdline'"
    ereturn local depvar "`y'"
    ereturn local idvar "`id'"
    ereturn local groupvar "`g'"
    ereturn local timevar "`period'"
    ereturn local vce "`vce_type'"
    ereturn scalar threshold_specified = `threshold_specified'
    
    if "`cluster'" != "" {
        ereturn local clustvar "`cluster'"
    }
    
    ereturn local placebo_names "`effective_placebo_names'"
    
end
