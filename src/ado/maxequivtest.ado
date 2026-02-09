*! maxequivtest.ado - Maximum equivalence test for pre-trends in DiD estimation
*!
*! Implements the equivalence hypothesis (3.1) from Dette and Schumann (2024):
*!   H0: ||beta||_inf >= delta  vs.  H1: ||beta||_inf < delta
*! where ||beta||_inf := max_{l=1,...,T} |beta_l| is the supremum norm of placebo
*! coefficients. Rejection provides statistical evidence that the maximum deviation
*! from parallel trends is bounded by the equivalence threshold delta.
*!
*! Three testing procedures are available:
*!   - IU:   Intersection-union test using folded normal quantiles (Eq. 4.4)
*!   - Boot: Parametric bootstrap under constrained estimation (Theorem 1)
*!   - Wild: Wild cluster bootstrap for clustered/heteroskedastic errors (Remark 1c)

program define maxequivtest, eclass
    version 16.0
    
    // Initialize Mata functions (compile on first use if needed)
    equitrends_init
    
    local cmdline "maxequivtest `0'"
    
    // -------------------------------------------------------------------------
    // Syntax Parsing and Input Validation
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
        [Method(string)] [Reps(integer 1000)] [Seed(integer 0)] ///
        [VCE(string)] [NODISplay] [NODOTS]
    
    local Y `varlist'
    
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
    
    _maxequivtest_parse, y(`Y') id(`id') g(`g') period(`period') ///
        x(`x') cluster(`cluster') pretreatment(`pretreatment') ///
        baseperiod(`baseperiod') threshold(`threshold') alpha(`alpha') ///
        method(`method') reps(`reps') seed(`seed') vce(`vce')
    
    local method = r(method)
    local alpha = r(alpha)
    local reps = r(reps)
    local seed = r(seed)
    local threshold_specified = r(threshold_specified)
    local threshold_val = r(threshold)
    local vce_type = r(vce_type)
    local cluster_var "`r(cluster)'"
    
    if "`cluster_var'" != "" & "`cluster_var'" != "." {
        local cluster "`cluster_var'"
    }
    
    // -------------------------------------------------------------------------
    // Data Preparation
    // -------------------------------------------------------------------------
    // Placebo indicators W_{i,t} = G_i * D_l(t) are constructed and double-demeaned
    // for two-way fixed effects estimation following the TWFE model (Eq. 2.5-2.6).
    
    local dataproc_cmd "_equitrends_dataproc, y(`Y') id(`id') g(`g') period(`period')"
    
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
    
    local dp_n = r(n)
    local dp_n_t = r(n_t)
    local dp_N_obs = r(N_obs)
    local dp_baseperiod = r(baseperiod)
    local dp_no_placebos = r(no_placebos)
    local dp_balanced = r(balanced_panel)
    local dp_placebo_names = "`r(placebo_names)'"
    local dp_no_periods = "`r(no_periods)'"
    
    // Extract period range for unbalanced panel display
    local dp_t_min = -999999
    local dp_t_max = -999999
    if `dp_balanced' == 0 {
        local dp_t_min : word 1 of `dp_no_periods'
        local dp_t_max : word 2 of `dp_no_periods'
    }
    
    // -------------------------------------------------------------------------
    // Test Computation
    // -------------------------------------------------------------------------
    // The test statistic ||beta_hat||_inf is computed and compared against the
    // appropriate critical value based on the selected method:
    //   IU:        Folded normal alpha-quantile Q_{N_F(delta, Sigma_tt/n)}(alpha)
    //   Bootstrap: Empirical alpha-quantile from constrained bootstrap distribution
    
    local core_cmd "_maxequivtest_core, y(`Y') id(`id') g(`g') period(`period')"
    local core_cmd "`core_cmd' noplacebos(`dp_no_placebos')"
    local core_cmd "`core_cmd' thresholdspecified(`threshold_specified') threshold(`threshold_val')"
    local core_cmd "`core_cmd' alpha(`alpha') method(`method') reps(`reps') seed(`seed')"
    local core_cmd "`core_cmd' vce(`vce_type')"
    
    if "`x'" != "" {
        local core_cmd "`core_cmd' x(`x')"
    }
    if "`cluster'" != "" {
        local core_cmd "`core_cmd' cluster(`cluster')"
    }
    if "`pretreatment'" != "" {
        local core_cmd "`core_cmd' pretreatment(`pretreatment')"
    }
    if `baseperiod' != -999999 {
        local core_cmd "`core_cmd' baseperiod(`baseperiod')"
    }
    if "`nodots'" != "" {
        local core_cmd "`core_cmd' nodots"
    }
    
    `core_cmd'
    
    local max_abs_coef = r(max_abs_coef)
    local critical_value = r(critical_value)
    local reject = r(reject)
    local min_threshold = r(min_threshold)
    local sigma_hathat_c = r(sigma_hathat_c)
    
    // Placebo count may be reduced if multicollinearity is detected
    local effective_no_placebos = r(no_placebos)
    local effective_placebo_names "`r(placebo_names)'"
    
    tempname placebo_coefs_mat placebo_se_mat vcov_placebo_mat crit_values_mat min_thresholds_mat
    matrix `placebo_coefs_mat' = r(placebo_coefs)
    
    if "`method'" == "IU" {
        matrix `placebo_se_mat' = r(placebo_se)
        matrix `vcov_placebo_mat' = r(vcov_placebo)
        if `threshold_specified' == 1 {
            matrix `crit_values_mat' = r(crit_values)
        }
        else {
            matrix `min_thresholds_mat' = r(min_thresholds)
        }
    }
    
    // -------------------------------------------------------------------------
    // Output Display
    // -------------------------------------------------------------------------
    
    if "`nodisplay'" == "" {
        local placebo_coefs_str ""
        forvalues i = 1/`effective_no_placebos' {
            local coef_val = `placebo_coefs_mat'[`i', 1]
            local placebo_coefs_str "`placebo_coefs_str' `coef_val'"
        }
        local placebo_coefs_str = strtrim("`placebo_coefs_str'")
        
        _maxequivtest_display, method(`method') ///
            threshold_specified(`threshold_specified') threshold(`threshold_val') ///
            alpha(`alpha') reps(`reps') ///
            max_abs_coef(`max_abs_coef') critical_value(`critical_value') ///
            reject(`reject') min_threshold(`min_threshold') ///
            n(`dp_n') n_t(`dp_n_t') n_total(`dp_N_obs') ///
            no_placebos(`effective_no_placebos') placebo_names(`effective_placebo_names') ///
            placebo_coefs(`placebo_coefs_str') ///
            balanced(`dp_balanced') baseperiod(`dp_baseperiod') ///
            t_min(`dp_t_min') t_max(`dp_t_max')
    }
    
    // -------------------------------------------------------------------------
    // Return Values
    // -------------------------------------------------------------------------
    
    ereturn clear
    
    ereturn matrix b_placebo = `placebo_coefs_mat'
    
    if "`method'" == "IU" {
        ereturn matrix se_placebo = `placebo_se_mat'
        ereturn matrix V_placebo = `vcov_placebo_mat'
        if `threshold_specified' == 1 {
            ereturn matrix crit_values = `crit_values_mat'
        }
        else {
            ereturn matrix min_thresholds = `min_thresholds_mat'
        }
    }
    
    ereturn scalar max_abs_coef = `max_abs_coef'
    ereturn scalar alpha = `alpha'
    ereturn scalar N_g = `dp_n'
    ereturn scalar T = `dp_n_t'
    ereturn scalar N = `dp_N_obs'
    ereturn scalar no_placebos = `effective_no_placebos'
    ereturn scalar base_period = `dp_baseperiod'
    ereturn scalar is_balanced = `dp_balanced'
    
    // Store T_min and T_max for unbalanced panels
    if `dp_balanced' == 0 {
        ereturn scalar T_min = `dp_t_min'
        ereturn scalar T_max = `dp_t_max'
    }
    
    if `threshold_specified' == 1 {
        ereturn scalar threshold = `threshold_val'
        ereturn scalar critical_value = `critical_value'
        ereturn scalar reject = `reject'
    }
    else {
        ereturn scalar min_threshold = `min_threshold'
    }
    
    if "`method'" == "Boot" | "`method'" == "Wild" {
        ereturn scalar B = `reps'
        if `seed' > 0 {
            ereturn scalar seed = `seed'
        }
        else {
            ereturn scalar seed = .
        }
        ereturn scalar sigma_hathat_c = `sigma_hathat_c'
    }
    
    ereturn local method "`method'"
    ereturn local cmd "maxequivtest"
    ereturn local cmdline "`cmdline'"
    ereturn local depvar "`Y'"
    ereturn local idvar "`id'"
    ereturn local groupvar "`g'"
    ereturn local timevar "`period'"
    
    if "`cluster'" != "" {
        ereturn local clustvar "`cluster'"
    }
    if "`vce_type'" != "" {
        ereturn local vce "`vce_type'"
    }
    
    ereturn local placebo_names "`effective_placebo_names'"
    ereturn scalar threshold_specified = `threshold_specified'
    
end
