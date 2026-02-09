*! equivtest.ado - Unified equivalence testing for pre-trends in DiD estimation
*!
*! Implements equivalence tests for assessing the plausibility of the parallel
*! trends assumption (PTA) in Difference-in-Differences estimation. Three
*! distinct hypotheses are supported, each imposing bounds on placebo treatment
*! effects in pre-treatment periods:
*!
*!   type(max)  : H_0: ||beta||_inf >= delta  vs  H_1: ||beta||_inf < delta
*!   type(mean) : H_0: |beta_bar| >= tau      vs  H_1: |beta_bar| < tau
*!   type(rms)  : H_0: beta_RMS >= zeta       vs  H_1: beta_RMS < zeta
*!
*! Rejection of the null hypothesis provides statistical evidence in favor of
*! negligible trend differences, increasing the credibility of the PTA.

program define equivtest, eclass
    version 16.0
    
    // Initialize Mata library
    equitrends_init
    
    // Store full command line for result replication
    local cmdline "equivtest `0'"
    
    // Parse command options
    // Note: Group() accepts g/gr/gro/grou/group; Time() accepts t/ti/tim/time
    // Period() is an alias for Time() for backward compatibility
    syntax varname(numeric), TYPE(string) ID(varname numeric) ///
        [Group(varname numeric) G(varname numeric)] ///
        [Time(varname numeric) Period(varname numeric)] ///
        [X(varlist numeric)] [THRESHold(real -999999)] [PREtreatment(numlist)] ///
        [BASEperiod(real -999999)] [Method(string)] [VCE(string)] ///
        [Alpha(real 0.05)] [Nboot(integer 1000)] [NOLambda(integer 5)] ///
        [SEED(integer -999999)] [Cluster(varname numeric)] [NODOTS]
    
    local depvar `varlist'
    
    // -------------------------------------------------------------------------
    // Handle Option Aliases
    // -------------------------------------------------------------------------
    
    // Group/G alias handling: accept both group() and g()
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
    
    // Time/Period alias handling: accept both time() and period()
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
    // Input Validation
    // -------------------------------------------------------------------------
    
    // Construct parameter validation command with panel structure identifiers
    local parse_cmd "_equivtest_parse, type(`type') period(`period') g(`g')"
    
    if "`method'" != "" {
        local parse_cmd "`parse_cmd' method(`method')"
    }
    if "`vce'" != "" {
        local parse_cmd "`parse_cmd' vce(`vce')"
    }
    local parse_cmd "`parse_cmd' alpha(`alpha')"
    local parse_cmd "`parse_cmd' nboot(`nboot')"
    local parse_cmd "`parse_cmd' nolambda(`nolambda')"
    if "`cluster'" != "" {
        local parse_cmd "`parse_cmd' cluster(`cluster')"
    }
    if `threshold' != -999999 {
        local parse_cmd "`parse_cmd' threshold(`threshold')"
    }
    if "`pretreatment'" != "" {
        local parse_cmd "`parse_cmd' pretreatment(`pretreatment')"
    }
    if `baseperiod' != -999999 {
        local parse_cmd "`parse_cmd' baseperiod(`baseperiod')"
    }
    
    // Validate parameters and resolve defaults
    `parse_cmd'
    
    // Extract validated parameter values
    local type_val "`r(type)'"
    local method_val "`r(method)'"
    local vce_val "`r(vce)'"
    local alpha_val = r(alpha)
    local nboot_val = r(nboot)
    local nolambda_val = r(nolambda)
    local threshold_specified = r(threshold_specified)
    local threshold_val = r(threshold)
    
    // -------------------------------------------------------------------------
    // Route to Test Implementation
    // -------------------------------------------------------------------------
    
    if "`type_val'" == "max" {
        // Maximum test for hypothesis: H_0: ||beta||_inf >= delta
        // Bounds the largest absolute placebo coefficient across pre-treatment periods
        local subcmd "maxequivtest `depvar', id(`id') g(`g') period(`period') nodisplay"
        
        // Append optional parameters to subcommand
        if "`x'" != "" {
            local subcmd "`subcmd' x(`x')"
        }
        if `threshold_specified' == 1 {
            local subcmd "`subcmd' threshold(`threshold_val')"
        }
        if "`pretreatment'" != "" {
            local subcmd "`subcmd' pretreatment(`pretreatment')"
        }
        if `baseperiod' != -999999 {
            local subcmd "`subcmd' baseperiod(`baseperiod')"
        }
        local subcmd "`subcmd' alpha(`alpha_val')"
        
        // Method-specific inference options
        if "`method_val'" == "iu" {
            // Intersection-union (IU) test: rejects when all individual tests reject
            // Test statistic follows folded normal distribution under H_0
            local subcmd "`subcmd' method(IU)"
            if "`vce_val'" != "" {
                local subcmd "`subcmd' vce(`vce_val')"
            }
            if "`cluster'" != "" {
                local subcmd "`subcmd' cluster(`cluster')"
            }
        }
        else if "`method_val'" == "boot" {
            // Parametric bootstrap: generates samples under ||beta||_inf = delta
            // More powerful than IU but assumes spherical errors
            local subcmd "`subcmd' method(Boot)"
            local subcmd "`subcmd' reps(`nboot_val')"
            if `seed' != -999999 {
                local subcmd "`subcmd' seed(`seed')"
            }
        }
        else if "`method_val'" == "wild" {
            // Wild cluster bootstrap: robust to heteroskedasticity and serial correlation
            // Uses Rademacher weights to preserve within-cluster dependence structure
            local subcmd "`subcmd' method(Wild)"
            local subcmd "`subcmd' reps(`nboot_val')"
            if `seed' != -999999 {
                local subcmd "`subcmd' seed(`seed')"
            }
        }
        
        // Pass nodots option to suppress bootstrap progress display
        if "`nodots'" != "" {
            local subcmd "`subcmd' nodots"
        }
        
        // For bootstrap methods, allow progress display unless nodots is specified
        // Don't use quietly so progress bar can be displayed
        // Note: nodisplay option already suppresses the final result output
        if ("`method_val'" == "boot" | "`method_val'" == "wild") & "`nodots'" == "" {
            `subcmd'
        }
        else {
            quietly `subcmd'
        }
        
        // Extract test results and panel structure
        local max_abs_coef = e(max_abs_coef)
        local n_obs = e(N)
        local n_g = e(N_g)
        local n_t = e(T)
        local no_placebos = e(no_placebos)
        
        // Retrieve placebo coefficient estimates: beta_l = gamma_l - gamma_{T+1}
        tempname b_placebo se_placebo V_placebo iu_crit_values iu_min_thresholds
        matrix `b_placebo' = e(b_placebo)
        
        // IU method provides standard errors and variance matrix
        if "`method_val'" == "iu" {
            capture matrix `se_placebo' = e(se_placebo)
            capture matrix `V_placebo' = e(V_placebo)
            if `threshold_specified' == 1 {
                capture matrix `iu_crit_values' = e(crit_values)
            }
            else {
                capture matrix `iu_min_thresholds' = e(min_thresholds)
            }
        }
        
        if `threshold_specified' == 1 {
            local reject = e(reject)
            local critical_value = e(critical_value)
        }
        else {
            local min_threshold = e(min_threshold)
        }
        
        // Extract panel structure metadata
        local baseperiod_val = e(base_period)
        local is_balanced = e(is_balanced)
        if `is_balanced' == . {
            local is_balanced = 1
        }
        
        local placebo_names_val "`e(placebo_names)'"
    }
    else if "`type_val'" == "mean" {
        // Mean test for hypothesis: H_0: |beta_bar| >= tau
        // beta_bar = (1/T) * sum(beta_l) is the average placebo coefficient
        // Sensitive to cancellation effects when coefficients have opposing signs
        local subcmd "meanequivtest `depvar', id(`id') g(`g') period(`period') nodisplay"
        
        // Append optional parameters to subcommand
        if "`x'" != "" {
            local subcmd "`subcmd' x(`x')"
        }
        if `threshold_specified' == 1 {
            local subcmd "`subcmd' threshold(`threshold_val')"
        }
        if "`pretreatment'" != "" {
            local subcmd "`subcmd' pretreatment(`pretreatment')"
        }
        if `baseperiod' != -999999 {
            local subcmd "`subcmd' baseperiod(`baseperiod')"
        }
        local subcmd "`subcmd' alpha(`alpha_val')"
        
        // Variance estimation specification for standard error computation
        if "`vce_val'" != "" {
            local subcmd "`subcmd' vce(`vce_val')"
        }
        if "`cluster'" != "" {
            local subcmd "`subcmd' cluster(`cluster')"
        }
        
        quietly `subcmd'
        
        // Extract test statistic: sqrt(n) * 1'(beta_hat - beta) -> N(0, 1'*Sigma*1)
        local abs_mean_placebo = e(abs_mean_placebo)
        local var_mean_placebo = e(var_mean_placebo)
        local se_mean_placebo = e(se_mean_placebo)
        local n_obs = e(N)
        local n_g = e(N_g)
        local n_t = e(T)
        local no_placebos = e(no_placebos)
        local baseperiod_val = e(base_period)
        local is_balanced = e(is_balanced)
        
        // Retrieve coefficient vector and variance-covariance matrix
        tempname b_placebo V_placebo
        matrix `b_placebo' = e(b_placebo)
        matrix `V_placebo' = e(V_placebo)
        
        if `threshold_specified' == 1 {
            local reject = e(reject)
            local p_value = e(p_value)
            local mean_critical = e(critical_value)
        }
        else {
            local min_threshold = e(min_threshold)
        }
        
        local placebo_names_val "`e(placebo_names)'"
    }
    else if "`type_val'" == "rms" {
        // RMS test for hypothesis: H_0: beta_RMS >= zeta
        // beta_RMS = ||beta||/sqrt(T) is the root mean square of placebo coefficients
        // Self-normalized test statistic: M_n = (beta_RMS^2 - zeta^2) / V_n
        // V_n is computed via subsampling without requiring variance estimation
        local subcmd "rmsequivtest `depvar', id(`id') g(`g') period(`period') nodisplay"
        
        // Append optional parameters
        if "`x'" != "" {
            local subcmd "`subcmd' x(`x')"
        }
        if `threshold_specified' == 1 {
            local subcmd "`subcmd' threshold(`threshold_val')"
        }
        if "`pretreatment'" != "" {
            local subcmd "`subcmd' pretreatment(`pretreatment')"
        }
        if `baseperiod' != -999999 {
            local subcmd "`subcmd' baseperiod(`baseperiod')"
        }
        local subcmd "`subcmd' alpha(`alpha_val')"
        // Number of subsample fractions lambda_k for variance computation
        // V_n = sqrt(sum_k (beta_RMS^2(lambda_k) - beta_RMS^2(1))^2 / K)
        local subcmd "`subcmd' nolambda(`nolambda_val')"
        
        if `seed' != -999999 {
            local subcmd "`subcmd' seed(`seed')"
        }
        
        quietly `subcmd'
        
        // Extract RMS statistic and panel structure
        local rms_placebo = e(RMS_placebo)
        local n_obs = e(N)
        local n_g = e(N_g)
        local n_t = e(T)
        local no_placebos = e(no_placebos)
        local baseperiod_val = e(base_period)
        local is_balanced = e(is_balanced)
        
        // Retrieve placebo coefficient vector for reference
        tempname b_placebo
        matrix `b_placebo' = e(b_placebo)
        
        if `threshold_specified' == 1 {
            local reject = e(reject)
            local rms_critical = e(RMS_critical)
        }
        else {
            local min_threshold = e(min_threshold)
        }
        
        local placebo_names_val "`e(placebo_names)'"
    }

    // -------------------------------------------------------------------------
    // Output Display
    // -------------------------------------------------------------------------
    
    // Assemble display command with test-specific parameters
    local disp_cmd "_equivtest_display, type(`type_val')"
    local disp_cmd "`disp_cmd' threshold_specified(`threshold_specified')"
    local disp_cmd "`disp_cmd' alpha(`alpha_val')"
    local disp_cmd "`disp_cmd' n(`n_obs') n_g(`n_g') no_placebos(`no_placebos')"
    local disp_cmd "`disp_cmd' baseperiod(`baseperiod_val') balanced(`is_balanced')"
    
    // Include variance estimation method in output
    if "`vce_val'" != "" {
        local disp_cmd "`disp_cmd' vce(`vce_val')"
    }
    
    // Handle balanced vs unbalanced panel time structure
    if `is_balanced' == 1 {
        local disp_cmd "`disp_cmd' tcount(`n_t')"
    }
    else {
        // Unbalanced panels: report minimum and maximum period counts
        capture local t_min = e(T_min)
        capture local t_max = e(T_max)
        if "`t_min'" != "" & "`t_max'" != "" {
            local disp_cmd "`disp_cmd' tmin(`t_min') tmax(`t_max')"
        }
        else {
            local disp_cmd "`disp_cmd' tcount(`n_t')"
        }
    }
    
    // Threshold specification determines output format:
    // - When specified: display rejection decision at given equivalence threshold
    // - When omitted: report minimum threshold for rejection (equivalence confidence bound)
    if `threshold_specified' == 1 {
        local disp_cmd "`disp_cmd' threshold(`threshold_val') reject(`reject')"
    }
    else {
        local disp_cmd "`disp_cmd' min_threshold(`min_threshold')"
    }
    
    // Test-specific display parameters
    if "`type_val'" == "max" {
        local disp_cmd "`disp_cmd' method(`method_val')"
        local disp_cmd "`disp_cmd' max_abs_coef(`max_abs_coef')"
        
        if "`method_val'" == "iu" {
            // Format coefficient-level results for tabular display
            local placebo_coefs_str ""
            local se_placebo_str ""
            local cv_or_mt_str ""
            
            forvalues i = 1/`no_placebos' {
                local coef_val = `b_placebo'[`i', 1]
                local placebo_coefs_str "`placebo_coefs_str' `coef_val'"
                
                capture local se_val = `se_placebo'[`i', 1]
                if _rc == 0 {
                    local se_placebo_str "`se_placebo_str' `se_val'"
                }
                
                if `threshold_specified' == 1 {
                    capture local cv_val = `iu_crit_values'[`i', 1]
                    if _rc == 0 {
                        local cv_or_mt_str "`cv_or_mt_str' `cv_val'"
                    }
                }
                else {
                    capture local mt_val = `iu_min_thresholds'[`i', 1]
                    if _rc == 0 {
                        local cv_or_mt_str "`cv_or_mt_str' `mt_val'"
                    }
                }
            }
            
            local placebo_coefs_str = strtrim("`placebo_coefs_str'")
            local se_placebo_str = strtrim("`se_placebo_str'")
            local cv_or_mt_str = strtrim("`cv_or_mt_str'")
            
            local disp_cmd "`disp_cmd' placebo_coefs(`placebo_coefs_str')"
            local disp_cmd "`disp_cmd' se_placebo(`se_placebo_str')"
            
            if `threshold_specified' == 1 {
                local disp_cmd "`disp_cmd' critical_values(`cv_or_mt_str')"
            }
            else {
                local disp_cmd "`disp_cmd' min_thresholds(`cv_or_mt_str')"
            }
        }
        else {
            // Bootstrap methods: report replication count and empirical critical value
            local disp_cmd "`disp_cmd' nboot(`nboot_val')"
            if `threshold_specified' == 1 {
                local disp_cmd "`disp_cmd' boot_critical(`critical_value')"
            }
        }
    }
    else if "`type_val'" == "mean" {
        local disp_cmd "`disp_cmd' abs_mean_placebo(`abs_mean_placebo')"
        local disp_cmd "`disp_cmd' se_mean_placebo(`se_mean_placebo')"
        if `threshold_specified' == 1 {
            local disp_cmd "`disp_cmd' p_value(`p_value')"
            local disp_cmd "`disp_cmd' mean_critical(`mean_critical')"
        }
    }
    else if "`type_val'" == "rms" {
        local disp_cmd "`disp_cmd' rms_placebo(`rms_placebo')"
        if `threshold_specified' == 1 {
            local disp_cmd "`disp_cmd' rms_critical(`rms_critical')"
        }
        if `seed' != -999999 {
            local disp_cmd "`disp_cmd' seed(`seed')"
        }
    }
    
    // Include period labels for coefficient table
    if "`placebo_names_val'" != "" {
        local disp_cmd "`disp_cmd' placebo_names(`placebo_names_val')"
    }
    
    `disp_cmd'
    
    // -------------------------------------------------------------------------
    // Store Estimation Results
    // -------------------------------------------------------------------------
    
    // Construct ereturn command for post-estimation access
    local eret_cmd "_equivtest_ereturn, type(`type_val')"
    local eret_cmd "`eret_cmd' cmdline(`"`cmdline'"')"
    local eret_cmd "`eret_cmd' depvar(`depvar') idvar(`id') groupvar(`g') timevar(`period')"
    local eret_cmd "`eret_cmd' nobs(`n_obs') n_g(`n_g') no_placebos(`no_placebos')"
    local eret_cmd "`eret_cmd' alpha(`alpha_val') baseperiod(`baseperiod_val') balanced(`is_balanced')"
    local eret_cmd "`eret_cmd' threshold_specified(`threshold_specified')"
    local eret_cmd "`eret_cmd' preperiods(`"`pretreatment'"') placebo_names(`"`placebo_names_val'"')"
    
    // Store panel time dimension
    if `is_balanced' == 1 {
        local eret_cmd "`eret_cmd' tcount(`n_t')"
    }
    else {
        if "`t_min'" != "" & "`t_max'" != "" {
            local eret_cmd "`eret_cmd' tmin(`t_min') tmax(`t_max')"
        }
    }
    
    // Store test decision or minimum equivalence threshold
    if `threshold_specified' == 1 {
        local eret_cmd "`eret_cmd' threshold(`threshold_val') reject(`reject')"
    }
    else {
        local eret_cmd "`eret_cmd' min_threshold(`min_threshold')"
    }
    
    // Store variance estimation specification
    if "`vce_val'" != "" {
        local eret_cmd "`eret_cmd' vce(`vce_val')"
    }
    if "`cluster'" != "" {
        local eret_cmd "`eret_cmd' clustvar(`cluster')"
    }
    
    // Store test-specific estimation results
    if "`type_val'" == "max" {
        local eret_cmd "`eret_cmd' method(`method_val')"
        local eret_cmd "`eret_cmd' max_abs_coef(`max_abs_coef')"
        local eret_cmd "`eret_cmd' b_placebo(`b_placebo')"
        
        if "`method_val'" == "iu" {
            // IU method: coefficient-level inference quantities
            local eret_cmd "`eret_cmd' se_placebo(`se_placebo')"
            local eret_cmd "`eret_cmd' v_placebo(`V_placebo')"
            if `threshold_specified' == 1 {
                local eret_cmd "`eret_cmd' iu_critical_values(`iu_crit_values')"
            }
            else {
                local eret_cmd "`eret_cmd' min_equiv_thresholds(`iu_min_thresholds')"
            }
        }
        else {
            // Bootstrap methods: store replication count and bootstrap critical value
            local eret_cmd "`eret_cmd' nboot(`nboot_val')"
            if "`method_val'" == "wild" {
                local eret_cmd "`eret_cmd' wild(1)"
            }
            else {
                local eret_cmd "`eret_cmd' wild(0)"
            }
            if `threshold_specified' == 1 {
                local eret_cmd "`eret_cmd' boot_critical(`critical_value')"
            }
        }
    }
    else if "`type_val'" == "mean" {
        local eret_cmd "`eret_cmd' b_placebo(`b_placebo') v_placebo(`V_placebo')"
        local eret_cmd "`eret_cmd' abs_mean_placebo(`abs_mean_placebo')"
        local eret_cmd "`eret_cmd' var_mean_placebo(`var_mean_placebo')"
        local eret_cmd "`eret_cmd' se_mean_placebo(`se_mean_placebo')"
        if `threshold_specified' == 1 {
            local eret_cmd "`eret_cmd' p_value(`p_value')"
            local eret_cmd "`eret_cmd' mean_critical_value(`mean_critical')"
        }
    }
    else if "`type_val'" == "rms" {
        local eret_cmd "`eret_cmd' b_placebo(`b_placebo')"
        local eret_cmd "`eret_cmd' rms_placebo_coefs(`rms_placebo')"
        local eret_cmd "`eret_cmd' nolambda(`nolambda_val')"
        if `threshold_specified' == 1 {
            local eret_cmd "`eret_cmd' rms_critical_value(`rms_critical')"
        }
    }
    
    `eret_cmd'
    
end
