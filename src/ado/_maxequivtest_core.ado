*! _maxequivtest_core.ado - Core computation engine for the max equivalence test
*!
*! Implements the test statistic and critical value computation for hypothesis (3.1):
*!   H0: ||beta||_inf >= delta  vs  H1: ||beta||_inf < delta
*! Supports both the IU (intersection-union) method using the folded normal distribution
*! and bootstrap-based methods (standard and wild cluster bootstrap).

program define _maxequivtest_core, rclass
    version 16.0
    
    syntax, Y(varname numeric) ID(varname numeric) G(varname numeric) Period(varname numeric) ///
        NOplacebos(integer) ///
        THRESHOLDspecified(integer) THRESHold(real) ///
        Alpha(real) Method(string) Reps(integer) Seed(integer) ///
        [X(varlist numeric)] [Cluster(varname numeric)] [VCE(string)] ///
        [PREtreatment(numlist)] [BASEperiod(real -999999)] [NODOTS]
    
    * Store syntax arguments in standardized local variables
    local no_placebos = `noplacebos'
    local threshold_specified = `thresholdspecified'
    
    * ========================================================================
    * Progress display control
    * Set scalar to control bootstrap progress display in Mata
    * ========================================================================
    
    if "`nodots'" == "" {
        * Show progress by default for bootstrap methods
        scalar __eqt_show_dots = 1
    }
    else {
        * User requested no progress display
        scalar __eqt_show_dots = 0
    }
    
    * ========================================================================
    * Pretreatment period specification
    * ========================================================================
    
    * Construct pretreatment period matrix for Mata functions
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
        * Use sentinel value -999999 to indicate all pre-base periods
        * Note: Stata does not support zero-row matrices in J()
        matrix `pretreat_mat' = J(1, 1, -999999)
    }
    
    * ========================================================================
    * Test statistic computation
    * ========================================================================
    
    * Prepare covariate specification
    local x_vars "`x'"
    
    if "`method'" == "IU" {
        * Intersection-union test using folded normal quantiles (Section 4.2.1)
        * Test rejects H0 when |beta_t| < Q_N_F(threshold, Sigma_tt/n)(alpha) for all t
        capture noisily mata: _maxequivtest_core_iu_mata("`y'", "`x_vars'", "`id'", "`g'", "`period'", ///
            `no_placebos', `threshold_specified', `threshold', `alpha', "`vce'", "`cluster'", ///
            st_matrix("`pretreat_mat'"), `baseperiod')
        
        * Verify successful execution
        if _rc {
            display as error "maxequivtest: Mata computation failed (IU method)"
            display as error "  Error code: " _rc
            exit _rc
        }
        
        * Validate returned matrices
        capture confirm matrix r(placebo_coefs)
        if _rc {
            display as error "maxequivtest: Mata function did not return placebo_coefs matrix"
            exit 498
        }
        
        * Validate test statistic
        if missing(r(max_abs_coef)) {
            display as error "maxequivtest: Mata function returned invalid max_abs_coef"
            display as error "  Possible causes:"
            display as error "    - Missing values in outcome or treatment variables"
            display as error "    - Insufficient variation in treatment assignment"
            display as error "    - Singular variance-covariance matrix"
            exit 498
        }
        
        local max_abs_coef = r(max_abs_coef)
        local critical_value = r(critical_value)
        local reject = r(reject)
        local min_threshold = r(min_delta)
        local sigma_hathat_c = r(sigma_hathat_c)
        
        * Store scalar results before matrix assignments clear r()
        local effective_no_placebos = r(no_placebos)
        
        matrix placebo_coefs = r(placebo_coefs)
        matrix placebo_se = r(placebo_se)
        matrix vcov_placebo = r(vcov_placebo)
        
        * Build placebo coefficient names from effective periods
        tempname placebo_periods_mat
        matrix `placebo_periods_mat' = r(placebo_periods)
        local effective_placebo_names ""
        forvalues i = 1/`effective_no_placebos' {
            local period_val = `placebo_periods_mat'[`i', 1]
            local effective_placebo_names "`effective_placebo_names' placebo_`period_val'"
        }
        local effective_placebo_names = strtrim("`effective_placebo_names'")
        
        * Retrieve method-specific results for IU test
        if `threshold_specified' == 1 {
            matrix crit_values = r(crit_values)
        }
        else {
            matrix min_thresholds = r(min_thresholds)
        }
    }
    else {
        * Bootstrap-based test using constrained estimation (Section 4.2.1, equation 4.7)
        * Test rejects H0 when ||beta||_inf < Q_alpha^* (bootstrap quantile)
        capture noisily mata: _maxequivtest_core_bootstrap("`y'", "`x_vars'", "`id'", "`g'", "`period'", ///
            `no_placebos', `threshold_specified', `threshold', `alpha', `reps', ///
            "`method'", `seed', st_matrix("`pretreat_mat'"), `baseperiod')
        
        * Verify successful execution
        if _rc {
            display as error "maxequivtest: Mata computation failed (`method' method)"
            display as error "  Error code: " _rc
            exit _rc
        }
        
        * Validate returned matrices
        capture confirm matrix r(placebo_coefs)
        if _rc {
            display as error "maxequivtest: Mata function did not return placebo_coefs matrix"
            exit 498
        }
        
        * Validate test statistic
        if missing(r(max_abs_coef)) {
            display as error "maxequivtest: Mata function returned invalid max_abs_coef"
            display as error "  Possible causes:"
            display as error "    - Missing values in outcome or treatment variables"
            display as error "    - Constrained optimization failed to converge"
            display as error "    - Insufficient bootstrap replications"
            exit 498
        }
        
        local max_abs_coef = r(max_abs_coef)
        local critical_value = r(critical_value)
        local reject = r(reject)
        local min_threshold = r(min_delta)
        local sigma_hathat_c = r(sigma_hathat_c)
        
        * Store scalar results before matrix assignments clear r()
        local effective_no_placebos = r(no_placebos)
        
        matrix placebo_coefs = r(placebo_coefs)
        
        * Build placebo coefficient names from effective periods
        tempname placebo_periods_mat
        matrix `placebo_periods_mat' = r(placebo_periods)
        local effective_placebo_names ""
        forvalues i = 1/`effective_no_placebos' {
            local period_val = `placebo_periods_mat'[`i', 1]
            local effective_placebo_names "`effective_placebo_names' placebo_`period_val'"
        }
        local effective_placebo_names = strtrim("`effective_placebo_names'")
    }
    
    * ========================================================================
    * Return test results to calling program
    * ========================================================================
    
    * Core test results
    return matrix placebo_coefs = placebo_coefs
    return scalar max_abs_coef = `max_abs_coef'
    return scalar critical_value = `critical_value'
    return scalar reject = `reject'
    return scalar min_threshold = `min_threshold'
    return scalar sigma_hathat_c = `sigma_hathat_c'
    
    * Effective sample characteristics after data processing
    return scalar no_placebos = `effective_no_placebos'
    return local placebo_names = "`effective_placebo_names'"
    
    * Additional results for IU method
    if "`method'" == "IU" {
        return matrix placebo_se = placebo_se
        return matrix vcov_placebo = vcov_placebo
        if `threshold_specified' == 1 {
            return matrix crit_values = crit_values
        }
        else {
            return matrix min_thresholds = min_thresholds
        }
    }
    
    * ========================================================================
    * Cleanup: remove progress display control scalar
    * ========================================================================
    
    capture scalar drop __eqt_show_dots
    
end
