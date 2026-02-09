*! _equivtest_ereturn.ado - Post-estimation result storage for equivalence tests
*!
*! Stores estimation results in e() for post-estimation access and replay.
*! Organizes results by test type following the equivalence testing framework
*! for assessing pre-trends in Difference-in-Differences estimation.
*!
*! Hypotheses stored:
*!   - Max test: H0: ||beta||_inf >= delta (maximum placebo coefficient bound)
*!   - Mean test: H0: |beta_bar| >= tau (average placebo coefficient bound)
*!   - RMS test: H0: beta_RMS >= zeta (root mean square bound)

program define _equivtest_ereturn, eclass
    version 16.0
    
    syntax, TYPE(string) CMDLINE(string) ///
        DEPVAR(string) IDVAR(string) GROUPVAR(string) TIMEVAR(string) ///
        Nobs(integer) N_g(integer) NO_placebos(integer) ///
        Alpha(real) BASEperiod(real) BALanced(integer) ///
        THRESHOLD_specified(integer) ///
        [PREperiods(string)] [PLACEBO_names(string)] ///
        [Method(string)] [THRESHold(real -999999)] [REJect(integer -999999)] ///
        [MIN_threshold(real -999999)] [VCE(string)] [CLUSTvar(string)] ///
        [Tmin(integer -999999)] [Tmax(integer -999999)] [Tcount(integer -999999)] ///
        [B_placebo(name)] [SE_placebo(name)] [V_placebo(name)] ///
        [IU_critical_values(name)] [MIN_equiv_thresholds(name)] ///
        [MAX_abs_coef(real -999999)] [BOOT_critical(real -999999)] ///
        [Nboot(integer 1000)] [WILD(integer 0)] ///
        [ABS_mean_placebo(real -999999)] [VAR_mean_placebo(real -999999)] ///
        [SE_mean_placebo(real -999999)] [P_value(real -999999)] ///
        [MEAN_critical_value(real -999999)] ///
        [RMS_placebo_coefs(real -999999)] [RMS_critical_value(real -999999)] ///
        [NOLambda(integer 5)]
    
    // -------------------------------------------------------------------------
    // Initialize e() class storage
    // -------------------------------------------------------------------------
    ereturn clear
    
    // -------------------------------------------------------------------------
    // Common scalars across all test types
    // -------------------------------------------------------------------------
    ereturn scalar N = `nobs'
    ereturn scalar N_g = `n_g'
    ereturn scalar no_placebos = `no_placebos'
    ereturn scalar alpha = `alpha'
    ereturn scalar base_period = `baseperiod'
    ereturn scalar is_balanced = `balanced'
    
    // Pre-treatment periods: T for balanced panels, [T_min, T_max] for unbalanced
    if `balanced' == 1 {
        if `tcount' != -999999 {
            ereturn scalar T_pre = `tcount'
        }
    }
    else {
        if `tmin' != -999999 {
            ereturn scalar T_min = `tmin'
        }
        if `tmax' != -999999 {
            ereturn scalar T_max = `tmax'
        }
    }
    
    // -------------------------------------------------------------------------
    // Equivalence threshold storage
    // When threshold is specified: stores threshold value and rejection indicator
    // When threshold is not specified: stores minimum threshold for rejection
    // -------------------------------------------------------------------------
    ereturn scalar threshold_specified = `threshold_specified'
    
    if `threshold_specified' == 1 {
        if `threshold' != -999999 {
            ereturn scalar threshold = `threshold'
        }
        if `reject' != -999999 {
            ereturn scalar reject = `reject'
        }
    }
    else {
        if `min_threshold' != -999999 {
            ereturn scalar min_threshold = `min_threshold'
        }
    }
    
    // -------------------------------------------------------------------------
    // Test-specific results storage
    // -------------------------------------------------------------------------
    
    if "`type'" == "max" {
        if "`method'" == "iu" | "`method'" == "" {
            // -----------------------------------------------------------------
            // Intersection-union (IU) test for H0: ||beta||_inf >= delta
            // Rejection rule: reject if |beta_t| < Q_{FN}(alpha) for all t
            // where Q_{FN} is the alpha-quantile of the folded normal
            // -----------------------------------------------------------------
            
            if "`b_placebo'" != "" {
                ereturn matrix b_placebo = `b_placebo'
            }
            
            if "`se_placebo'" != "" {
                ereturn matrix se_placebo = `se_placebo'
            }
            
            if "`v_placebo'" != "" {
                ereturn matrix V_placebo = `v_placebo'
            }
            
            if `max_abs_coef' != -999999 {
                ereturn scalar max_abs_coef = `max_abs_coef'
            }
            
            if `threshold_specified' == 1 {
                if "`iu_critical_values'" != "" {
                    ereturn matrix IU_critical_values = `iu_critical_values'
                }
            }
            else {
                if "`min_equiv_thresholds'" != "" {
                    ereturn matrix min_equiv_thresholds = `min_equiv_thresholds'
                }
            }
        }
        else {
            // -----------------------------------------------------------------
            // Bootstrap test for H0: ||beta||_inf >= delta
            // Rejection rule: reject if ||beta||_inf < Q_alpha^*
            // where Q_alpha^* is the alpha-quantile from bootstrap samples
            // -----------------------------------------------------------------
            
            // Store placebo coefficients for equivtest_plot compatibility
            if "`b_placebo'" != "" {
                ereturn matrix b_placebo = `b_placebo'
            }
            
            if `max_abs_coef' != -999999 {
                ereturn scalar max_abs_coef = `max_abs_coef'
            }
            
            ereturn scalar nboot = `nboot'
            ereturn scalar wild = `wild'
            
            if `threshold_specified' == 1 {
                if `boot_critical' != -999999 {
                    ereturn scalar boot_critical = `boot_critical'
                }
            }
            
            // Bootstrap inference uses empirical quantiles; standard errors not applicable
        }
    }
    else if "`type'" == "mean" {
        // ---------------------------------------------------------------------
        // Mean test for H0: |beta_bar| >= tau
        // Test statistic: beta_bar = (1/T) * sum(beta_t)
        // Rejection rule: reject if |beta_bar| < Q_{FN}(tau, sigma^2)(alpha)
        // ---------------------------------------------------------------------
        
        if "`b_placebo'" != "" {
            ereturn matrix b_placebo = `b_placebo'
        }
        
        if "`v_placebo'" != "" {
            ereturn matrix V_placebo = `v_placebo'
        }
        
        if `abs_mean_placebo' != -999999 {
            ereturn scalar abs_mean_placebo = `abs_mean_placebo'
        }
        if `var_mean_placebo' != -999999 {
            ereturn scalar var_mean_placebo = `var_mean_placebo'
        }
        if `se_mean_placebo' != -999999 {
            ereturn scalar se_mean_placebo = `se_mean_placebo'
        }
        
        if `threshold_specified' == 1 {
            if `p_value' != -999999 {
                ereturn scalar p_value = `p_value'
            }
            if `mean_critical_value' != -999999 {
                ereturn scalar mean_critical_value = `mean_critical_value'
            }
        }
    }
    else if "`type'" == "rms" {
        // ---------------------------------------------------------------------
        // RMS (root mean square) test for H0: beta_RMS >= zeta
        // Test statistic: beta_RMS = sqrt((1/T) * sum(beta_t^2))
        // Uses self-normalized inference with subsample variance estimation
        // ---------------------------------------------------------------------
        
        if "`b_placebo'" != "" {
            // Ensure b_placebo is stored as column vector for consistency
            // with other test types (required for equivtest_plot)
            tempname b_placebo_col
            if colsof(`b_placebo') > rowsof(`b_placebo') {
                matrix `b_placebo_col' = `b_placebo''
            }
            else {
                matrix `b_placebo_col' = `b_placebo'
            }
            ereturn matrix b_placebo = `b_placebo_col'
        }
        
        if `rms_placebo_coefs' != -999999 {
            ereturn scalar rms_placebo_coefs = `rms_placebo_coefs'
        }
        
        ereturn scalar nolambda = `nolambda'
        
        if `threshold_specified' == 1 {
            if `rms_critical_value' != -999999 {
                ereturn scalar rms_critical_value = `rms_critical_value'
            }
        }
    }
    
    // -------------------------------------------------------------------------
    // Common macros across all test types
    // -------------------------------------------------------------------------
    ereturn local cmd "equivtest"
    ereturn local cmdline `"`cmdline'"'
    ereturn local type "`type'"
    
    if "`method'" != "" {
        ereturn local method "`method'"
    }
    
    ereturn local depvar "`depvar'"
    ereturn local idvar "`idvar'"
    ereturn local groupvar "`groupvar'"
    ereturn local timevar "`timevar'"
    
    if "`vce'" != "" {
        ereturn local vce "`vce'"
    }
    if "`clustvar'" != "" {
        ereturn local clustvar "`clustvar'"
    }
    
    ereturn local placebo_coef_names "`placebo_names'"
    ereturn local preperiods "`preperiods'"
    ereturn local properties "b"
    
end
