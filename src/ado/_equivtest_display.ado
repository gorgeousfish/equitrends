*! _equivtest_display.ado - Unified output display for equivalence tests
*! Version 2.0 - Professional Stata-style output with SMCL formatting
*!
*! Dispatches result display to type-specific subroutines based on test type.
*! Supports maximum (IU and bootstrap), mean, and RMS tests for equivalence
*! of pre-trends in Difference-in-Differences estimation.
*!
*! Output structure:
*!   1. Header with title and test configuration
*!   2. Type-specific test results
*!   3. Panel information footer

program define _equivtest_display
    version 16.0
    
    // =========================================================================
    // Syntax definition
    // =========================================================================
    
    syntax, TYPE(string) ///
        THRESHOLD_specified(integer) Alpha(real) ///
        N(integer) N_g(integer) NO_placebos(integer) ///
        BASEperiod(real) BALanced(integer) ///
        [Method(string)] [THRESHold(real -999999)] [REJect(integer -999999)] ///
        [MIN_threshold(real -999999)] [Nboot(integer 1000)] ///
        [Tmin(integer -999999)] [Tmax(integer -999999)] [Tcount(integer -999999)] ///
        [PLACEBO_coefs(string)] [SE_placebo(string)] ///
        [CRITICAL_values(string)] [MIN_thresholds(string)] ///
        [MAX_abs_coef(real -999999)] [BOOT_critical(real -999999)] ///
        [ABS_mean_placebo(real -999999)] [SE_mean_placebo(real -999999)] ///
        [P_value(real -999999)] [MEAN_critical(real -999999)] ///
        [RMS_placebo(real -999999)] [RMS_critical(real -999999)] ///
        [PLACEBO_names(string)] [VCE(string)] [SEED(integer -999999)]
    
    // =========================================================================
    // Output Header with Test Configuration
    // =========================================================================
    
    // Build title command with configuration parameters
    local title_cmd "_equivtest_display_title, type(`type')"
    local title_cmd "`title_cmd' alpha(`alpha')"
    local title_cmd "`title_cmd' threshold_specified(`threshold_specified')"
    
    if "`method'" != "" {
        local title_cmd "`title_cmd' method(`method')"
    }
    if `threshold_specified' == 1 & `threshold' != -999999 {
        local title_cmd "`title_cmd' threshold(`threshold')"
    }
    if "`vce'" != "" {
        local title_cmd "`title_cmd' vce(`vce')"
    }
    if "`method'" == "boot" | "`method'" == "wild" {
        local title_cmd "`title_cmd' nboot(`nboot')"
    }
    
    `title_cmd'
    
    // =========================================================================
    // Type-specific Result Display
    // =========================================================================
    
    if "`type'" == "max" {
        // Maximum equivalence test: H0: ||beta||_inf >= delta
        
        if "`method'" == "iu" | "`method'" == "" {
            // Intersection-union (IU) test based on folded normal quantiles
            _equivtest_display_iu, ///
                threshold_specified(`threshold_specified') ///
                threshold(`threshold') alpha(`alpha') reject(`reject') ///
                min_threshold(`min_threshold') ///
                placebo_coefs(`placebo_coefs') se_placebo(`se_placebo') ///
                critical_values(`critical_values') min_thresholds(`min_thresholds') ///
                placebo_names(`placebo_names')
        }
        else {
            // Bootstrap-based test (parametric or wild cluster bootstrap)
            _equivtest_display_boot, ///
                method(`method') threshold_specified(`threshold_specified') ///
                threshold(`threshold') alpha(`alpha') reject(`reject') ///
                min_threshold(`min_threshold') nboot(`nboot') ///
                max_abs_coef(`max_abs_coef') boot_critical(`boot_critical')
        }
    }
    else if "`type'" == "mean" {
        // Mean equivalence test: H0: |mean(beta)| >= tau
        _equivtest_display_mean, ///
            threshold_specified(`threshold_specified') ///
            threshold(`threshold') alpha(`alpha') reject(`reject') ///
            min_threshold(`min_threshold') ///
            abs_mean_placebo(`abs_mean_placebo') se_mean_placebo(`se_mean_placebo') ///
            p_value(`p_value') mean_critical(`mean_critical')
    }
    else if "`type'" == "rms" {
        // RMS equivalence test: H0: beta_RMS >= zeta (self-normalized)
        local rms_disp_cmd "_equivtest_display_rms, threshold_specified(`threshold_specified')"
        local rms_disp_cmd "`rms_disp_cmd' threshold(`threshold') alpha(`alpha')"
        local rms_disp_cmd "`rms_disp_cmd' reject(`reject') min_threshold(`min_threshold')"
        local rms_disp_cmd "`rms_disp_cmd' rms_placebo(`rms_placebo') rms_critical(`rms_critical')"
        if `seed' != -999999 {
            local rms_disp_cmd "`rms_disp_cmd' seed(`seed')"
        }
        `rms_disp_cmd'
    }
    
    // =========================================================================
    // Panel Information Footer
    // =========================================================================
    
    _equivtest_display_panel, type(`type') method(`method') ///
        no_placebos(`no_placebos') baseperiod(`baseperiod') ///
        balanced(`balanced') tmin(`tmin') tmax(`tmax') tcount(`tcount') ///
        n_g(`n_g') nobs(`n')
    
end
