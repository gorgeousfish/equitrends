*! _maxequivtest_display.ado - Display results for maximum equivalence test
*! Version 2.0 - Professional Stata-style output with SMCL formatting
*!
*! This module formats and displays inference results for testing the hypothesis
*! H0: ||beta||_inf >= delta vs H1: ||beta||_inf < delta, where ||.||_inf denotes
*! the supremum norm. Output includes test configuration, placebo coefficients,
*! and test decisions in a unified Stata-style format.

program define _maxequivtest_display
    version 16.0
    
    syntax, Method(string) ///
        THRESHOLD_specified(integer) THRESHold(real) ///
        Alpha(real) Reps(integer) ///
        MAX_abs_coef(real) CRITICAL_value(real) ///
        Reject(real) MIN_threshold(real) ///
        N(integer) N_t(integer) N_total(integer) ///
        NO_placebos(integer) [PLACEBO_names(string)] ///
        [PLACEBO_coefs(string)] ///
        [BALanced(integer -999999)] [BASEperiod(real -999999)] ///
        [T_min(integer -999999)] [T_max(integer -999999)]
    
    // =========================================================================
    // Header Section
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    display as text _col(18) "Maximum Equivalence Test Results"
    display as text "{hline 78}"
    
    // =========================================================================
    // Test Configuration
    // =========================================================================
    
    display ""
    display as text "Test Configuration:"
    
    // Right-align all values in 25-char field starting at column 30
    if "`method'" == "IU" {
        display as text "  Test type:" _col(30) as result %25s "Maximum (Intersection-Union)"
    }
    else if "`method'" == "Boot" {
        display as text "  Test type:" _col(30) as result %25s "Maximum (Bootstrap)"
    }
    else if "`method'" == "Wild" {
        display as text "  Test type:" _col(30) as result %25s "Maximum (Wild Cluster Bootstrap)"
    }
    else {
        display as text "  Test type:" _col(30) as result %25s "Maximum (`method')"
    }
    
    display as text "  Significance level:" _col(30) as result %25.4f `alpha'
    
    if `threshold_specified' == 1 {
        display as text "  Equiv. threshold:" _col(30) as result %25.6f `threshold'
    }
    else {
        display as text "  Equiv. threshold:" _col(30) as result %25s "(searching for minimum)"
    }
    
    // Bootstrap replications (if applicable)
    if "`method'" == "Boot" | "`method'" == "Wild" {
        display as text "  Bootstrap reps:" _col(30) as result %25.0fc `reps'
    }
    
    // =========================================================================
    // Hypothesis Statement
    // =========================================================================
    
    display ""
    display as text "Hypothesis Test:"
    display as text "  H0: max|placebo effect| >= threshold  (non-equivalence)"
    display as text "  H1: max|placebo effect| <  threshold  (equivalence)"
    
    // =========================================================================
    // Placebo Coefficients Table
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    display as text "Estimated Placebo Coefficients"
    display as text "{hline 78}"
    
    // Table header with vertical separator
    display as text _col(3) "Period" _col(16) "{c |}" ///
        _col(20) "Coefficient" ///
        _col(38) "Abs. Value"
    display as text "{hline 15}{c +}{hline 62}"
    
    // Display coefficients from string parameter or stored matrix
    if "`placebo_coefs'" != "" {
        local i = 1
        foreach coef of local placebo_coefs {
            local period_label : word `i' of `placebo_names'
            if "`period_label'" == "" {
                local period_label = "placebo_`i'"
            }
            local abs_coef = abs(`coef')
            display as text _col(3) %12s abbrev("`period_label'", 12) _col(16) "{c |}" ///
                as result _col(20) %12.6f `coef' ///
                _col(38) %12.6f `abs_coef'
            local i = `i' + 1
        }
    }
    else {
        capture matrix placebo = r(placebo_coefs)
        if _rc == 0 {
            forvalues i = 1/`no_placebos' {
                local coef = placebo[`i', 1]
                local period_label : word `i' of `placebo_names'
                if "`period_label'" == "" {
                    local period_label = "placebo_`i'"
                }
                local abs_coef = abs(`coef')
                display as text _col(3) %12s abbrev("`period_label'", 12) _col(16) "{c |}" ///
                    as result _col(20) %12.6f `coef' ///
                    _col(38) %12.6f `abs_coef'
            }
        }
        else {
            display as text _col(3) "(coefficients not available)"
        }
    }
    
    display as text "{hline 78}"
    
    // =========================================================================
    // Test Results Table
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    
    if `threshold_specified' == 1 {
        // ---------------------------------------------------------------------
        // Mode 1: User-specified threshold
        // ---------------------------------------------------------------------
        
        display as text "Test Results (Equivalence Threshold = " as result %9.6f `threshold' as text ")"
        display as text "{hline 78}"
        
        // Table header with vertical separator (right-aligned)
        display as text _col(3) "Statistic" _col(16) "{c |}" ///
            _col(22) %12s "Max |Coef|" ///
            _col(38) %12s "Crit. Value" ///
            _col(56) %9s "Reject H0"
        display as text "{hline 15}{c +}{hline 62}"
        
        // Rejection indicator
        if `reject' == 1 {
            local reject_str "Yes"
        }
        else {
            local reject_str "No"
        }
        
        // Data row
        display as text _col(3) "Max Placebo" _col(16) "{c |}" ///
            as result _col(22) %12.6f `max_abs_coef' ///
            _col(38) %12.6f `critical_value' ///
            as text _col(56) %9s "`reject_str'"
        
        display as text "{hline 78}"
        
        // Summary decision
        display ""
        if `reject' == 1 {
            display as result "Decision: REJECT H0" as text " at alpha = " as result %5.3f `alpha'
            display as text "Conclusion: Evidence supports equivalence of pre-trends"
        }
        else {
            display as result "Decision: FAIL TO REJECT H0" as text " at alpha = " as result %5.3f `alpha'
            display as text "Conclusion: Insufficient evidence for equivalence"
        }
    }
    else {
        // ---------------------------------------------------------------------
        // Mode 2: Minimum threshold search
        // ---------------------------------------------------------------------
        
        display as text "Minimum Equivalence Threshold Search"
        display as text "{hline 78}"
        
        // Table header with vertical separator (right-aligned)
        display as text _col(3) "Statistic" _col(16) "{c |}" ///
            _col(22) %12s "Max |Coef|" ///
            _col(42) %23s "Min. Threshold"
        display as text "{hline 15}{c +}{hline 62}"
        
        // Data row
        display as text _col(3) "Max Placebo" _col(16) "{c |}" ///
            as result _col(22) %12.6f `max_abs_coef' ///
            _col(42) %23.6f `min_threshold'
        
        display as text "{hline 78}"
        
        // Interpretation
        display ""
        display as text "Minimum equivalence threshold = " as result %10.6f `min_threshold'
        display ""
        display as text "Interpretation:"
        display as text "  The minimum threshold is the smallest value at which H0 can be rejected"
        display as text "  at the " as result %4.1f `alpha'*100 as text "% significance level."
    }
    
    // =========================================================================
    // Panel Information Footer
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    display as text "Panel Information"
    display as text "{hline 78}"
    
    display as text "  Observations:            " as result %12.0fc `n_total'
    display as text "  Number of groups:        " as result %12.0fc `n'
    
    // Pre-treatment periods (handle balanced/unbalanced)
    if `balanced' != -999999 {
        if `balanced' == 1 {
            display as text "  Pre-treatment periods:   " as result %12.0fc `n_t'
        }
        else {
            if `t_min' != -999999 & `t_max' != -999999 {
                display as text "  Pre-treatment periods:        " as result `t_min' as text " - " as result `t_max' as text " (range)"
            }
            else {
                display as text "  Pre-treatment periods:   " as result %12.0fc `n_t'
            }
        }
    }
    else {
        display as text "  Pre-treatment periods:   " as result %12.0fc `n_t'
    }
    
    display as text "  Placebo coefficients:    " as result %12.0fc `no_placebos'
    
    if `baseperiod' != -999999 {
        display as text "  Base period:             " as result %12.0g `baseperiod'
    }
    
    // Panel type
    if `balanced' != -999999 {
        if `balanced' == 1 {
            display as text "  Panel type:              " as result %12s "Balanced"
        }
        else {
            display as text "  Panel type:              " as result %12s "Unbalanced"
        }
    }
    
    display as text "{hline 78}"
    display ""
    
end
