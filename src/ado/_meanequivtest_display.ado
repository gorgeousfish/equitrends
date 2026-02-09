*! _meanequivtest_display.ado - Display output for the mean equivalence test
*! Version 2.0 - Professional Stata-style output with SMCL formatting
*!
*! This module formats and displays results for hypothesis (3.2) from
*! Dette & Schumann (2024):
*!   H0: |beta_bar| >= tau  vs  H1: |beta_bar| < tau
*! where beta_bar = (1/T) * sum_{t=1}^{T} beta_t is the mean placebo effect.
*!
*! Two output modes are supported:
*!   (1) Threshold specified: displays test decision and p-value
*!   (2) Threshold not specified: reports minimum threshold tau*

program define _meanequivtest_display
    version 16.0
    
    syntax, ABS_mean_placebo(real) SE_mean_placebo(real) ///
        THRESHOLD_specified(integer) ///
        Alpha(real) N(integer) N_g(integer) N_t(integer) NO_placebos(integer) ///
        BASEperiod(real) BALanced(integer) VCE(string) ///
        [THRESHold(real -999999)] [CRITICAL_value(real -999999)] ///
        [P_value(real -999999)] [REJect(integer -999999)] ///
        [MIN_threshold(real -999999)] ///
        [T_min(integer -999999)] [T_max(integer -999999)]
    
    // =========================================================================
    // Header Section
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    display as text _col(20) "Mean Equivalence Test Results"
    display as text "{hline 78}"
    
    // =========================================================================
    // Test Configuration
    // =========================================================================
    
    display ""
    display as text "Test Configuration:"
    display as text "  Test type:" _col(30) as result %25s "Mean Placebo Effect"
    display as text "  Significance level:" _col(30) as result %25.4f `alpha'
    
    if `threshold_specified' == 1 {
        display as text "  Equiv. threshold (tau):" _col(30) as result %25.6f `threshold'
    }
    else {
        display as text "  Equiv. threshold (tau):" _col(30) as result %25s "(searching for minimum)"
    }
    
    // VCE type (right-aligned in 25-char field)
    if "`vce'" == "ols" {
        display as text "  Variance estimator:" _col(30) as result %25s "OLS (homoskedastic)"
    }
    else if "`vce'" == "hc1" {
        display as text "  Variance estimator:" _col(30) as result %25s "HC1 (heteroskedasticity-robust)"
    }
    else if "`vce'" == "hc3" {
        display as text "  Variance estimator:" _col(30) as result %25s "HC3 (heteroskedasticity-robust)"
    }
    else if "`vce'" == "cr0" {
        display as text "  Variance estimator:" _col(30) as result %25s "CR0 (cluster-robust)"
    }
    else {
        display as text "  Variance estimator:" _col(30) as result %25s "`vce'"
    }
    
    // =========================================================================
    // Hypothesis Statement
    // =========================================================================
    
    display ""
    display as text "Hypothesis Test:"
    display as text "  H0: |mean placebo effect| >= tau  (non-equivalence)"
    display as text "  H1: |mean placebo effect| <  tau  (equivalence)"
    
    // =========================================================================
    // Results Table
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    
    if `threshold_specified' == 1 {
        // ---------------------------------------------------------------------
        // Mode 1: User-specified threshold
        // ---------------------------------------------------------------------
        
        display as text "Test Results (Equivalence Threshold tau = " as result %9.6f `threshold' as text ")"
        display as text "{hline 78}"
        
        // Table header with vertical separator (right-aligned)
        display as text _col(3) "Statistic" _col(16) "{c |}" ///
            _col(22) %10s "Abs. Mean" ///
            _col(37) %10s "Std. Error" ///
            _col(52) %11s "Crit. Value" ///
            _col(67) %10s "p-value"
        display as text "{hline 15}{c +}{hline 62}"
        
        // Data row
        display as text _col(3) "Mean Placebo" _col(16) "{c |}" ///
            as result _col(22) %10.6f `abs_mean_placebo' ///
            _col(37) %10.6f `se_mean_placebo' ///
            _col(52) %11.6f `critical_value' ///
            _col(67) %10.6f `p_value'
        
        display as text "{hline 78}"
        
        // Summary decision
        display ""
        if `reject' == 1 {
            display as result "Decision: REJECT H0" as text " at alpha = " as result %5.3f `alpha'
            display as text "Conclusion: Evidence supports equivalence of pre-trends"
            display as text "            (mean placebo effect is negligibly small)"
        }
        else {
            display as result "Decision: FAIL TO REJECT H0" as text " at alpha = " as result %5.3f `alpha'
            display as text "Conclusion: Insufficient evidence for equivalence"
            display as text "            (cannot conclude mean placebo effect is small)"
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
            _col(22) %10s "Abs. Mean" ///
            _col(37) %10s "Std. Error" ///
            _col(54) %14s "Min. Threshold"
        display as text "{hline 15}{c +}{hline 62}"
        
        // Data row
        display as text _col(3) "Mean Placebo" _col(16) "{c |}" ///
            as result _col(22) %10.6f `abs_mean_placebo' ///
            _col(37) %10.6f `se_mean_placebo' ///
            _col(54) %14.6f `min_threshold'
        
        display as text "{hline 78}"
        
        // Interpretation
        display ""
        display as text "Minimum equivalence threshold tau* = " as result %10.6f `min_threshold'
        display ""
        display as text "Interpretation:"
        display as text "  tau* is the smallest threshold at which H0 can be rejected"
        display as text "  at the " as result %4.1f `alpha'*100 as text "% significance level."
    }
    
    // =========================================================================
    // Panel Information Footer
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    display as text "Panel Information"
    display as text "{hline 78}"
    
    display as text "  Observations:            " as result %12.0fc `n'
    display as text "  Number of groups:        " as result %12.0fc `n_g'
    
    // Pre-treatment periods (handle balanced/unbalanced)
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
    
    display as text "  Placebo coefficients:    " as result %12.0fc `no_placebos'
    display as text "  Base period:             " as result %12.0g `baseperiod'
    
    // Panel type
    if `balanced' == 1 {
        display as text "  Panel type:              " as result %12s "Balanced"
    }
    else {
        display as text "  Panel type:              " as result %12s "Unbalanced"
    }
    
    display as text "{hline 78}"
    display ""
    
end
