*! _equivtest_display_boot.ado - Display results for bootstrap equivalence test
*! Version 2.0 - Professional Stata-style output with SMCL formatting
*!
*! Formats and displays bootstrap-based equivalence test output.
*! Supports both spherical errors bootstrap and cluster wild bootstrap.
*!
*! Two display modes:
*!   1. Threshold specified: displays max|coef|, bootstrap critical value, decision
*!   2. Minimum threshold: displays max|coef|, minimum threshold

program define _equivtest_display_boot
    version 16.0
    
    syntax, Method(string) THRESHOLD_specified(integer) Alpha(real) ///
        [THRESHold(real -999999)] [REJect(integer -999999)] ///
        [MIN_threshold(real -999999)] [Nboot(integer 1000)] ///
        [MAX_abs_coef(real -999999)] [BOOT_critical(real -999999)]
    
    // =========================================================================
    // Hypothesis Statement
    // =========================================================================
    
    display ""
    display as text "Hypothesis Test:"
    display as text "  H0: max|placebo effect| >= delta  (non-equivalence)"
    display as text "  H1: max|placebo effect| <  delta  (equivalence)"
    display ""
    
    // Bootstrap method description
    if "`method'" == "boot" {
        display as text "Method: Parametric bootstrap for spherical errors"
    }
    else {
        display as text "Method: Wild cluster bootstrap"
    }
    display as text "Bootstrap replications: " as result `nboot'
    
    // =========================================================================
    // Results Table
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    
    if `threshold_specified' == 1 {
        // ---------------------------------------------------------------------
        // Mode 1: User-specified threshold
        // ---------------------------------------------------------------------
        
        display as text "Test Results (Equivalence Threshold delta = " as result %9.6f `threshold' as text ")"
        display as text "{hline 78}"
        
        // Table header (right-aligned)
        display as text _col(3) "Statistic" _col(16) "{c |}" ///
            _col(22) %12s "Max |Coef|" ///
            _col(38) %16s "Boot. Crit. Val." ///
            _col(60) %9s "Reject H0"
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
            _col(38) %16.6f `boot_critical' ///
            as text _col(60) %9s "`reject_str'"
        
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
        
        // Table header (right-aligned)
        display as text _col(3) "Statistic" _col(16) "{c |}" ///
            _col(22) %12s "Max |Coef|" ///
            _col(42) %23s "Min. Threshold (delta*)"
        display as text "{hline 15}{c +}{hline 62}"
        
        // Data row
        display as text _col(3) "Max Placebo" _col(16) "{c |}" ///
            as result _col(22) %12.6f `max_abs_coef' ///
            _col(42) %23.6f `min_threshold'
        
        display as text "{hline 78}"
        
        // Interpretation
        display ""
        display as text "Minimum equivalence threshold delta* = " as result %10.6f `min_threshold'
        display ""
        display as text "Interpretation:"
        display as text "  delta* is the smallest threshold at which H0 can be rejected"
        display as text "  at the " as result %4.1f `alpha'*100 as text "% significance level."
    }
    
end
