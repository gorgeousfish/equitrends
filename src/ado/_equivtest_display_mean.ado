*! _equivtest_display_mean.ado - Display results for mean placebo effect test
*! Version 2.0 - Professional Stata-style output with SMCL formatting
*!
*! Formats and displays the mean placebo effect equivalence test output.
*! Tests H0: |beta_bar| >= tau vs H1: |beta_bar| < tau
*!
*! Two display modes:
*!   1. Threshold specified: displays estimate, SE, p-value, critical value
*!   2. Minimum threshold: displays estimate, SE, minimum threshold

program define _equivtest_display_mean
    version 16.0
    
    syntax, THRESHOLD_specified(integer) Alpha(real) ///
        [THRESHold(real -999999)] [REJect(integer -999999)] ///
        [MIN_threshold(real -999999)] ///
        [ABS_mean_placebo(real -999999)] [SE_mean_placebo(real -999999)] ///
        [P_value(real -999999)] [MEAN_critical(real -999999)]
    
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
        
        // Table header (right-aligned to match data columns)
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
            _col(52) %11.6f `mean_critical' ///
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
        
        // Table header (right-aligned to match data columns)
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
    
end
