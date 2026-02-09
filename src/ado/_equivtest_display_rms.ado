*! _equivtest_display_rms.ado - Display results for RMS equivalence test
*! Version 2.0 - Professional Stata-style output with SMCL formatting
*!
*! Formats and displays the root mean square (RMS) placebo effect test.
*! Tests H0: beta_RMS >= zeta vs H1: beta_RMS < zeta
*!
*! The RMS statistic: beta_RMS = ||beta|| / sqrt(T)
*!
*! Two display modes:
*!   1. Threshold specified: displays RMS, simulated critical value, decision
*!   2. Minimum threshold: displays RMS, minimum threshold

program define _equivtest_display_rms
    version 16.0
    
    syntax, THRESHOLD_specified(integer) Alpha(real) ///
        [THRESHold(real -999999)] [REJect(integer -999999)] ///
        [MIN_threshold(real -999999)] ///
        [RMS_placebo(real -999999)] [RMS_critical(real -999999)] ///
        [SEED(integer -999999)]
    
    // =========================================================================
    // Hypothesis Statement
    // =========================================================================
    
    display ""
    display as text "Hypothesis Test:"
    display as text "  H0: RMS placebo effect >= zeta  (non-equivalence)"
    display as text "  H1: RMS placebo effect <  zeta  (equivalence)"
    display ""
    display as text "Note: RMS test uses self-normalized subsampling for inference."
    display as text "      Results may vary between runs. Set seed() for reproducibility."
    if `seed' != -999999 {
        display as text "      Random seed used: " as result `seed'
    }
    
    // =========================================================================
    // Results Table
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    
    if `threshold_specified' == 1 {
        // ---------------------------------------------------------------------
        // Mode 1: User-specified threshold
        // ---------------------------------------------------------------------
        
        display as text "Test Results (Equivalence Threshold zeta = " as result %9.6f `threshold' as text ")"
        display as text "{hline 78}"
        
        // Table header (right-aligned)
        display as text _col(3) "Statistic" _col(16) "{c |}" ///
            _col(22) %12s "RMS Value" ///
            _col(38) %16s "Simul. Crit. Val." ///
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
        display as text _col(3) "RMS Placebo" _col(16) "{c |}" ///
            as result _col(22) %12.6f `rms_placebo' ///
            _col(38) %16.6f `rms_critical' ///
            as text _col(60) %9s "`reject_str'"
        
        display as text "{hline 78}"
        
        // Summary decision
        display ""
        if `reject' == 1 {
            display as result "Decision: REJECT H0" as text " at alpha = " as result %5.3f `alpha'
            display as text "Conclusion: Evidence supports equivalence of pre-trends"
            display as text "            (RMS placebo effect is negligibly small)"
        }
        else {
            display as result "Decision: FAIL TO REJECT H0" as text " at alpha = " as result %5.3f `alpha'
            display as text "Conclusion: Insufficient evidence for equivalence"
            display as text "            (cannot conclude RMS placebo effect is small)"
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
            _col(22) %12s "RMS Value" ///
            _col(42) %22s "Min. Threshold (zeta*)"
        display as text "{hline 15}{c +}{hline 62}"
        
        // Data row
        display as text _col(3) "RMS Placebo" _col(16) "{c |}" ///
            as result _col(22) %12.6f `rms_placebo' ///
            _col(42) %22.6f `min_threshold'
        
        display as text "{hline 78}"
        
        // Interpretation
        display ""
        display as text "Minimum equivalence threshold zeta* = " as result %10.6f `min_threshold'
        display ""
        display as text "Interpretation:"
        display as text "  zeta* is the smallest threshold at which H0 can be rejected"
        display as text "  at the " as result %4.1f `alpha'*100 as text "% significance level."
    }
    
end
