*! _rmsequivtest_display.ado - Display output for RMS equivalence test
*! Version 2.0 - Professional Stata-style output with SMCL formatting
*!
*! Formats and displays results from the RMS equivalence test for pre-trends.
*! The RMS test uses self-normalized subsampling for inference.
*!
*! Two display modes:
*!   1. Threshold specified: displays RMS, simulated critical value, decision
*!   2. Minimum threshold: displays RMS, minimum threshold

program define _rmsequivtest_display
    version 16.0
    
    syntax, RMS_placebo(real) THRESHOLD_specified(integer) ///
        Alpha(real) N(integer) N_g(integer) N_t(integer) NO_placebos(integer) ///
        BASEperiod(real) BALanced(integer) ///
        RMS_ci_lower(real) RMS_ci_upper(real) ///
        MS_ci_lower(real) MS_ci_upper(real) LEVel(real) ///
        [THRESHold(real -999999)] [RMS_critical(real -999999)] ///
        [REJect(integer -999999)] [MIN_threshold(real -999999)] ///
        [SEED(integer -999999)] [T_min(integer -999999)] [T_max(integer -999999)]
    
    // =========================================================================
    // Header
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    display as text _col(14) "Equivalence Tests for Pre-trends in DiD Estimation"
    display as text "{hline 78}"
    
    // =========================================================================
    // Test Configuration
    // =========================================================================
    
    display ""
    display as text "Test Configuration:"
    display as text "  Test type:" _col(30) as result %25s "Root Mean Square (RMS)"
    display as text "  Significance level:" _col(30) as result %25.4f `alpha'
    display as text "  Confidence level:" _col(30) as result %24.0f `level' as result "%"
    
    if `threshold_specified' == 1 & `threshold' != -999999 {
        display as text "  Equiv. threshold (zeta):" _col(30) as result %25.6f `threshold'
    }
    else {
        display as text "  Equiv. threshold (zeta):" _col(30) as result %25s "(searching for minimum)"
    }
    
    // =========================================================================
    // Important Note about Subsampling
    // =========================================================================
    
    display ""
    display as result "Note: " as text "RMS test uses self-normalized subsampling for inference."
    display as text "      Results may vary between runs. Set seed() for reproducibility."
    if `seed' != -999999 {
        display as text "      Random seed used: " as result `seed'
    }
    
    // =========================================================================
    // Hypothesis Statement
    // =========================================================================
    
    display ""
    display as text "Hypothesis Test:"
    display as text "  H0: RMS placebo effect >= zeta  (non-equivalence)"
    display as text "  H1: RMS placebo effect <  zeta  (equivalence)"
    
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
            _col(24) %10s "RMS Value" ///
            _col(40) %16s "Simul. Crit. Val." ///
            _col(62) %9s "Reject H0"
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
            as result _col(24) %10.6f `rms_placebo' ///
            _col(40) %16.6f `rms_critical' ///
            as text _col(62) %9s "`reject_str'"
        
        display as text "{hline 78}"
        
        // Confidence Interval display
        display ""
        display as text "Confidence Intervals (based on self-normalized inference):"
        display as text "  RMS Placebo Effect:" _col(30) as result %10.6f `rms_placebo' ///
            as text "  [" as result %6.0f `level' as text "% CI: " ///
            as result %9.6f `rms_ci_lower' as text ", " as result %9.6f `rms_ci_upper' as text "]"
        display as text "  MS Placebo Effect:" _col(30) as result %10.6f `rms_placebo'^2 ///
            as text "  [" as result %6.0f `level' as text "% CI: " ///
            as result %9.6f `ms_ci_lower' as text ", " as result %9.6f `ms_ci_upper' as text "]"
        
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
        
        // Table header
        display as text _col(3) "Statistic" _col(16) "{c |}" ///
            _col(20) "RMS Value" ///
            _col(40) "Min. Threshold (zeta*)"
        display as text "{hline 15}{c +}{hline 62}"
        
        // Data row
        display as text _col(3) "RMS Placebo" _col(16) "{c |}" ///
            as result _col(20) %12.6f `rms_placebo' ///
            _col(40) %14.6f `min_threshold'
        
        display as text "{hline 78}"
        
        // Confidence Interval display
        display ""
        display as text "Confidence Intervals (based on self-normalized inference):"
        display as text "  RMS Placebo Effect:" _col(30) as result %10.6f `rms_placebo' ///
            as text "  [" as result %6.0f `level' as text "% CI: " ///
            as result %9.6f `rms_ci_lower' as text ", " as result %9.6f `rms_ci_upper' as text "]"
        display as text "  MS Placebo Effect:" _col(30) as result %10.6f `rms_placebo'^2 ///
            as text "  [" as result %6.0f `level' as text "% CI: " ///
            as result %9.6f `ms_ci_lower' as text ", " as result %9.6f `ms_ci_upper' as text "]"
        
        // Interpretation
        display ""
        display as text "Minimum equivalence threshold zeta* = " as result %10.6f `min_threshold'
        display ""
        display as text "Interpretation:"
        display as text "  zeta* is the smallest threshold at which H0 can be rejected"
        display as text "  at the " as result %4.1f `alpha'*100 as text "% significance level."
    }
    
    // =========================================================================
    // Panel Information
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    display as text "Panel Information"
    display as text "{hline 78}"
    
    display as text "  Observations:            " as result %12.0fc `n'
    display as text "  Number of groups:        " as result %12.0fc `n_g'
    
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
    
    if `balanced' == 1 {
        display as text "  Panel type:              " as result %12s "Balanced"
    }
    else {
        display as text "  Panel type:              " as result %12s "Unbalanced"
    }
    
    display as text "{hline 78}"
    display ""
    
end
