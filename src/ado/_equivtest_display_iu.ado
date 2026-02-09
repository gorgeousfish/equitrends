*! _equivtest_display_iu.ado - Display results for Intersection-Union (IU) test
*! Version 2.0 - Professional Stata-style output with SMCL formatting
*!
*! Formats and displays the IU equivalence test output using standard Stata
*! table formatting with proper column alignment and SMCL directives.
*!
*! Two display modes:
*!   1. Threshold specified: displays estimates, SEs, critical values, decision
*!   2. Minimum threshold: displays estimates, SEs, minimum thresholds

program define _equivtest_display_iu
    version 16.0
    
    syntax, THRESHOLD_specified(integer) Alpha(real) ///
        [THRESHold(real -999999)] [REJect(integer -999999)] ///
        [MIN_threshold(real -999999)] ///
        [PLACEBO_coefs(string)] [SE_placebo(string)] ///
        [CRITICAL_values(string)] [MIN_thresholds(string)] ///
        [PLACEBO_names(string)]
    
    // =========================================================================
    // Hypothesis Statement
    // =========================================================================
    
    display ""
    display as text "Hypothesis Test:"
    display as text "  H0: max|placebo effect| >= delta  (non-equivalence)"
    display as text "  H1: max|placebo effect| <  delta  (equivalence)"
    
    // =========================================================================
    // Results Table
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    
    // Get number of placebo coefficients
    local nrows : word count `placebo_coefs'
    if `nrows' == 0 {
        local nrows = 1
    }
    
    if `threshold_specified' == 1 {
        // ---------------------------------------------------------------------
        // Mode 1: User-specified threshold - show critical values
        // ---------------------------------------------------------------------
        
        display as text "Test Results (Equivalence Threshold delta = " as result %9.6f `threshold' as text ")"
        display as text "{hline 78}"
        
        // Table header with vertical separator (right-aligned)
        display as text _col(3) "Period" _col(16) "{c |}" ///
            _col(22) %13s "Abs. Estimate" ///
            _col(38) %12s "Std. Error" ///
            _col(54) %12s "Crit. Value" ///
            _col(70) %6s "Reject"
        display as text "{hline 15}{c +}{hline 62}"
        
        // Data rows
        forvalues i = 1/`nrows' {
            local coef : word `i' of `placebo_coefs'
            local se : word `i' of `se_placebo'
            local cv : word `i' of `critical_values'
            local pname : word `i' of `placebo_names'
            
            // Default period name if not provided
            if "`pname'" == "" {
                local pname "placebo_`i'"
            }
            
            // Absolute value of coefficient
            local abs_coef = abs(`coef')
            
            // Determine rejection for this coefficient
            if `abs_coef' < `cv' {
                local row_reject "Yes"
            }
            else {
                local row_reject "No"
            }
            
            display as text _col(3) %12s abbrev("`pname'", 12) _col(16) "{c |}" ///
                as result _col(22) %13.6f `abs_coef' ///
                _col(38) %12.6f `se' ///
                _col(54) %12.6f `cv' ///
                as text _col(70) %6s "`row_reject'"
        }
        
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
        display as text _col(3) "Period" _col(16) "{c |}" ///
            _col(22) %13s "Abs. Estimate" ///
            _col(38) %12s "Std. Error" ///
            _col(54) %14s "Min. Threshold"
        display as text "{hline 15}{c +}{hline 62}"
        
        // Data rows
        forvalues i = 1/`nrows' {
            local coef : word `i' of `placebo_coefs'
            local se : word `i' of `se_placebo'
            local mt : word `i' of `min_thresholds'
            local pname : word `i' of `placebo_names'
            
            // Default period name if not provided
            if "`pname'" == "" {
                local pname "placebo_`i'"
            }
            
            // Absolute value of coefficient
            local abs_coef = abs(`coef')
            
            display as text _col(3) %12s abbrev("`pname'", 12) _col(16) "{c |}" ///
                as result _col(22) %13.6f `abs_coef' ///
                _col(38) %12.6f `se' ///
                _col(54) %14.6f `mt'
        }
        
        display as text "{hline 78}"
        
        // Overall minimum threshold
        display ""
        display as text "Minimum equivalence threshold delta* = " as result %10.6f `min_threshold'
        display ""
        display as text "Interpretation:"
        display as text "  delta* is the smallest threshold at which H0 can be rejected"
        display as text "  at the " as result %4.1f `alpha'*100 as text "% significance level."
    }
    
end
