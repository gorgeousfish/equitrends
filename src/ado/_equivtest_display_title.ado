*! _equivtest_display_title.ado - Display formatted header for equivalence tests
*! Version 2.0 - Professional Stata-style output with SMCL formatting
*!
*! Displays a formatted header using standard Stata SMCL directives.
*! Accepts optional parameters to display test configuration in the header.

program define _equivtest_display_title
    version 16.0
    
    syntax [, TYPE(string) Method(string) Alpha(real -999999) ///
        THRESHold(real -999999) THRESHOLD_specified(integer 0) ///
        VCE(string) Nboot(integer -999999)]
    
    // =========================================================================
    // Header Section
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    display as text _col(14) "Equivalence Tests for Pre-trends in DiD Estimation"
    display as text "{hline 78}"
    
    // =========================================================================
    // Test Configuration (if parameters provided)
    // =========================================================================
    
    if "`type'" != "" {
        display ""
        display as text "Test Configuration:"
        
        // Use _col() to position at column 30, then right-justify values in 25-char field
        if "`type'" == "max" {
            if "`method'" == "iu" | "`method'" == "" {
                display as text "  Test type:" _col(30) as result %25s "Maximum (Intersection-Union)"
            }
            else if "`method'" == "boot" {
                display as text "  Test type:" _col(30) as result %25s "Maximum (Bootstrap)"
            }
            else if "`method'" == "wild" {
                display as text "  Test type:" _col(30) as result %25s "Maximum (Wild Cluster Bootstrap)"
            }
        }
        else if "`type'" == "mean" {
            display as text "  Test type:" _col(30) as result %25s "Mean Placebo Effect"
        }
        else if "`type'" == "rms" {
            display as text "  Test type:" _col(30) as result %25s "Root Mean Square (RMS)"
        }
        
        // Significance level
        if `alpha' != -999999 {
            display as text "  Significance level:" _col(30) as result %25.4f `alpha'
        }
        
        // Equivalence threshold with mathematical symbol
        if `threshold_specified' == 1 & `threshold' != -999999 {
            if "`type'" == "max" {
                display as text "  Equiv. threshold (delta):" _col(30) as result %25.6f `threshold'
            }
            else if "`type'" == "mean" {
                display as text "  Equiv. threshold (tau):" _col(30) as result %25.6f `threshold'
            }
            else if "`type'" == "rms" {
                display as text "  Equiv. threshold (zeta):" _col(30) as result %25.6f `threshold'
            }
        }
        else if `threshold_specified' == 0 {
            if "`type'" == "max" {
                display as text "  Equiv. threshold (delta):" _col(30) as result %25s "(searching for minimum)"
            }
            else if "`type'" == "mean" {
                display as text "  Equiv. threshold (tau):" _col(30) as result %25s "(searching for minimum)"
            }
            else if "`type'" == "rms" {
                display as text "  Equiv. threshold (zeta):" _col(30) as result %25s "(searching for minimum)"
            }
            else {
                display as text "  Equiv. threshold:" _col(30) as result %25s "(searching for minimum)"
            }
        }
        
        // VCE type
        if "`vce'" != "" {
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
        }
        
        // Bootstrap replications
        if `nboot' != -999999 & inlist("`method'", "boot", "wild") {
            display as text "  Bootstrap reps:" _col(30) as result %25.0fc `nboot'
        }
    }
    
end
