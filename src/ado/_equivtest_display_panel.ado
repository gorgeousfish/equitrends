*! _equivtest_display_panel.ado - Display panel information for equivalence tests
*! Version 2.0 - Professional Stata-style output with SMCL formatting
*!
*! Displays panel structure information including:
*!   - Number of observations
*!   - Number of individuals/groups
*!   - Number of time periods
*!   - Number of placebo coefficients
*!   - Base period
*!   - Panel balance status

program define _equivtest_display_panel
    version 16.0
    
    syntax, NO_placebos(integer) BASEperiod(real) ///
        BALanced(integer) N_g(integer) Nobs(integer) ///
        [TYPE(string)] [Method(string)] ///
        [Tmin(integer -999999)] [Tmax(integer -999999)] ///
        [Tcount(integer -999999)]
    
    // =========================================================================
    // Panel Information Section
    // =========================================================================
    
    display ""
    display as text "{hline 78}"
    display as text "Panel Information"
    display as text "{hline 78}"
    
    // Number of observations
    display as text "  Observations:            " as result %12.0fc `nobs'
    
    // Number of individuals/groups
    display as text "  Number of groups:        " as result %12.0fc `n_g'
    
    // Number of pre-treatment periods
    if `balanced' == 1 {
        if `tcount' != -999999 {
            display as text "  Pre-treatment periods:   " as result %12.0fc `tcount'
        }
    }
    else {
        if `tmin' != -999999 & `tmax' != -999999 {
            display as text "  Pre-treatment periods:        " as result `tmin' as text " - " as result `tmax' as text " (range)"
        }
        else if `tcount' != -999999 {
            display as text "  Pre-treatment periods:   " as result %12.0fc `tcount'
        }
    }
    
    // Number of placebo coefficients estimated
    display as text "  Placebo coefficients:    " as result %12.0fc `no_placebos'
    
    // Base period
    display as text "  Base period:             " as result %12.0g `baseperiod'
    
    // Panel type (balanced/unbalanced)
    if `balanced' == 1 {
        display as text "  Panel type:              " as result %12s "Balanced"
    }
    else {
        display as text "  Panel type:              " as result %12s "Unbalanced"
    }
    
    display as text "{hline 78}"
    display ""
    
end
