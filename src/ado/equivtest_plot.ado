*! equivtest_plot.ado - Visualization of placebo coefficient estimates
*!
*! Generates coefficient plots for pre-treatment placebo effects from
*! equivalence test results. Point estimates are displayed with optional
*! confidence intervals and equivalence threshold bounds (+/- delta).

program define equivtest_plot
    version 16.0
    
    syntax [, THRESHold(real -999999) CI Level(real 95) ///
        CONnect NOLine NOBase NOTHRESHold ///
        MSYMbol(string) MSIze(string) MCOlor(string) ///
        BASEMSymbol(string) BASEMColor(string) ///
        CILColor(string) CILWidth(string) ///
        THRESHLColor(string) THRESHLWidth(string) THRESHLPattern(string) ///
        SCHeme(string) TItle(string) SUBtitle(string) ///
        XTItle(string) YTItle(string) XLAbel(string) YLAbel(string) ///
        LEGend(string asis) NOTE(string asis) ///
        SAVing(string) REPLACE NAME(string)]
    
    // -------------------------------------------------------------------------
    // Validation
    // -------------------------------------------------------------------------
    
    if "`e(cmd)'" != "equivtest" {
        di as error "equivtest_plot requires equivtest results"
        exit 301
    }
    
    capture confirm matrix e(b_placebo)
    if _rc != 0 {
        di as error "placebo coefficients not found in e()"
        exit 301
    }
    
    // -------------------------------------------------------------------------
    // Extraction of stored estimation results
    // -------------------------------------------------------------------------
    
    tempname b_placebo
    matrix `b_placebo' = e(b_placebo)
    local no_placebos = e(no_placebos)
    
    local base_period = e(base_period)
    local preperiods "`e(preperiods)'"
    
    local test_type "`e(type)'"
    local test_method "`e(method)'"
    
    // Period-specific standard errors are available only for IU tests
    local se_available = 0
    if "`test_type'" == "max" & "`test_method'" == "iu" {
        capture confirm matrix e(se_placebo)
        if _rc == 0 {
            tempname se_placebo
            matrix `se_placebo' = e(se_placebo)
            local se_available = 1
        }
    }
    
    // -------------------------------------------------------------------------
    // Confidence interval option processing
    // -------------------------------------------------------------------------
    
    // Confidence intervals require period-specific standard errors (IU test only)
    if "`ci'" != "" & `se_available' == 0 {
        di as text "Warning: ci option ignored. Confidence intervals are only available for type(max) method(iu)."
        di as text "         Bootstrap, Mean, and RMS tests do not provide period-specific standard errors."
        local ci ""
    }
    
    if `level' <= 0 | `level' >= 100 {
        di as error "level() must be between 0 and 100"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Equivalence threshold determination
    // -------------------------------------------------------------------------
    
    local plot_threshold = .
    local has_threshold = 0
    local thresh_source ""
    
    if "`nothreshold'" == "" {
        // Priority: user-specified > e(threshold) > e(min_threshold)
        if `threshold' != -999999 {
            if `threshold' <= 0 {
                di as error "equiv_threshold must be strictly positive (> 0)"
                exit 198
            }
            local plot_threshold = `threshold'
            local has_threshold = 1
            local thresh_source "user"
        }
        else {
            local thresh_spec = e(threshold_specified)
            if `thresh_spec' == 1 {
                capture scalar _temp_thresh = e(threshold)
                if _rc == 0 {
                    local plot_threshold = e(threshold)
                    local has_threshold = 1
                    local thresh_source "specified"
                }
            }
            else {
                capture scalar _temp_thresh = e(min_threshold)
                if _rc == 0 {
                    local plot_threshold = e(min_threshold)
                    local has_threshold = 1
                    local thresh_source "minimum"
                }
            }
        }
    }
    
    // -------------------------------------------------------------------------
    // Extraction of test-specific summary statistics for annotations
    // -------------------------------------------------------------------------
    
    // Retrieve the test statistic value used in the equivalence hypothesis
    local test_stat_value = .
    local test_stat_label ""
    local thresh_symbol ""
    
    if "`test_type'" == "max" {
        capture scalar _temp = e(max_abs_coef)
        if _rc == 0 {
            local test_stat_value = e(max_abs_coef)
        }
        local test_stat_label "max|{&beta}{sub:l}|"
        local thresh_symbol "{&delta}"
    }
    else if "`test_type'" == "mean" {
        capture scalar _temp = e(abs_mean_placebo)
        if _rc == 0 {
            local test_stat_value = e(abs_mean_placebo)
        }
        local test_stat_label "|{&beta}{sub:avg}|"
        local thresh_symbol "{&tau}"
    }
    else if "`test_type'" == "rms" {
        capture scalar _temp = e(rms_placebo_coefs)
        if _rc == 0 {
            local test_stat_value = e(rms_placebo_coefs)
        }
        local test_stat_label "{&beta}{sub:RMS}"
        local thresh_symbol "{&zeta}"
    }
    
    // Retrieve rejection decision (available only when threshold was specified)
    local reject_decision = .
    local thresh_spec_flag = e(threshold_specified)
    if `thresh_spec_flag' == 1 {
        capture scalar _temp = e(reject)
        if _rc == 0 {
            local reject_decision = e(reject)
        }
    }
    
    // -------------------------------------------------------------------------
    // Computation of relative time indices
    // -------------------------------------------------------------------------
    
    // Period information is extracted from e(preperiods) when available;
    // otherwise periods are parsed from e(placebo_coef_names)
    
    local n_preperiods : word count `preperiods'
    local placebo_idx = 0
    
    if `n_preperiods' > 0 {
        forvalues i = 1/`n_preperiods' {
            local period : word `i' of `preperiods'
            local rel_time = `period' - `base_period'
            
            // Base period (relative time = 0) is excluded from placebo estimates
            if `period' != `base_period' {
                local placebo_idx = `placebo_idx' + 1
                local rel_time_`placebo_idx' = `rel_time'
            }
        }
    }
    else {
        // Fallback: periods are parsed from coefficient names (format: placebo_X)
        local placebo_names "`e(placebo_coef_names)'"
        local n_placebo_names : word count `placebo_names'
        
        if `n_placebo_names' > 0 {
            forvalues i = 1/`n_placebo_names' {
                local pname : word `i' of `placebo_names'
                local period_str = subinstr("`pname'", "placebo_", "", 1)
                local period = real("`period_str'")
                
                if `period' != . {
                    local placebo_idx = `placebo_idx' + 1
                    local rel_time_`placebo_idx' = `period' - `base_period'
                }
            }
        }
        else {
            // Index-based relative time assignment when period information is unavailable
            di as text "Warning: Unable to determine exact pre-treatment periods."
            di as text "         Using index as relative time (may not be accurate)."
            forvalues i = 1/`no_placebos' {
                local placebo_idx = `placebo_idx' + 1
                local rel_time_`placebo_idx' = `i' - `no_placebos' - 1
            }
        }
    }
    
    // -------------------------------------------------------------------------
    // Computation of confidence intervals
    // -------------------------------------------------------------------------
    
    if "`ci'" != "" & `se_available' == 1 {
        // CI bounds: beta +/- z_{alpha/2} * SE
        local z_alpha = invnormal((100 + `level')/200)
        
        forvalues j = 1/`no_placebos' {
            local coef = `b_placebo'[`j', 1]
            local se = `se_placebo'[`j', 1]
            local ci_lower_`j' = `coef' - `z_alpha' * `se'
            local ci_upper_`j' = `coef' + `z_alpha' * `se'
        }
    }
    
    // -------------------------------------------------------------------------
    // Preparation of plot data
    // -------------------------------------------------------------------------
    
    preserve
    clear
    
    // Observations include placebo coefficients plus base period reference
    local total_points = `no_placebos' + 1
    quietly set obs `total_points'
    
    gen relative_time = .
    gen coef = .
    gen ci_lower = .
    gen ci_upper = .
    gen is_base = 0
    
    // Populate placebo coefficient data
    forvalues j = 1/`no_placebos' {
        quietly replace relative_time = `rel_time_`j'' in `j'
        quietly replace coef = `b_placebo'[`j', 1] in `j'
        
        if "`ci'" != "" & `se_available' == 1 {
            quietly replace ci_lower = `ci_lower_`j'' in `j'
            quietly replace ci_upper = `ci_upper_`j'' in `j'
        }
    }
    
    // Base period reference point: coefficient normalized to zero by construction
    local base_obs = `no_placebos' + 1
    quietly replace relative_time = 0 in `base_obs'
    quietly replace coef = 0 in `base_obs'
    quietly replace is_base = 1 in `base_obs'
    
    sort relative_time
    
    // -------------------------------------------------------------------------
    // Default style configuration
    // -------------------------------------------------------------------------
    
    // Point style defaults
    if "`msymbol'" == "" local msymbol "O"
    if "`msize'" == "" local msize "medium"
    if "`mcolor'" == "" local mcolor "navy"
    
    // Base period style defaults
    if "`basemsymbol'" == "" local basemsymbol "S"
    if "`basemcolor'" == "" local basemcolor "`mcolor'"
    
    // Confidence interval style defaults
    if "`cilcolor'" == "" local cilcolor "navy"
    if "`cilwidth'" == "" local cilwidth "medium"
    
    // Threshold line style defaults
    if "`threshlcolor'" == "" local threshlcolor "red"
    if "`threshlwidth'" == "" local threshlwidth "medium"
    if "`threshlpattern'" == "" local threshlpattern "dash"
    
    // Title defaults
    if `"`title'"' == "" local title "Pre-trend Placebo Coefficients"
    if `"`xtitle'"' == "" local xtitle "Period relative to treatment"
    if `"`ytitle'"' == "" local ytitle "Coefficient"
    
    // -------------------------------------------------------------------------
    // Auto-generated subtitle based on test type (when user does not override)
    // -------------------------------------------------------------------------
    
    if `"`subtitle'"' == "" & `has_threshold' == 1 {
        // Format threshold value for display
        local thresh_fmt : di %9.4f `plot_threshold'
        local thresh_fmt = strtrim("`thresh_fmt'")
        
        if "`test_type'" == "max" {
            local subtitle `"Equivalence bound: {&plusmn}`thresh_symbol' = {&plusmn}`thresh_fmt'  (all |{&beta}{sub:l}| < `thresh_symbol')"'
        }
        else if "`test_type'" == "mean" {
            local subtitle `"Equivalence bound: {&plusmn}`thresh_symbol' = {&plusmn}`thresh_fmt'  (`test_stat_label' < `thresh_symbol', conservative visual)"'
        }
        else if "`test_type'" == "rms" {
            local subtitle `"Equivalence bound: {&plusmn}`thresh_symbol' = {&plusmn}`thresh_fmt'  (`test_stat_label' < `thresh_symbol', conservative visual)"'
        }
    }
    
    // -------------------------------------------------------------------------
    // Auto-generated note with test result summary (when user does not override)
    // -------------------------------------------------------------------------
    
    // Note: Auto-generated notes are built as separate line locals (note_line1,
    // note_line2) and assembled directly into note_opt to avoid compound-quote
    // nesting issues. User-provided notes via note() option are handled separately.
    
    local auto_note = 0
    local note_line1 ""
    local note_line2 ""
    
    if `"`note'"' == "" {
        local auto_note = 1
        
        // Line 1: Test type and method identification
        if "`test_type'" == "max" {
            local method_label ""
            if "`test_method'" == "iu" {
                local method_label "Intersection-Union"
            }
            else if "`test_method'" == "boot" {
                local method_label "Bootstrap"
            }
            else if "`test_method'" == "wild" {
                local method_label "Wild Cluster Bootstrap"
            }
            local note_line1 "Test: Maximum (H{sub:0}: ||{&beta}||{sub:{&infin}} {&ge} `thresh_symbol')  Method: `method_label'"
        }
        else if "`test_type'" == "mean" {
            local note_line1 "Test: Mean (H{sub:0}: |{&beta}{sub:avg}| {&ge} `thresh_symbol')  Bound applies to average, not individual coefficients"
        }
        else if "`test_type'" == "rms" {
            local note_line1 "Test: RMS (H{sub:0}: {&beta}{sub:RMS} {&ge} `thresh_symbol')  Bound applies to root mean square, not individual coefficients"
        }
        
        // Line 2: Test statistic value and decision
        if `test_stat_value' != . {
            local stat_fmt : di %9.4f `test_stat_value'
            local stat_fmt = strtrim("`stat_fmt'")
            
            local alpha_fmt : di %4.2f e(alpha)
            local alpha_fmt = strtrim("`alpha_fmt'")
            
            if `has_threshold' == 1 {
                if `reject_decision' == 1 {
                    local note_line2 "`test_stat_label' = `stat_fmt'  Reject H{sub:0}: Yes (equivalence supported at {&alpha} = `alpha_fmt')"
                }
                else if `reject_decision' == 0 {
                    local note_line2 "`test_stat_label' = `stat_fmt'  Reject H{sub:0}: No"
                }
                else {
                    // min_threshold mode: no rejection decision available
                    local note_line2 "`test_stat_label' = `stat_fmt'  Min. equiv. threshold: `thresh_fmt'"
                }
            }
            else {
                local note_line2 "`test_stat_label' = `stat_fmt'"
            }
        }
    }
    
    // -------------------------------------------------------------------------
    // Construction of graph command
    // -------------------------------------------------------------------------
    
    local tw_cmd ""
    local plot_num = 0
    
    // Confidence intervals (if enabled)
    if "`ci'" != "" & `se_available' == 1 {
        local plot_num = `plot_num' + 1
        local tw_cmd "`tw_cmd' (rcap ci_lower ci_upper relative_time if is_base==0, lcolor(`cilcolor') lwidth(`cilwidth'))"
    }
    
    // Point estimates
    local plot_num = `plot_num' + 1
    if "`connect'" != "" {
        local tw_cmd "`tw_cmd' (connected coef relative_time if is_base==0, msymbol(`msymbol') msize(`msize') mcolor(`mcolor') lcolor(`mcolor'))"
    }
    else {
        local tw_cmd "`tw_cmd' (scatter coef relative_time if is_base==0, msymbol(`msymbol') msize(`msize') mcolor(`mcolor'))"
    }
    
    // Base period reference point
    if "`nobase'" == "" {
        local plot_num = `plot_num' + 1
        local tw_cmd "`tw_cmd' (scatter coef relative_time if is_base==1, msymbol(`basemsymbol') msize(`msize') mcolor(`basemcolor'))"
    }
    
    // Zero reference line
    local yline_opt ""
    if "`noline'" == "" {
        local yline_opt "yline(0, lcolor(gs10) lwidth(thin))"
    }
    
    // Equivalence threshold lines (+/- delta)
    local thresh_opt ""
    local yscale_opt ""
    if `has_threshold' == 1 {
        local neg_thresh = -`plot_threshold'
        local thresh_opt "yline(`plot_threshold', lcolor(`threshlcolor') lwidth(`threshlwidth') lpattern(`threshlpattern')) yline(`neg_thresh', lcolor(`threshlcolor') lwidth(`threshlwidth') lpattern(`threshlpattern'))"
        
        // Y-axis range is adjusted to ensure threshold bounds are visible
        quietly summarize coef
        local coef_min = r(min)
        local coef_max = r(max)
        
        if "`ci'" != "" & `se_available' == 1 {
            quietly summarize ci_lower
            if r(min) < `coef_min' {
                local coef_min = r(min)
            }
            quietly summarize ci_upper
            if r(max) > `coef_max' {
                local coef_max = r(max)
            }
        }
        
        // Padding is included around threshold bounds for visual clarity
        local y_min = min(`coef_min', `neg_thresh') - 0.05 * abs(`plot_threshold')
        local y_max = max(`coef_max', `plot_threshold') + 0.05 * abs(`plot_threshold')
        local yscale_opt "yscale(range(`y_min' `y_max'))"
    }
    
    // -------------------------------------------------------------------------
    // Assembly of graph options
    // -------------------------------------------------------------------------
    
    local title_opts `"title(`"`title'"')"'
    if `"`subtitle'"' != "" {
        local title_opts `"`title_opts' subtitle(`"`subtitle'"')"'
    }
    
    local axis_opts `"xtitle(`"`xtitle'"') ytitle(`"`ytitle'"')"'
    
    if `"`xlabel'"' != "" {
        local axis_opts `"`axis_opts' xlabel(`xlabel')"'
    }
    if `"`ylabel'"' != "" {
        local axis_opts `"`axis_opts' ylabel(`ylabel')"'
    }
    
    local legend_opt ""
    if `"`legend'"' != "" {
        local legend_opt `"legend(`legend')"'
    }
    else {
        local legend_opt "legend(off)"
    }
    
    local note_opt ""
    if `auto_note' == 1 {
        // Auto-generated note: build note_opt directly with proper quoting
        // Each line is a separate quoted string inside note()
        if `"`note_line2'"' != "" {
            local note_opt `"note("`note_line1'" "`note_line2'", size(vsmall))"'
        }
        else if `"`note_line1'"' != "" {
            local note_opt `"note("`note_line1'", size(vsmall))"'
        }
    }
    else if `"`note'"' != "" {
        // User-provided note
        local note_opt `"note(`"`note'"', size(vsmall))"'
    }
    
    local scheme_opt ""
    if "`scheme'" != "" {
        local scheme_opt "scheme(`scheme')"
    }
    
    local name_opt ""
    if "`name'" != "" {
        // name() is passed through as-is from user input (may include , replace)
        local name_opt `"name(`name')"'
    }
    
    // -------------------------------------------------------------------------
    // Graph execution
    // -------------------------------------------------------------------------
    
    local full_cmd `"twoway `tw_cmd', `yline_opt' `thresh_opt' `yscale_opt' `title_opts' `axis_opts' `legend_opt' `note_opt' `scheme_opt' `name_opt'"'
    
    `full_cmd'
    
    if `"`saving'"' != "" {
        if "`replace'" != "" {
            graph save `"`saving'"', replace
        }
        else {
            graph save `"`saving'"'
        }
    }
    
    restore
    
end
