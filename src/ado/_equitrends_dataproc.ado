*! _equitrends_dataproc.ado - Panel data processing for equivalence testing
*!
*! Prepares panel data for pre-trend equivalence tests by performing
*! data binding, missing value handling, pretreatment period filtering,
*! base period selection, and placebo variable construction.

program define _equitrends_dataproc, rclass
    version 16.0
    
    syntax, Y(varname numeric) ID(varname numeric) G(varname numeric) Period(varname numeric) ///
        [X(varlist numeric)] [Cluster(varname numeric)] ///
        [PREtreatment(numlist)] [BASEperiod(real -999999)]
    
    * Note: Stata's syntax command stores parameter values in lowercase locals
    * So Y(varname) stores the value in local `y', not `Y'
    
    * -------------------------------------------------------------------------
    * Data loading
    * -------------------------------------------------------------------------
    
    * Get data vectors
    tempname Y_vec ID_vec G_vec period_vec X_mat cluster_vec pretreat_vec
    
    mata: st_matrix("`Y_vec'", st_data(., "`y'"))
    mata: st_matrix("`ID_vec'", st_data(., "`id'"))
    mata: st_matrix("`G_vec'", st_data(., "`g'"))
    mata: st_matrix("`period_vec'", st_data(., "`period'"))
    
    * Handle X variables (preserve names for diagnostic messages)
    local has_X = 0
    local orig_x_names ""
    if "`x'" != "" {
        local has_X = 1
        local orig_x_names "`x'"
        mata: st_matrix("`X_mat'", st_data(., "`x'"))
    }
    
    * Handle cluster variable
    local has_cluster = 0
    if "`cluster'" != "" {
        local has_cluster = 1
        mata: st_matrix("`cluster_vec'", st_data(., "`cluster'"))
    }
    
    * Handle pretreatment periods
    local has_pretreat = 0
    if "`pretreatment'" != "" {
        local has_pretreat = 1
        local n_pretreat : word count `pretreatment'
        matrix `pretreat_vec' = J(`n_pretreat', 1, .)
        local i = 1
        foreach val of local pretreatment {
            matrix `pretreat_vec'[`i', 1] = `val'
            local i = `i' + 1
        }
    }
    
    * Handle base period (default sentinel value indicates auto-selection)
    local has_baseperiod = 0
    if `baseperiod' != -999999 {
        local has_baseperiod = 1
    }
    
    * -------------------------------------------------------------------------
    * Data processing pipeline
    * -------------------------------------------------------------------------
    
    mata: _equitrends_run_dataproc( ///
        st_matrix("`Y_vec'"), ///
        st_matrix("`ID_vec'"), ///
        st_matrix("`G_vec'"), ///
        st_matrix("`period_vec'"), ///
        `has_X', ///
        `has_X' ? st_matrix("`X_mat'") : J(0, 0, .), ///
        `has_cluster', ///
        `has_cluster' ? st_matrix("`cluster_vec'") : J(0, 1, .), ///
        `has_pretreat', ///
        `has_pretreat' ? st_matrix("`pretreat_vec'") : J(0, 1, .), ///
        `has_baseperiod', ///
        `baseperiod', ///
        "`orig_x_names'" ///
    )
    
    * -------------------------------------------------------------------------
    * Result extraction
    * -------------------------------------------------------------------------
    return scalar baseperiod = `_dp_baseperiod'
    return scalar balanced_panel = `_dp_is_balanced'
    return scalar no_placebos = `_dp_no_placebos'
    return scalar n = `_dp_n'
    return scalar n_t = `_dp_n_t'
    return scalar N_obs = `_dp_N_obs'
    return scalar na_omitted = `_dp_na_omitted'
    
    return local orig_names = "`_dp_orig_names'"
    return local placebo_names = "`_dp_placebo_names'"
    return local no_periods = "`_dp_no_periods'"
    
    * Display warning if NA rows were omitted
    if `_dp_na_omitted' > 0 {
        di as text "Warning: `_dp_na_omitted' rows of pre-treatment data omitted due to NAs."
    }
    
end

// ============================================================================
// _equitrends_run_dataproc()
// Bridge function between Stata ado-file and Mata data processor class.
//
// Instantiates EquitrendsDataProcessor, executes the processing pipeline,
// and stores results in Stata locals for return by the calling program.
// ============================================================================

mata:

void _equitrends_run_dataproc(
    real colvector Y,
    real colvector ID,
    real colvector G,
    real colvector period,
    real scalar has_X,
    real matrix X,
    real scalar has_cluster,
    real colvector cluster,
    real scalar has_pretreat,
    real colvector pretreat,
    real scalar has_baseperiod,
    real scalar baseperiod,
    string scalar orig_x_names_str
)
{
    class EquitrendsDataProcessor scalar dp
    string vector placebo_names, orig_names, orig_x_names
    real rowvector no_periods_info
    string scalar placebo_str, orig_str, no_periods_str
    real scalar i, bp_final
    
    // Parse covariate names for diagnostic output
    if (orig_x_names_str != "") {
        orig_x_names = tokens(orig_x_names_str)
    }
    else {
        orig_x_names = J(1, 0, "")
    }
    
    // Initialize data processor and load data
    dp = EquitrendsDataProcessor()
    bp_final = has_baseperiod ? baseperiod : .
    dp.load_data(Y, ID, G, period, X, cluster, pretreat, bp_final)
    
    // Set covariate names for diagnostic output
    if (cols(orig_x_names) > 0) {
        dp.set_orig_x_names(orig_x_names)
    }
    
    // Run processing pipeline
    dp.process()
    
    // Store results in Stata locals
    st_local("_dp_baseperiod", strofreal(dp.get_base_period()))
    st_local("_dp_is_balanced", strofreal(dp.get_is_balanced()))
    st_local("_dp_no_placebos", strofreal(dp.get_no_placebos()))
    st_local("_dp_n", strofreal(dp.get_n()))
    st_local("_dp_n_t", strofreal(dp.get_n_t()))
    st_local("_dp_N_obs", strofreal(dp.get_N_obs()))
    st_local("_dp_na_omitted", strofreal(dp.get_na_omitted()))
    
    // Build placebo names string
    placebo_names = dp.get_placebo_names()
    placebo_str = ""
    for (i = 1; i <= cols(placebo_names); i++) {
        if (i > 1) placebo_str = placebo_str + " "
        placebo_str = placebo_str + placebo_names[i]
    }
    st_local("_dp_placebo_names", placebo_str)
    
    // Build original X names string
    orig_names = dp.get_orig_X_names()
    orig_str = ""
    for (i = 1; i <= cols(orig_names); i++) {
        if (i > 1) orig_str = orig_str + " "
        orig_str = orig_str + orig_names[i]
    }
    st_local("_dp_orig_names", orig_str)
    
    // Build no_periods string
    no_periods_info = dp.get_no_periods_info()
    if (cols(no_periods_info) == 1) {
        no_periods_str = strofreal(no_periods_info[1])
    }
    else {
        no_periods_str = strofreal(no_periods_info[1]) + " " + strofreal(no_periods_info[2])
    }
    st_local("_dp_no_periods", no_periods_str)
}

end

