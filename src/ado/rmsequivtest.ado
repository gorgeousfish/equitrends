*! rmsequivtest.ado - RMS equivalence test for pre-trends in DiD estimation
*!
*! Implements the root mean square (RMS) equivalence test for assessing
*! the parallel trends assumption in difference-in-differences estimation.
*! A self-normalized inference procedure is used that does not require
*! explicit variance estimation.
*!
*! Hypothesis:
*!   H0: beta_RMS >= zeta  vs  H1: beta_RMS < zeta
*!   where beta_RMS = ||beta|| / sqrt(T) = sqrt((1/T) * sum(beta_l^2))
*!
*! Rejection criterion:
*!   Reject H0 if beta_RMS^2 < zeta^2 + Q_W(alpha) * V_n
*!   where Q_W(alpha) is the alpha-quantile of the limiting W distribution
*!   and V_n is the self-normalized variance estimator computed from
*!   sequential subsamples.

program define rmsequivtest, eclass
    version 16.0
    
    // Load required Mata library
    equitrends_init
    
    local cmdline "rmsequivtest `0'"
    
    // -------------------------------------------------------------------------
    // Syntax Parsing
    // -------------------------------------------------------------------------
    
    // Parse command options
    // Note: Group() accepts g/gr/gro/grou/group; Time() accepts t/ti/tim/time
    // Period() is an alias for Time() for backward compatibility
    syntax varname(numeric), ID(varname numeric) ///
        [Group(varname numeric) G(varname numeric)] ///
        [Time(varname numeric) Period(varname numeric)] ///
        [X(varlist numeric)] [THRESHold(real -999999)] ///
        [PREtreatment(numlist)] [BASEperiod(real -999999)] ///
        [NOLambda(integer 5)] [Alpha(real 0.05)] [LEVel(integer 95)] ///
        [SEED(integer -999999)] [NODISplay]
    
    local y `varlist'
    
    // Handle Group/G alias
    if "`group'" != "" & "`g'" != "" {
        display as error "may not specify both group() and g()"
        exit 198
    }
    if "`group'" == "" & "`g'" == "" {
        display as error "option group() or g() required"
        exit 198
    }
    if "`group'" != "" {
        local g `group'
    }
    
    // Handle Time/Period alias
    if "`time'" != "" & "`period'" != "" {
        display as error "may not specify both time() and period()"
        exit 198
    }
    if "`time'" == "" & "`period'" == "" {
        display as error "option time() or period() required"
        exit 198
    }
    if "`time'" != "" {
        local period `time'
    }
    
    local id_var `id'
    local g_var `g'
    local period_var `period'
    
    // -------------------------------------------------------------------------
    // Parameter Validation
    // -------------------------------------------------------------------------
    local parse_cmd "_rmsequivtest_parse, y(`y') id(`id_var') g(`g_var') period(`period_var')"
    local parse_cmd "`parse_cmd' alpha(`alpha') nolambda(`nolambda') level(`level')"
    if `threshold' != -999999 {
        local parse_cmd "`parse_cmd' threshold(`threshold')"
    }
    if "`pretreatment'" != "" {
        local parse_cmd "`parse_cmd' pretreatment(`pretreatment')"
    }
    if `baseperiod' != -999999 {
        local parse_cmd "`parse_cmd' baseperiod(`baseperiod')"
    }
    
    quietly `parse_cmd'
    
    local threshold_specified = r(threshold_specified)
    local threshold_val = r(threshold)
    local alpha_val = r(alpha)
    local level_val = r(level)
    
    // -------------------------------------------------------------------------
    // Data Processing
    // -------------------------------------------------------------------------
    local dataproc_cmd "_equitrends_dataproc, y(`y') id(`id_var') g(`g_var') period(`period_var')"
    if "`x'" != "" {
        local dataproc_cmd "`dataproc_cmd' x(`x')"
    }
    if "`pretreatment'" != "" {
        local dataproc_cmd "`dataproc_cmd' pretreatment(`pretreatment')"
    }
    if `baseperiod' != -999999 {
        local dataproc_cmd "`dataproc_cmd' baseperiod(`baseperiod')"
    }
    
    quietly `dataproc_cmd'
    
    local baseperiod_val = r(baseperiod)
    local is_balanced = r(balanced_panel)
    local no_placebos = r(no_placebos)
    local n_individuals = r(n)
    local n_periods = r(n_t)
    local N_obs = r(N_obs)
    local dp_placebo_names "`r(placebo_names)'"
    local dp_no_periods "`r(no_periods)'"
    
    // Extract T_min and T_max for unbalanced panels
    local dp_t_min = -999999
    local dp_t_max = -999999
    if `is_balanced' == 0 {
        local dp_t_min : word 1 of `dp_no_periods'
        local dp_t_max : word 2 of `dp_no_periods'
    }
    
    // -------------------------------------------------------------------------
    // Data Validation
    // -------------------------------------------------------------------------
    
    if `n_individuals' == 0 {
        display as error "No valid individuals found in data."
        display as error "Check that ID variable contains non-missing values."
        exit 2000
    }
    
    // Check sample size requirement for sequential subsampling
    if `n_individuals' < `nolambda' {
        display as error "Number of individuals (N=`n_individuals') must be at least nolambda (`nolambda')."
        display as error "Please reduce nolambda or use data with more individuals."
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // RMS Test Computation
    // -------------------------------------------------------------------------
    local core_cmd "_rmsequivtest_core, y(`y') id(`id_var') g(`g_var') period(`period_var')"
    local core_cmd "`core_cmd' noplacebos(`no_placebos')"
    local core_cmd "`core_cmd' thresholdspecified(`threshold_specified')"
    local core_cmd "`core_cmd' threshold(`threshold_val')"
    local core_cmd "`core_cmd' alpha(`alpha_val') nolambda(`nolambda') level(`level_val')"
    
    if "`x'" != "" {
        local core_cmd "`core_cmd' x(`x')"
    }
    if "`pretreatment'" != "" {
        local core_cmd "`core_cmd' preperiod(`pretreatment')"
    }
    if `baseperiod_val' != -999999 {
        local core_cmd "`core_cmd' baseperiod(`baseperiod_val')"
    }
    if `seed' != -999999 {
        local core_cmd "`core_cmd' seed(`seed')"
    }
    
    `core_cmd'
    
    local MS_placebo = r(MS_placebo)
    local RMS_placebo = r(RMS_placebo)
    local V_n = r(V_n)
    local Q_W = r(Q_W)
    
    // CI results
    local RMS_ci_lower = r(RMS_ci_lower)
    local RMS_ci_upper = r(RMS_ci_upper)
    local MS_ci_lower = r(MS_ci_lower)
    local MS_ci_upper = r(MS_ci_upper)
    local level_ret = r(level)
    
    tempname b_placebo MS_lambda
    matrix `b_placebo' = r(placebo_coefs)
    matrix `MS_lambda' = r(MS_lambda)
    
    if `threshold_specified' == 1 {
        local MS_critical = r(MS_critical)
        local RMS_critical = r(RMS_critical)
        local reject = r(reject)
    }
    else {
        local min_threshold = r(min_threshold)
    }
    
    // -------------------------------------------------------------------------
    // Display Results (unless nodisplay option specified)
    // -------------------------------------------------------------------------
    if "`nodisplay'" == "" {
        local disp_cmd "_rmsequivtest_display, rms_placebo(`RMS_placebo')"
        local disp_cmd "`disp_cmd' threshold_specified(`threshold_specified')"
        local disp_cmd "`disp_cmd' alpha(`alpha_val') n(`N_obs') n_g(`n_individuals')"
        local disp_cmd "`disp_cmd' n_t(`n_periods') no_placebos(`no_placebos')"
        local disp_cmd "`disp_cmd' baseperiod(`baseperiod_val') balanced(`is_balanced')"
        local disp_cmd "`disp_cmd' t_min(`dp_t_min') t_max(`dp_t_max')"
        local disp_cmd "`disp_cmd' rms_ci_lower(`RMS_ci_lower') rms_ci_upper(`RMS_ci_upper')"
        local disp_cmd "`disp_cmd' ms_ci_lower(`MS_ci_lower') ms_ci_upper(`MS_ci_upper')"
        local disp_cmd "`disp_cmd' level(`level_ret')"
        
        if `threshold_specified' == 1 {
            local disp_cmd "`disp_cmd' threshold(`threshold_val')"
            local disp_cmd "`disp_cmd' rms_critical(`RMS_critical') reject(`reject')"
        }
        else {
            local disp_cmd "`disp_cmd' min_threshold(`min_threshold')"
        }
        
        if `seed' != -999999 {
            local disp_cmd "`disp_cmd' seed(`seed')"
        }
        
        `disp_cmd'
    }
    
    // -------------------------------------------------------------------------
    // Store Estimation Results
    // -------------------------------------------------------------------------
    
    ereturn clear
    
    // Matrices
    ereturn matrix b_placebo = `b_placebo'
    ereturn matrix MS_lambda = `MS_lambda'
    
    // Scalars
    ereturn scalar N = `N_obs'
    ereturn scalar N_g = `n_individuals'
    ereturn scalar T = `n_periods'
    ereturn scalar no_placebos = `no_placebos'
    ereturn scalar MS_placebo = `MS_placebo'
    ereturn scalar RMS_placebo = `RMS_placebo'
    ereturn scalar V_n = `V_n'
    ereturn scalar Q_W = `Q_W'
    ereturn scalar alpha = `alpha_val'
    ereturn scalar nolambda = `nolambda'
    ereturn scalar base_period = `baseperiod_val'
    ereturn scalar is_balanced = `is_balanced'
    
    // Confidence interval results
    ereturn scalar RMS_ci_lower = `RMS_ci_lower'
    ereturn scalar RMS_ci_upper = `RMS_ci_upper'
    ereturn scalar MS_ci_lower = `MS_ci_lower'
    ereturn scalar MS_ci_upper = `MS_ci_upper'
    ereturn scalar level = `level_ret'
    
    // Store T_min and T_max for unbalanced panels
    if `is_balanced' == 0 {
        ereturn scalar T_min = `dp_t_min'
        ereturn scalar T_max = `dp_t_max'
    }
    
    if `seed' != -999999 {
        ereturn scalar seed = `seed'
    }
    else {
        ereturn scalar seed = .
    }
    
    if `threshold_specified' == 1 {
        ereturn scalar threshold_specified = 1
        ereturn scalar threshold = `threshold_val'
        ereturn scalar MS_critical = `MS_critical'
        ereturn scalar RMS_critical = `RMS_critical'
        ereturn scalar reject = `reject'
    }
    else {
        ereturn scalar threshold_specified = 0
        ereturn scalar min_threshold = `min_threshold'
    }
    
    // Macros
    ereturn local cmd "rmsequivtest"
    ereturn local cmdline "`cmdline'"
    ereturn local depvar "`y'"
    ereturn local idvar "`id_var'"
    ereturn local groupvar "`g_var'"
    ereturn local timevar "`period_var'"
    ereturn local properties "b"
    
    ereturn local placebo_names "`dp_placebo_names'"
    
end
