*! _rmsequivtest_parse.ado - Parameter parsing for RMS equivalence test
*!
*! Parses and validates user-supplied arguments for the rmsequivtest command.
*! This module implements input validation for the root mean square (RMS)
*! equivalence test, which tests H_0: beta_RMS >= zeta vs H_1: beta_RMS < zeta.
*! Common validations are delegated to _equitrends_validate_common.ado.

program define _rmsequivtest_parse, rclass
    version 16.0
    
    syntax, Y(varname numeric) ID(varname numeric) G(varname numeric) Period(varname numeric) ///
        Alpha(real) NOLambda(integer) ///
        [THRESHold(real -999999)] [PREtreatment(numlist)] [BASEperiod(real -999999)] ///
        [LEVel(real 95)]
    
    // -------------------------------------------------------------------------
    // Common validation
    // -------------------------------------------------------------------------
    
    local validate_cmd "_equitrends_validate_common, g(`g') period(`period')"
    if "`pretreatment'" != "" {
        local validate_cmd "`validate_cmd' pretreatment(`pretreatment')"
    }
    if `baseperiod' != -999999 {
        local validate_cmd "`validate_cmd' baseperiod(`baseperiod')"
    }
    `validate_cmd'
    
    // -------------------------------------------------------------------------
    // Validate alpha
    // The self-normalized statistic converges to a pivotal limiting
    // distribution W. Pre-computed quantiles Q_W(alpha) are available for
    // alpha in {0.01, 0.025, 0.05, 0.1, 0.2}.
    // -------------------------------------------------------------------------
    
    local valid_alpha "0.01 0.025 0.05 0.1 0.2"
    local alpha_ok = 0
    local tol = 1e-10
    local alpha_normalized = `alpha'
    
    foreach a of local valid_alpha {
        if abs(`alpha' - `a') < `tol' {
            local alpha_ok = 1
            local alpha_normalized = `a'
        }
    }
    
    if !`alpha_ok' {
        display as error "alpha must be one of 0.01, 0.025, 0.05, 0.1 or 0.2"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Validate nolambda
    // Number of support points for the discrete measure nu used in the
    // self-normalization variance estimator V_n. A minimum of 2 points
    // is required for non-degenerate variance estimation.
    // -------------------------------------------------------------------------
    
    if `nolambda' <= 0 | `nolambda' != round(`nolambda') {
        display as error "no_lambda must be a positive integer"
        exit 198
    }
    
    if `nolambda' < 2 {
        display as error "nolambda must be at least 2 (otherwise V_n cannot be computed)"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Validate threshold
    // The equivalence threshold zeta > 0 defines the null hypothesis
    // H_0: beta_RMS >= zeta vs H_1: beta_RMS < zeta.
    // -------------------------------------------------------------------------
    
    local threshold_specified = 0
    local threshold_val = -999999
    
    if `threshold' != -999999 {
        if `threshold' <= 0 {
            display as error "equiv_threshold must be positive"
            exit 198
        }
        local threshold_specified = 1
        local threshold_val = `threshold'
    }
    
    // -------------------------------------------------------------------------
    // Validate level
    // Confidence level for confidence intervals. Must be one of the supported
    // values determined by available W distribution quantiles.
    // Supported levels: 80, 90, 95, 98, 99 (corresponding to α = 0.2, 0.1, 0.05, 0.02, 0.01)
    // -------------------------------------------------------------------------
    
    if `level' <= 0 | `level' >= 100 {
        display as error "level must be between 0 and 100 (exclusive)"
        exit 198
    }
    
    // Check if level is supported (requires both lower and upper quantiles)
    // Supported levels: 80, 90, 95, 98, 99
    local valid_levels "80 90 95 98 99"
    local level_ok = 0
    local level_normalized = `level'
    
    foreach l of local valid_levels {
        if abs(`level' - `l') < `tol' {
            local level_ok = 1
            local level_normalized = `l'
        }
    }
    
    if !`level_ok' {
        display as error "level must be one of 80, 90, 95, 98 or 99"
        display as error "  (Other levels require additional W distribution quantiles)"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Return results
    // -------------------------------------------------------------------------
    
    return scalar threshold_specified = `threshold_specified'
    return scalar threshold = `threshold_val'
    return scalar alpha = `alpha_normalized'
    return scalar level = `level_normalized'
    
end
