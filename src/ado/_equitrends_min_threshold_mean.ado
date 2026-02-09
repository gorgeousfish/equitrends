*! _equitrends_min_threshold_mean.ado - Compute minimum threshold for mean equivalence test
*!
*! Computes the smallest equivalence threshold tau such that the null hypothesis
*! H_0: |beta_bar| >= tau can be rejected at significance level alpha, based on
*! the folded normal distribution. This implements the equivalence test for the
*! average placebo treatment effect described in Section 4.2.2 of the paper.

capture program drop _equitrends_min_threshold_mean
program define _equitrends_min_threshold_mean, rclass
    version 16.0
    
    // -------------------------------------------------------------------------
    // Parameter parsing and validation
    // -------------------------------------------------------------------------
    syntax anything(name=args), Result(name)
    
    // Extract positional arguments: coef, sd, alpha
    tokenize `args'
    local coef `1'
    local sd `2'
    local alpha `3'
    
    // Validate all three parameters are provided
    if "`coef'" == "" | "`sd'" == "" | "`alpha'" == "" {
        display as error "Error: _equitrends_min_threshold_mean requires 3 arguments: coef sd alpha"
        display as error "Usage: _equitrends_min_threshold_mean coef sd alpha, result(name)"
        exit 198
    }
    
    // Validate coef is numeric
    capture confirm number `coef'
    if _rc != 0 {
        display as error "Error: coef must be a numeric value, got '`coef''"
        exit 198
    }
    
    // Validate sd is numeric
    capture confirm number `sd'
    if _rc != 0 {
        display as error "Error: sd must be a numeric value, got '`sd''"
        exit 198
    }
    
    // Validate alpha is numeric
    capture confirm number `alpha'
    if _rc != 0 {
        display as error "Error: alpha must be a numeric value, got '`alpha''"
        exit 198
    }
    
    // Validate coef >= 0 (absolute value of mean coefficient)
    if `coef' < 0 {
        display as error "Error: coef must be non-negative, got `coef'"
        exit 198
    }
    
    // Validate sd > 0 (standard error must be positive)
    if `sd' <= 0 {
        display as error "Error: sd must be strictly positive, got `sd'"
        exit 198
    }
    
    // Validate alpha in (0, 1) (significance level)
    if `alpha' <= 0 | `alpha' >= 1 {
        display as error "Error: alpha must be strictly between 0 and 1, got `alpha'"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Compute minimum threshold via Mata
    // -------------------------------------------------------------------------
    // Minimum threshold tau: smallest value such that pfoldnorm(|beta_bar|, tau, sigma) = alpha
    // where pfoldnorm is the CDF of the folded normal distribution
    tempname mata_result
    
    mata: st_numscalar("`mata_result'", _eqt_mean_min_threshold(`coef', `sd', `alpha'))
    
    // Check for missing result (indicates numerical computation failure)
    if scalar(`mata_result') == . {
        display as error "Error: Mata computation failed for min_threshold_mean"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Store and return result
    // -------------------------------------------------------------------------
    scalar `result' = scalar(`mata_result')
    return scalar `result' = scalar(`mata_result')
end
