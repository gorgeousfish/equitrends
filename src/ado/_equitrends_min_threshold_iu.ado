*! _equitrends_min_threshold_iu.ado - Minimum equivalence threshold for IU test
*!
*! Computes the minimum equivalence threshold delta such that the IU test
*! rejects at significance level alpha. Solves pfoldnorm(|beta|, delta, se) = alpha.
*!
*! Arguments:
*!   coef   : real - Absolute coefficient value |beta| >= 0
*!   sd     : real - Standard error of coefficient, > 0
*!   alpha  : real - Significance level, in (0, 1)
*!   result : name - Scalar name to store result
*!
*! Returns:
*!   Scalar containing minimum threshold delta
*!
*! Implementation:
*!   Delegates numerical computation to Mata function _eqt_iu_min_threshold_single()

capture program drop _equitrends_min_threshold_iu
program define _equitrends_min_threshold_iu, rclass
    version 16.0
    
    // -------------------------------------------------------------------------
    // Parameter Parsing and Validation
    // -------------------------------------------------------------------------
    syntax anything(name=args), Result(name)
    
    // Parse positional arguments
    tokenize `args'
    local coef `1'
    local sd `2'
    local alpha `3'
    
    // Validate all three parameters are provided
    if "`coef'" == "" | "`sd'" == "" | "`alpha'" == "" {
        display as error "Error: _equitrends_min_threshold_iu requires 3 arguments: coef sd alpha"
        display as error "Usage: _equitrends_min_threshold_iu coef sd alpha, result(name)"
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
    
    // Validate coef >= 0
    if `coef' < 0 {
        display as error "Error: coef must be non-negative, got `coef'"
        exit 198
    }
    
    // Validate sd > 0
    if `sd' <= 0 {
        display as error "Error: sd must be strictly positive, got `sd'"
        exit 198
    }
    
    // Validate alpha in (0, 1)
    if `alpha' <= 0 | `alpha' >= 1 {
        display as error "Error: alpha must be strictly between 0 and 1, got `alpha'"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Mata Function Call
    // -------------------------------------------------------------------------
    tempname mata_result
    
    mata: st_numscalar("`mata_result'", _eqt_iu_min_threshold_single(`coef', `sd', `alpha'))
    
    // Check for missing result (indicates Mata error)
    if scalar(`mata_result') == . {
        display as error "Error: Mata computation failed for min_threshold_iu"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Store and Return Result
    // -------------------------------------------------------------------------
    scalar `result' = scalar(`mata_result')
    return scalar `result' = scalar(`mata_result')
end
