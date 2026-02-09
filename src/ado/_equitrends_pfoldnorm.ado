*! _equitrends_pfoldnorm.ado - Folded normal CDF interface
*!
*! Computes the cumulative distribution function (CDF) of the folded normal
*! distribution. If X ~ N(μ, σ²), then |X| follows N_F(μ, σ²).
*!
*! Syntax: _equitrends_pfoldnorm x mean sd, result(name)
*!
*! Arguments:
*!   x      - Quantile value (x >= 0)
*!   mean   - Location parameter of the underlying normal distribution
*!   sd     - Scale parameter (sd > 0)
*!   result - Name of scalar to store the result
*!
*! Returns:
*!   Scalar containing P(|X| <= x) where X ~ N(mean, sd^2)
*!
*! Mathematical formula:
*!   F(x; mu, sigma) = Phi((x-mu)/sigma) + Phi((x+mu)/sigma) - 1

capture program drop _equitrends_pfoldnorm
program define _equitrends_pfoldnorm, rclass
    version 16.0
    
    // -------------------------------------------------------------------------
    // Parameter Parsing and Validation
    // -------------------------------------------------------------------------
    syntax anything(name=args), Result(name)
    
    // Extract parameters using tokenize
    tokenize `args'
    local x `1'
    local mean `2'
    local sd `3'
    
    // Validate all three parameters are provided
    if "`x'" == "" | "`mean'" == "" | "`sd'" == "" {
        display as error "Error: _equitrends_pfoldnorm requires 3 arguments: x mean sd"
        display as error "Usage: _equitrends_pfoldnorm x mean sd, result(name)"
        exit 198
    }
    
    // Validate x is numeric
    capture confirm number `x'
    if _rc != 0 {
        display as error "Error: x must be a numeric value, got '`x''"
        exit 198
    }
    
    // Validate mean is numeric
    capture confirm number `mean'
    if _rc != 0 {
        display as error "Error: mean must be a numeric value, got '`mean''"
        exit 198
    }
    
    // Validate sd is numeric
    capture confirm number `sd'
    if _rc != 0 {
        display as error "Error: sd must be a numeric value, got '`sd''"
        exit 198
    }
    
    // Validate x >= 0 (folded normal has support [0, ∞))
    if `x' < 0 {
        display as error "Error: x must be non-negative for folded normal CDF, got `x'"
        exit 198
    }
    
    // Validate sd > 0
    if `sd' <= 0 {
        display as error "Error: sd must be strictly positive, got `sd'"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Call Mata Function
    // -------------------------------------------------------------------------
    tempname mata_result
    
    mata: st_numscalar("`mata_result'", _eqt_pfoldnorm(`x', `mean', `sd'))
    
    // Check for missing result (indicates Mata error)
    if scalar(`mata_result') == . {
        display as error "Error: Mata computation failed for pfoldnorm"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Store and Return Result
    // -------------------------------------------------------------------------
    scalar `result' = scalar(`mata_result')
    return scalar `result' = scalar(`mata_result')
end
