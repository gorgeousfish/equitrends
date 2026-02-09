*! _equitrends_qfoldnorm.ado - Folded normal quantile function interface
*!
*! Syntax: _equitrends_qfoldnorm p mean sd, result(name)
*!
*! Parameters:
*!   p      - Probability value, 0 < p < 1
*!   mean   - Mean of the underlying normal distribution μ
*!   sd     - Standard deviation of the underlying normal distribution σ > 0
*!   result - Name of scalar to store the result
*!
*! Returns:
*!   r(result) - The p-th quantile of the folded normal distribution
*!
*! Exit codes:
*!   0   - Success
*!   198 - Parameter error
*!
*! Mathematical definition:
*!   If X ~ N(μ, σ²), then |X| follows a folded normal distribution.
*!   This function computes the p-th quantile of |X|.
*!
*! Implementation:
*!   Uses Mata function _eqt_qfoldnorm() for numerical computation.

capture program drop _equitrends_qfoldnorm
program define _equitrends_qfoldnorm, rclass
    version 16.0
    
    * ========================================
    * Layer 1: Parameter parsing and validation
    * ========================================
    syntax anything(name=args), Result(name)
    
    * Extract parameters using tokenize
    tokenize `args'
    local p `1'
    local mean `2'
    local sd `3'
    
    * Validate all three parameters are provided
    if "`p'" == "" | "`mean'" == "" | "`sd'" == "" {
        display as error "Error: _equitrends_qfoldnorm requires 3 arguments: p mean sd"
        display as error "Usage: _equitrends_qfoldnorm p mean sd, result(name)"
        exit 198
    }
    
    * Validate p is numeric
    capture confirm number `p'
    if _rc != 0 {
        display as error "Error: p must be a numeric value, got '`p''"
        exit 198
    }
    
    * Validate mean is numeric
    capture confirm number `mean'
    if _rc != 0 {
        display as error "Error: mean must be a numeric value, got '`mean''"
        exit 198
    }
    
    * Validate sd is numeric
    capture confirm number `sd'
    if _rc != 0 {
        display as error "Error: sd must be a numeric value, got '`sd''"
        exit 198
    }
    
    * Validate p ∈ (0, 1)
    if `p' <= 0 | `p' >= 1 {
        display as error "Error: p must be strictly between 0 and 1, got `p'"
        exit 198
    }
    
    * Validate sd > 0
    if `sd' <= 0 {
        display as error "Error: sd must be strictly positive, got `sd'"
        exit 198
    }
    
    * ========================================
    * Layer 2: Call Mata function
    * ========================================
    tempname mata_result
    
    mata: st_numscalar("`mata_result'", _eqt_qfoldnorm(`p', `mean', `sd'))
    
    * Check for missing result (indicates Mata error)
    if scalar(`mata_result') == . {
        display as error "Error: Mata computation failed for qfoldnorm"
        exit 198
    }
    
    * ========================================
    * Store and return result
    * ========================================
    scalar `result' = scalar(`mata_result')
    return scalar `result' = scalar(`mata_result')
end
