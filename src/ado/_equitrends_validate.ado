*! _equitrends_validate.ado - Input validation for equivalence tests
*!
*! Validates all input parameters for the equivalence testing commands.
*! Performs three-phase validation: data structure, parameter values,
*! and test-specific constraints.

program define _equitrends_validate, rclass
    version 16.0
    
    syntax, Y(varname numeric) ID(varname numeric) G(varname numeric) Period(varname numeric) ///
        [X(varlist numeric)] [Cluster(varname numeric)] ///
        [PREtreatment(numlist)] [BASEperiod(real -999999)] ///
        [EQUIVthreshold(real 1)] [Alpha(real 0.05)] ///
        [Type(string)] [Method(string)] [B(integer 1000)]
    
    // Initialize return values
    return local error = "0"
    return local message = ""
    
    // Set defaults
    if "`type'" == "" {
        local type "max"
    }
    if "`method'" == "" {
        local method "IU"
    }
    
    // =========================================================================
    // Phase 1: Data Structure Validation
    // =========================================================================
    
    // Verify required variables are numeric
    capture confirm numeric variable `y'
    if _rc != 0 {
        return local error = "1"
        return local message = "Y must be a numeric variable."
        di as error "Y must be a numeric variable."
        exit 198
    }
    
    capture confirm numeric variable `id'
    if _rc != 0 {
        return local error = "1"
        return local message = "ID must be a numeric variable. Use 'encode' to convert string IDs."
        di as error "ID must be a numeric variable. Use 'encode' to convert string IDs."
        exit 198
    }
    
    capture confirm numeric variable `g'
    if _rc != 0 {
        return local error = "1"
        return local message = "G must be a numeric variable."
        di as error "G must be a numeric variable."
        exit 198
    }
    
    capture confirm numeric variable `period'
    if _rc != 0 {
        return local error = "1"
        return local message = "period must be a numeric variable."
        di as error "period must be a numeric variable."
        exit 198
    }
    
    // Validate cluster variable if provided
    if "`cluster'" != "" {
        capture confirm numeric variable `cluster'
        if _rc != 0 {
            return local error = "1"
            return local message = "cluster variable must be numeric. Use 'encode' to convert string cluster variables."
            di as error "cluster variable must be numeric. Use 'encode' to convert string cluster variables."
            exit 198
        }
    }
    
    // Validate covariate variables if provided
    if "`x'" != "" {
        foreach var of local x {
            capture confirm numeric variable `var'
            if _rc != 0 {
                return local error = "1"
                return local message = "X variable `var' must be numeric."
                di as error "X variable `var' must be numeric."
                exit 198
            }
        }
    }
    
    // Verify data contains observations
    quietly count
    if r(N) == 0 {
        return local error = "1"
        return local message = "No observations in data."
        di as error "Error: No observations in data."
        exit 2000
    }
    
    if r(N) == 1 {
        di as text "Warning: Only 1 observation in data."
    }
    
    // =========================================================================
    // Phase 2: Parameter Value Validation
    // =========================================================================
    
    // Equivalence threshold must be strictly positive: delta > 0
    mata: st_local("valid_thresh", strofreal(_equitrends_val_equiv_thresh(`equivthreshold')))
    if "`valid_thresh'" == "0" {
        return local error = "1"
        return local message = "equiv_threshold must be strictly positive (>0)."
        di as error "equiv_threshold must be strictly positive (>0)."
        exit 198
    }
    
    // Significance level must satisfy 0 < alpha < 1
    mata: st_local("valid_alpha", strofreal(_equitrends_validate_alpha_range(`alpha')))
    if "`valid_alpha'" == "0" {
        return local error = "1"
        return local message = "alpha must lie between 0 and 1."
        di as error "alpha must lie between 0 and 1."
        exit 198
    }
    
    // Treatment indicator must be binary (0/1)
    tempname G_vec
    mata: st_matrix("`G_vec'", st_data(., "`g'"))
    mata: st_local("valid_binary", strofreal(_equitrends_is_binary(st_matrix("`G_vec'"))))
    if "`valid_binary'" == "0" {
        return local error = "1"
        return local message = "Entries of G must either be logical (e.g. TRUE/FALSE) or binary (e.g. 0/1)."
        di as error "Entries of G must either be logical (e.g. TRUE/FALSE) or binary (e.g. 0/1)."
        exit 198
    }
    
    // Pre-treatment periods must be a subset of observed periods
    if "`pretreatment'" != "" {
        tempname period_vec pretreat_vec
        mata: st_matrix("`period_vec'", uniqrows(st_data(., "`period'")))
        
        local n_pretreat : word count `pretreatment'
        matrix `pretreat_vec' = J(`n_pretreat', 1, .)
        local i = 1
        foreach val of local pretreatment {
            matrix `pretreat_vec'[`i', 1] = `val'
            local i = `i' + 1
        }
        
        mata: st_local("valid_subset", strofreal(_equitrends_is_subset(st_matrix("`pretreat_vec'"), st_matrix("`period_vec'"))))
        if "`valid_subset'" == "0" {
            return local error = "1"
            return local message = "pretreatment_period must be a subset of period"
            di as error "pretreatment_period must be a subset of period"
            exit 198
        }
        
        // At least two unique pre-treatment periods are required
        local n_unique_pretreat : word count `pretreatment'
        local pretreat_unique : list uniq pretreatment
        local n_unique_pretreat : word count `pretreat_unique'
        if `n_unique_pretreat' < 2 {
            return local error = "1"
            return local message = "pretreatment_period must contain at least 2 unique periods."
            di as error "pretreatment_period must contain at least 2 unique periods."
            exit 198
        }
    }
    else {
        // Default: all periods used as pre-treatment; require at least two
        tempname period_vec_check
        mata: st_matrix("`period_vec_check'", uniqrows(st_data(., "`period'")))
        local n_unique_period = rowsof(`period_vec_check')
        if `n_unique_period' < 2 {
            return local error = "1"
            return local message = "pre-treatment period must have at least two unique periods."
            di as error "pre-treatment period must have at least two unique periods."
            exit 198
        }
    }
    
    // Base period must be an element of the pre-treatment periods
    if `baseperiod' != -999999 & "`pretreatment'" != "" {
        local base_found = 0
        foreach val of local pretreatment {
            if abs(`baseperiod' - `val') < 1e-10 {
                local base_found = 1
            }
        }
        if `base_found' == 0 {
            return local error = "1"
            return local message = "base_period must be an element of pretreatment_period."
            di as error "base_period must be an element of pretreatment_period."
            exit 198
        }
    }
    
    // =========================================================================
    // Phase 3: Test-Specific Validation
    // =========================================================================
    
    // Test type: max (hypothesis 3.1), mean (hypothesis 3.2), rms (hypothesis 3.3)
    mata: st_local("valid_type", strofreal(_equitrends_validate_type("`type'")))
    if "`valid_type'" == "0" {
        return local error = "1"
        return local message = "type must be max, mean, or rms"
        di as error "type must be max, mean, or rms"
        exit 198
    }
    
    // Method validation for maximum test: IU, Boot, or Wild bootstrap
    if "`type'" == "max" {
        mata: st_local("valid_method", strofreal(_equitrends_validate_method("`method'")))
        if "`valid_method'" == "0" {
            return local error = "1"
            return local message = "method must be IU, Boot, or Wild"
            di as error "method must be IU, Boot, or Wild"
            exit 198
        }
    }
    
    // RMS test requires specific alpha values for pre-computed critical values
    if "`type'" == "rms" {
        local valid_rms_alpha = 0
        foreach a in 0.01 0.025 0.05 0.1 0.2 {
            if abs(`alpha' - `a') < 1e-10 {
                local valid_rms_alpha = 1
            }
        }
        if `valid_rms_alpha' == 0 {
            return local error = "1"
            return local message = "For RMS equivalence test, alpha must be one of: 0.01, 0.025, 0.05, 0.1, 0.2"
            di as error "For RMS equivalence test, alpha must be one of: 0.01, 0.025, 0.05, 0.1, 0.2"
            exit 198
        }
    }
    
    // Bootstrap replications must be a positive integer
    if "`method'" == "Boot" | "`method'" == "Wild" {
        mata: st_local("valid_B", strofreal(_equitrends_validate_B(`b')))
        if "`valid_B'" == "0" {
            return local error = "1"
            return local message = "B must be a strictly positive integer scalar"
            di as error "B must be a strictly positive integer scalar"
            exit 198
        }
    }
    
    // =========================================================================
    // Validation Passed
    // =========================================================================
    
    return local error = "0"
    return local message = ""
    
end
