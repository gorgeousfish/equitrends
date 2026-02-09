*! _equitrends_validate_common.ado - Common input validation for equivalence tests
*!
*! This module provides centralized validation logic shared across all equivalence
*! test commands (maxequivtest, meanequivtest, rmsequivtest, equivtest).
*!
*! Validates:
*!   - Treatment indicator G is binary (0/1)
*!   - Pretreatment period is a valid subset of the time variable
*!   - At least two unique pretreatment periods exist
*!   - Base period is within the pretreatment period
*!
*! Syntax:
*!   _equitrends_validate_common, G(varname) Period(varname) ///
*!       [PREtreatment(numlist)] [BASEperiod(real)]
*!
*! Exit codes:
*!   0   - Validation passed
*!   198 - Validation failed (with error message)

program define _equitrends_validate_common
    version 16.0
    
    syntax, G(varname numeric) Period(varname numeric) ///
        [PREtreatment(numlist)] [BASEperiod(real -999999)]
    
    // -------------------------------------------------------------------------
    // Validate treatment indicator G is binary (0/1)
    // -------------------------------------------------------------------------
    
    tempname G_vec
    mata: st_matrix("`G_vec'", st_data(., "`g'"))
    mata: st_local("valid_binary", strofreal(_equitrends_is_binary(st_matrix("`G_vec'"))))
    if "`valid_binary'" == "0" {
        display as error "Entries of G must either be logical (e.g. TRUE/FALSE) or binary (e.g. 0/1)."
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Validate pretreatment period and base period
    // - Pretreatment period must be a subset of the time variable
    // - Base period must be within the pretreatment period
    // - At least two unique pretreatment periods are required
    // -------------------------------------------------------------------------
    
    // Get unique periods from data for validation
    quietly levelsof `period', local(data_periods)
    
    if "`pretreatment'" != "" {
        // Check pretreatment periods exist in data
        foreach p of numlist `pretreatment' {
            local p_found = 0
            foreach dp of local data_periods {
                if abs(`p' - `dp') < 1e-10 {
                    local p_found = 1
                }
            }
            if `p_found' == 0 {
                display as error "pretreatment_period must be a subset of period"
                display as error "  Period `p' not found in data"
                exit 198
            }
        }
        
        // Parse pretreatment periods to check count
        local pretreat_list "`pretreatment'"
        local pretreat_unique : list uniq pretreat_list
        local n_pretreat : word count `pretreat_unique'
        
        // Require at least two unique pretreatment periods
        if `n_pretreat' < 2 {
            display as error "pre-treatment period must have at least two unique periods."
            exit 198
        }
        
        // Validate base period is within pretreatment
        if `baseperiod' != -999999 {
            local base_in_pretreat = 0
            foreach p of numlist `pretreat_list' {
                if abs(`p' - `baseperiod') < 1e-10 {
                    local base_in_pretreat = 1
                }
            }
            if `base_in_pretreat' == 0 {
                display as error "base_period must be an element of pretreatment_period."
                display as error "  baseperiod = `baseperiod'"
                display as error "  pretreatment = `pretreatment'"
                exit 198
            }
        }
    }
    else {
        // Default: use all data periods as pretreatment
        // Validate base period is in data periods
        if `baseperiod' != -999999 {
            local base_in_data = 0
            foreach dp of local data_periods {
                if abs(`baseperiod' - `dp') < 1e-10 {
                    local base_in_data = 1
                }
            }
            if `base_in_data' == 0 {
                display as error "base_period must be an element of pretreatment_period."
                display as error "  baseperiod(`baseperiod') not found in data periods"
                exit 198
            }
            display as text "Note: baseperiod(`baseperiod') specified without pretreatment()."
            display as text "      All data periods will be used as pretreatment."
        }
        
        // Require at least two unique periods
        local n_data_periods : word count `data_periods'
        if `n_data_periods' < 2 {
            display as error "pre-treatment period must have at least two unique periods."
            exit 198
        }
    }
    
end
