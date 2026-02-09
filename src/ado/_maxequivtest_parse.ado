*! _maxequivtest_parse.ado - Parameter parsing for the maximum equivalence test
*!
*! This module parses and validates user-specified options for maxequivtest.
*! The maximum test examines hypothesis (3.1) from Dette & Schumann (2024):
*!   H0: ||beta||_inf >= delta  vs  H1: ||beta||_inf < delta
*! where ||beta||_inf = max_{l=1,...,T} |beta_l| is the sup-norm.
*!
*! Supported methods: IU (intersection-union), Boot (parametric bootstrap),
*! Wild (wild cluster bootstrap for non-spherical errors).

program define _maxequivtest_parse, rclass
    version 16.0
    
    syntax, Y(varname numeric) ID(varname numeric) G(varname numeric) Period(varname numeric) ///
        [X(varlist numeric)] [Cluster(varname numeric)] ///
        [PREtreatment(numlist)] [BASEperiod(real -999999)] ///
        [THRESHold(real -999999)] [Alpha(real 0.05)] ///
        [Method(string)] [Reps(integer 1000)] [Seed(integer 0)] ///
        [VCE(string)]
    
    // -------------------------------------------------------------------------
    // Common Validation
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
    // Method Specification
    // Three methods are available for the hypothesis (3.1):
    //   IU   - Intersection-union test based on folded normal quantiles
    //   Boot - Parametric bootstrap under spherical errors (Theorem 1)
    //   Wild - Wild cluster bootstrap for non-spherical errors (Remark 1c)
    // -------------------------------------------------------------------------
    
    if "`method'" == "" {
        local method "IU"
    }
    
    // Normalize method name to canonical form
    local method = upper("`method'")
    if "`method'" == "BOOT" {
        local method "Boot"
    }
    else if "`method'" == "WILD" {
        local method "Wild"
    }
    else if "`method'" == "IU" {
        local method "IU"
    }
    else {
        display as error "Invalid method: `method'"
        display as error "Valid methods are: IU, Boot, Wild"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Significance Level Validation
    // -------------------------------------------------------------------------
    
    if `alpha' <= 0 | `alpha' >= 1 {
        display as error "alpha must be between 0 and 1 (exclusive)"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Equivalence Threshold Validation
    // The threshold > 0 defines the boundary of negligible deviations
    // from the parallel trends assumption
    // -------------------------------------------------------------------------
    
    local threshold_specified = 0
    if `threshold' != -999999 {
        local threshold_specified = 1
        if `threshold' <= 0 {
            display as error "threshold must be strictly positive (> 0)"
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Bootstrap Replications
    // -------------------------------------------------------------------------
    
    if `reps' <= 0 {
        display as error "reps must be a positive integer"
        exit 198
    }
    
    if `reps' < 100 {
        display as text "Warning: reps < 100 may produce unreliable results"
    }
    
    // -------------------------------------------------------------------------
    // Random Seed
    // -------------------------------------------------------------------------
    
    if `seed' < 0 {
        display as error "seed must be non-negative"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Variance-Covariance Estimator Selection
    // Default is OLS (homoskedastic) variance estimator
    //
    // Custom Mata VCE types:
    //   ols:         Homoskedastic (sigma^2 * (X'X)^{-1})
    //   robust/hc1:  HC1 heteroskedasticity-robust (White 1980)
    //   hac:         HC3 Arellano panel-robust (for backward compatibility)
    //   cluster:     CR0 cluster-robust without adjustment (R package compatible)
    //   hc1_cluster: HC1 cluster-robust with small-sample adjustment
    //
    // Native Stata VCE types (via regress command):
    //   hc2:         HC2 leverage-adjusted (MacKinnon & White 1985)
    //   hc3:         HC3 more conservative leverage-adjusted
    //   cr0:         CR0 cluster-robust without adjustment
    //   cr1:         CR1 cluster-robust with DF adjustment (Stata default)
    // -------------------------------------------------------------------------
    
    local vce_type = "ols"
    if "`vce'" != "" {
        local vce_first : word 1 of `vce'
        local vce_second : word 2 of `vce'
        local vce_lower = lower("`vce_first'")
        
        // Map user-specified VCE type to internal representation
        if "`vce_lower'" == "ols" | "`vce_lower'" == "homoskedastic" | "`vce_lower'" == "hom" {
            // Classical OLS variance: V = sigma^2 * (X'X)^{-1}
            local vce_type = "ols"
        }
        else if "`vce_lower'" == "hc1" | "`vce_lower'" == "heteroskedastic" | "`vce_lower'" == "het" | "`vce_lower'" == "robust" | "`vce_lower'" == "hc" {
            // HC1: Heteroskedasticity-consistent with n/(n-k) adjustment (custom Mata)
            local vce_type = "robust"
        }
        else if "`vce_lower'" == "hc2" {
            // HC2: Leverage-adjusted heteroskedasticity-robust (native Stata)
            local vce_type = "hc2"
        }
        else if "`vce_lower'" == "hc3" {
            // HC3: More conservative leverage-adjusted (native Stata)
            local vce_type = "hc3"
        }
        else if "`vce_lower'" == "hac" {
            // HAC: Arellano panel-robust (custom Mata, for backward compatibility)
            local vce_type = "hac"
        }
        else if "`vce_lower'" == "hc1_cluster" | "`vce_lower'" == "hc1cluster" {
            // HC1 cluster-robust at the panel unit level
            local vce_type = "hc1_cluster"
            if "`vce_second'" != "" {
                local cluster "`vce_second'"
            }
            if "`cluster'" == "" {
                display as error "cluster variable required when vce(hc1_cluster) specified"
                display as error "Use either vce(hc1_cluster varname) or cluster(varname) with vce(hc1_cluster)"
                exit 198
            }
        }
        else if "`vce_lower'" == "cluster" | "`vce_lower'" == "cl" {
            // cluster: CR0 cluster-robust (R package compatible, custom Mata)
            local vce_type = "cluster"
            if "`vce_second'" != "" {
                local cluster "`vce_second'"
            }
            if "`cluster'" == "" {
                display as error "cluster variable required when vce(cluster) specified"
                display as error "Use either vce(cluster varname) or cluster(varname) with vce(cluster)"
                exit 198
            }
        }
        else if "`vce_lower'" == "cr0" {
            // CR0: Cluster-robust without adjustment (native Stata)
            local vce_type = "cr0"
            if "`vce_second'" != "" {
                local cluster "`vce_second'"
            }
            if "`cluster'" == "" {
                display as error "cluster variable required when vce(cr0) specified"
                display as error "Use either vce(cr0 varname) or cluster(varname) with vce(cr0)"
                exit 198
            }
        }
        else if "`vce_lower'" == "cr1" {
            // CR1: Cluster-robust with DF adjustment (native Stata, Stata default)
            local vce_type = "cr1"
            if "`vce_second'" != "" {
                local cluster "`vce_second'"
            }
            if "`cluster'" == "" {
                display as error "cluster variable required when vce(cr1) specified"
                display as error "Use either vce(cr1 varname) or cluster(varname) with vce(cr1)"
                exit 198
            }
        }
        else {
            display as error "Invalid vce option: `vce'"
            display as error "Valid options: ols, robust/hc1, hc2, hc3, hac, cluster, cr0, cr1, hc1_cluster"
            exit 198
        }
    }
    
    // Bootstrap methods compute variance internally; VCE option is not applicable
    if ("`method'" == "Boot" | "`method'" == "Wild") & "`vce'" != "" {
        display as text "Warning: vce() option is ignored for Bootstrap methods"
    }
    
    // -------------------------------------------------------------------------
    // Return Parsed Values
    // -------------------------------------------------------------------------
    
    return local method "`method'"
    return local vce_type "`vce_type'"
    
    if inlist("`vce_type'", "cluster", "cr0", "cr1", "hc1_cluster") & "`cluster'" != "" {
        return local cluster "`cluster'"
    }
    else {
        return local cluster ""
    }
    
    return scalar alpha = `alpha'
    return scalar reps = `reps'
    return scalar seed = `seed'
    return scalar threshold_specified = `threshold_specified'
    return scalar threshold = `threshold'
    
end
