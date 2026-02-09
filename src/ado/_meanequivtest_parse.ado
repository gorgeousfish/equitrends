*! _meanequivtest_parse.ado - Parameter parsing for mean equivalence test
*!
*! Parses and validates user-specified parameters for the meanequivtest command.
*! Implements input validation for the mean placebo effect test (hypothesis 3.2)
*! which tests H0: |β̄| ≥ τ vs H1: |β̄| < τ.
*!
*! Arguments:
*!   y       : varname  - Outcome variable
*!   id      : varname  - Panel unit identifier
*!   g       : varname  - Treatment group indicator
*!   period  : varname  - Time period variable
*!   x       : varlist  - Optional covariates
*!   cluster : varname  - Optional clustering variable
*!   pretreatment : numlist - Pre-treatment periods to include
*!   baseperiod   : real    - Reference period (omitted category)
*!   threshold    : real    - Equivalence threshold τ > 0
*!   alpha        : real    - Significance level (0, 1)
*!   vce          : string  - Variance-covariance estimator specification
*!
*! Returns:
*!   r(threshold_specified) : scalar - 1 if threshold specified, 0 otherwise
*!   r(threshold)           : scalar - Equivalence threshold (or -999999)
*!   r(alpha)               : scalar - Significance level
*!   r(vce_type)            : local  - VCE type identifier

version 16.0

program define _meanequivtest_parse, rclass
    
    syntax, Y(varname numeric) ID(varname numeric) G(varname numeric) Period(varname numeric) ///
        [X(varlist numeric)] [Cluster(varname numeric)] ///
        [PREtreatment(numlist)] [BASEperiod(real -999999)] ///
        [THRESHold(real -999999)] [Alpha(real 0.05)] ///
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
    // Significance Level Validation
    // -------------------------------------------------------------------------
    
    if `alpha' <= 0 | `alpha' >= 1 {
        display as error "Error: alpha must be strictly between 0 and 1"
        display as error "  Got: alpha = `alpha'"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Equivalence Threshold Validation
    // -------------------------------------------------------------------------
    
    local threshold_specified = 0
    if `threshold' != -999999 {
        local threshold_specified = 1
        if `threshold' <= 0 {
            display as error "Error: equiv_threshold must be strictly positive (> 0)"
            display as error "  Got: threshold = `threshold'"
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Variance-Covariance Estimator Specification
    // Default: classical homoskedastic (OLS) standard errors
    // -------------------------------------------------------------------------
    
    local vce_type = "ols"
    
    if "`vce'" != "" {
        local vce_first : word 1 of `vce'
        local vce_second : word 2 of `vce'
        local vce_lower = lower("`vce_first'")
        
        // OLS: classical homoskedastic standard errors
        if "`vce_lower'" == "ols" | "`vce_lower'" == "homoskedastic" | "`vce_lower'" == "hom" {
            local vce_type = "ols"
        }
        // HC1: heteroskedasticity-robust with N/(N-K) adjustment (custom Mata)
        else if "`vce_lower'" == "hc1" | "`vce_lower'" == "heteroskedastic" | "`vce_lower'" == "het" | "`vce_lower'" == "robust" | "`vce_lower'" == "hc" {
            local vce_type = "robust"
        }
        // HC2: leverage-adjusted heteroskedasticity-robust (native Stata)
        else if "`vce_lower'" == "hc2" {
            local vce_type = "hc2"
        }
        // HC3: more conservative leverage-adjusted (native Stata)
        else if "`vce_lower'" == "hc3" {
            local vce_type = "hc3"
        }
        // HAC: Arellano panel-robust (custom Mata, for backward compatibility)
        else if "`vce_lower'" == "hac" {
            local vce_type = "hac"
        }
        // HC1 with clustering on specified variable
        else if "`vce_lower'" == "hc1_cluster" | "`vce_lower'" == "hc1cluster" {
            local vce_type = "hc1_cluster"
            if "`vce_second'" != "" {
                local cluster "`vce_second'"
            }
            if "`cluster'" == "" {
                display as error "Error: cluster() option required when vce(hc1_cluster) specified"
                display as error "Use either vce(hc1_cluster varname) or cluster(varname) with vce(hc1_cluster)"
                exit 198
            }
        }
        // cluster: CR0 cluster-robust (R package compatible, custom Mata)
        else if "`vce_lower'" == "cluster" | "`vce_lower'" == "cl" {
            local vce_type = "cluster"
            if "`vce_second'" != "" {
                local cluster "`vce_second'"
            }
            if "`cluster'" == "" {
                display as error "Error: cluster() option required when vce(cluster) specified"
                exit 198
            }
        }
        // CR0: cluster-robust without adjustment (native Stata)
        else if "`vce_lower'" == "cr0" {
            local vce_type = "cr0"
            if "`vce_second'" != "" {
                local cluster "`vce_second'"
            }
            if "`cluster'" == "" {
                display as error "Error: cluster() option required when vce(cr0) specified"
                exit 198
            }
        }
        // CR1: cluster-robust with DF adjustment (native Stata, Stata default)
        else if "`vce_lower'" == "cr1" {
            local vce_type = "cr1"
            if "`vce_second'" != "" {
                local cluster "`vce_second'"
            }
            if "`cluster'" == "" {
                display as error "Error: cluster() option required when vce(cr1) specified"
                exit 198
            }
        }
        else {
            display as error "Error: invalid vce() specification: `vce'"
            display as error "  Valid options: ols, robust/hc1, hc2, hc3, hac, cluster, cr0, cr1, hc1_cluster"
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Cluster Variable Validation
    // -------------------------------------------------------------------------
    
    if "`cluster'" != "" {
        capture confirm variable `cluster'
        if _rc != 0 {
            display as error "Error: cluster variable '`cluster'' not found"
            exit 111
        }
        
        // Automatically select cluster-robust VCE when cluster variable is provided
        if !inlist("`vce_type'", "cluster", "cr0", "cr1", "hc1_cluster") {
            display as text "Note: cluster() specified, using vce(cluster)"
            local vce_type = "cluster"
        }
    }
    
    // -------------------------------------------------------------------------
    // Return Values
    // -------------------------------------------------------------------------
    
    return local vce_type "`vce_type'"
    
    if inlist("`vce_type'", "cluster", "cr0", "cr1", "hc1_cluster") & "`cluster'" != "" {
        return local cluster "`cluster'"
    }
    else {
        return local cluster ""
    }
    
    return scalar alpha = `alpha'
    return scalar threshold_specified = `threshold_specified'
    return scalar threshold = `threshold'
    
end
