*! _equivtest_parse.ado - Parameter parsing and validation for equivtest
*!
*! Parses and validates all parameters for the unified equivtest command.
*! Supports three equivalence test types corresponding to the hypotheses in
*! Dette and Schumann (2024, JBES):
*!   - max:  Test H0: ||beta||_inf >= delta (hypothesis 3.1)
*!   - mean: Test H0: |beta_bar| >= tau (hypothesis 3.2)
*!   - rms:  Test H0: beta_RMS >= zeta (hypothesis 3.3)
*!
*! Type-specific options are validated to ensure appropriate usage.
*! Common validation is delegated to _equitrends_validate_common.ado.

program define _equivtest_parse, rclass
    version 16.0
    
    syntax, TYPE(string) Period(varname numeric) G(varname numeric) ///
        [Method(string)] [VCE(string)] [Alpha(real 0.05)] ///
        [Nboot(integer 1000)] [NOLambda(integer 5)] ///
        [Cluster(varname numeric)] [THRESHold(real -999999)] ///
        [NBOOT_specified(integer 0)] [NOLAMBDA_specified(integer 0)] ///
        [PREtreatment(numlist)] [BASEperiod(real -999999)]
    
    // -------------------------------------------------------------------------
    // Common Validation
    // -------------------------------------------------------------------------
    // Validates treatment group indicator G (binary) and pre-treatment period
    // specification. Delegated to _equitrends_validate_common.ado.
    
    local validate_cmd "_equitrends_validate_common, g(`g') period(`period')"
    if "`pretreatment'" != "" {
        local validate_cmd "`validate_cmd' pretreatment(`pretreatment')"
    }
    if `baseperiod' != -999999 {
        local validate_cmd "`validate_cmd' baseperiod(`baseperiod')"
    }
    `validate_cmd'
    
    // -------------------------------------------------------------------------
    // Test Type Validation
    // -------------------------------------------------------------------------
    // Validates the equivalence test type:
    //   max:  Maximum test for H0: ||beta||_inf >= delta
    //   mean: Mean test for H0: |beta_bar| >= tau
    //   rms:  Root mean square test for H0: beta_RMS >= zeta
    
    local type = lower("`type'")
    if !inlist("`type'", "max", "mean", "rms") {
        display as error "type must be max, mean, or rms"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Method Validation (type=max only)
    // -------------------------------------------------------------------------
    // For the maximum test (hypothesis 3.1), three methods are available:
    //   iu:   Intersection-union test using folded normal quantiles
    //   boot: Parametric bootstrap with constrained estimator
    //   wild: Wild cluster bootstrap for non-spherical errors
    // Default method is 'iu' for computational efficiency.
    
    local method_validated ""
    if "`method'" != "" {
        if "`type'" != "max" {
            display as error "method() option only allowed with type(max)"
            exit 198
        }
        local method_validated = lower("`method'")
        if !inlist("`method_validated'", "iu", "boot", "wild") {
            display as error "method must be iu, boot, or wild"
            exit 198
        }
    }
    else if "`type'" == "max" {
        local method_validated "iu"
    }
    
    // -------------------------------------------------------------------------
    // Variance-Covariance Estimator Validation
    // -------------------------------------------------------------------------
    // Validates the variance-covariance estimator for Sigma (Assumption 2).
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
    //   hc3:         HC3 more conservative leverage-adjusted (Davidson & MacKinnon 1993)
    //   cr0:         CR0 cluster-robust without adjustment
    //   cr1:         CR1 cluster-robust with DF adjustment (Stata default)
    //
    // Note: RMS test uses self-normalization and does not require Sigma
    // estimation. Bootstrap methods handle variance internally.
    
    local vce_validated ""
    if "`vce'" != "" {
        // RMS test is self-normalized and does not require VCE specification
        if "`type'" == "rms" {
            display as error "vce() option not allowed with type(rms)"
            exit 198
        }
        // Bootstrap methods compute variance internally
        if "`type'" == "max" & inlist("`method_validated'", "boot", "wild") {
            display as error "vce() option not allowed with bootstrap methods"
            exit 198
        }
        
        local vce_lower = lower("`vce'")
        
        // OLS: homoskedastic variance
        if inlist("`vce_lower'", "ols", "homoskedastic", "hom") {
            local vce_validated = "ols"
        }
        // HC1: heteroskedasticity-robust with small-sample adjustment (custom Mata)
        else if inlist("`vce_lower'", "hc1", "robust", "heteroskedastic", "het", "hc") {
            local vce_validated = "robust"
        }
        // HC2: leverage-adjusted heteroskedasticity-robust (native Stata)
        else if "`vce_lower'" == "hc2" {
            local vce_validated = "hc2"
        }
        // HC3: more conservative leverage-adjusted (native Stata)
        else if "`vce_lower'" == "hc3" {
            local vce_validated = "hc3"
        }
        // HAC: Arellano panel-robust (custom Mata, for backward compatibility)
        else if "`vce_lower'" == "hac" {
            local vce_validated = "hac"
        }
        // CR0: cluster-robust without adjustment (native or custom based on alias)
        else if inlist("`vce_lower'", "cr0") {
            local vce_validated = "cr0"
        }
        // cluster: CR0 cluster-robust (R package compatible, custom Mata)
        else if inlist("`vce_lower'", "cluster", "cl") {
            local vce_validated = "cluster"
        }
        // CR1: cluster-robust with DF adjustment (native Stata)
        else if "`vce_lower'" == "cr1" {
            local vce_validated = "cr1"
        }
        // HC1 cluster-robust: cluster-robust with small-sample adjustment
        else if inlist("`vce_lower'", "hc1_cluster", "hc1cluster") {
            local vce_validated = "hc1_cluster"
        }
        else {
            display as error "vce must be: ols, robust/hc1, hc2, hc3, hac, cluster, cr0, cr1, or hc1_cluster"
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Cluster Variable Validation
    // -------------------------------------------------------------------------
    // Ensures cluster variable is specified if and only if cluster-robust
    // variance estimation is requested.
    
    if "`cluster'" != "" {
        if "`type'" == "rms" {
            display as error "cluster() option not allowed with type(rms)"
            exit 198
        }
        if !inlist("`vce_validated'", "cluster", "cr0", "cr1", "hc1_cluster") {
            display as error "cluster() requires vce(cluster), vce(cr0), vce(cr1), or vce(hc1_cluster)"
            exit 198
        }
    }
    else if inlist("`vce_validated'", "cluster", "cr0", "cr1", "hc1_cluster") {
        display as error "cluster() required with vce(cluster), vce(cr0), vce(cr1), or vce(hc1_cluster)"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Significance Level Validation
    // -------------------------------------------------------------------------
    // Alpha must be in (0, 1). For RMS test, only specific alpha values are
    // supported due to pre-computed critical values for the self-normalized
    // statistic W (Theorem 2).
    
    if `alpha' <= 0 | `alpha' >= 1 {
        display as error "alpha must lie between 0 and 1"
        exit 198
    }
    if "`type'" == "rms" {
        // RMS test requires pre-computed critical values for the limiting
        // distribution of W (equation 4.16)
        local alpha_valid = 0
        foreach a in 0.01 0.025 0.05 0.1 0.2 {
            if abs(`alpha' - `a') < 1e-10 {
                local alpha_valid = 1
            }
        }
        if `alpha_valid' == 0 {
            display as error "alpha must be one of 0.01, 0.025, 0.05, 0.1 or 0.2"
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Lambda Grid Size Validation (type=rms only)
    // -------------------------------------------------------------------------
    // The RMS test uses a discrete approximation to compute V_n (equation 4.15).
    // At least 2 lambda values are required to compute the variance estimator.
    
    if `nolambda_specified' == 1 {
        if "`type'" != "rms" {
            display as error "nolambda() option only allowed with type(rms)"
            exit 198
        }
    }
    if "`type'" == "rms" {
        if `nolambda' <= 0 | `nolambda' != round(`nolambda') {
            display as error "no_lambda must be a positive integer"
            exit 198
        }
        // V_n computation requires at least two lambda values
        if `nolambda' < 2 {
            display as error "nolambda must be at least 2 (otherwise V_n cannot be computed)"
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Bootstrap Replications Validation (type=max with method=boot|wild)
    // -------------------------------------------------------------------------
    // Number of bootstrap replications B for the maximum test bootstrap
    // procedures (Section 4.2.1).
    
    if `nboot_specified' == 1 {
        if "`type'" != "max" | !inlist("`method_validated'", "boot", "wild") {
            display as error "nboot() option only allowed with method(boot) or method(wild)"
            exit 198
        }
    }
    if "`type'" == "max" & inlist("`method_validated'", "boot", "wild") {
        if `nboot' <= 0 {
            display as error "B must be a strictly positive integer scalar"
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Equivalence Threshold Validation
    // -------------------------------------------------------------------------
    // The equivalence threshold (delta, tau, or zeta depending on test type)
    // must be strictly positive.
    
    local threshold_specified = 0
    if `threshold' != -999999 {
        local threshold_specified = 1
        if `threshold' <= 0 {
            display as error "equiv_threshold must be strictly positive (> 0)"
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Return Validated Parameters
    // -------------------------------------------------------------------------
    
    return local type "`type'"
    return local method "`method_validated'"
    return local vce "`vce_validated'"
    return local cluster "`cluster'"
    return scalar alpha = `alpha'
    return scalar nboot = `nboot'
    return scalar nolambda = `nolambda'
    return scalar threshold_specified = `threshold_specified'
    return scalar threshold = `threshold'
    
end
