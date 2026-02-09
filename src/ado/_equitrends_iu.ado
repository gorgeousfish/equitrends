*! _equitrends_iu.ado - Intersection-union test for maximum placebo coefficient
*!
*! Implements the intersection-union (IU) equivalence test for the maximum
*! absolute placebo coefficient. Tests the hypothesis:
*!     H0: ||beta||_inf >= delta  vs  H1: ||beta||_inf < delta
*! where ||beta||_inf = max_{l=1,...,T} |beta_l|.
*!
*! The IU test rejects H0 when all individual coefficients satisfy:
*!     |beta_l| < Q_alpha(delta, SE_l^2)  for all l = 1,...,T
*! where Q_alpha denotes the alpha-quantile of the folded normal distribution
*! N_F(delta, SE_l^2).
*!
*! The test is computationally efficient but tends to be conservative for T > 1,
*! as is characteristic of intersection-union procedures.

version 16.0

program define _equitrends_iu, rclass
    
    syntax [, Betas(string) SEs(string) Alpha(real 0.05) THreshold(real -999)]
    
    // =========================================================================
    // Input Validation
    // =========================================================================
    
    if "`betas'" == "" {
        display as error "betas() is required"
        exit 198
    }
    if "`ses'" == "" {
        display as error "ses() is required"
        exit 198
    }
    
    if `alpha' <= 0 | `alpha' >= 1 {
        display as error "alpha must be between 0 and 1"
        exit 198
    }
    
    capture confirm matrix `betas'
    if _rc {
        display as error "matrix `betas' not found"
        exit 111
    }
    
    capture confirm matrix `ses'
    if _rc {
        display as error "matrix `ses' not found"
        exit 111
    }
    
    local n_beta = rowsof(`betas')
    local n_se = rowsof(`ses')
    
    if `n_beta' != `n_se' {
        display as error "betas and ses must have same number of rows"
        exit 503
    }
    
    if `n_beta' == 0 {
        display as error "input matrices cannot be empty"
        exit 198
    }
    
    // =========================================================================
    // IU Test Execution
    // =========================================================================
    // 
    // Two operational modes:
    //   1. With threshold: Test H0 at specified delta, return rejection decision
    //   2. Without threshold: Compute minimum delta* for rejection at level alpha
    //
    // The minimum threshold delta* represents the smallest equivalence bound
    // for which the null hypothesis of non-negligible trend differences can
    // still be rejected at significance level alpha.
    
    if `threshold' != -999 {
        // ---------------------------------------------------------------------
        // Mode 1: Test at specified equivalence threshold
        // ---------------------------------------------------------------------
        if `threshold' <= 0 {
            display as error "equiv_threshold must be strictly positive (> 0)"
            exit 198
        }
        mata: _eqt_ado_iu_with_threshold("`betas'", "`ses'", `alpha', `threshold')
        return scalar reject = r(reject)
        return scalar equiv_threshold = r(equiv_threshold)
        return scalar alpha = r(alpha)
        return scalar n_coefs = r(n_coefs)
        return matrix critical_values = r(critical_values)
    }
    else {
        // ---------------------------------------------------------------------
        // Mode 2: Compute minimum equivalence threshold for rejection
        // ---------------------------------------------------------------------
        mata: _eqt_ado_iu_min_threshold("`betas'", "`ses'", `alpha')
        return scalar min_equiv_threshold = r(min_equiv_threshold)
        return scalar alpha = r(alpha)
        return scalar n_coefs = r(n_coefs)
        return matrix min_thresholds = r(min_thresholds)
    }
    
end
