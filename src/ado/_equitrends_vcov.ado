*! _equitrends_vcov.ado - Variance-covariance matrix computation interface
*!
*! This module provides a Stata interface to Mata VCE functions, implementing
*! variance-covariance estimators for the placebo coefficient vector beta.
*! Under Assumption 2, a consistent estimator of Sigma is required for inference
*! on the equivalence hypotheses (3.1), (3.2), and (3.3).
*!
*! Syntax:
*!   _equitrends_vcov, x_ddm(name) y_ddm(name) beta(name) id(name) ///
*!                     [vce(string) cluster(name) n(integer) n_t(integer) ///
*!                      no_placebos(integer) vcov_result(name) se_result(name)]
*!
*! VCE Types (Custom Mata Implementation):
*!   ols     - Homoskedastic VCE: Sigma = sigma^2 * (X'X)^{-1}
*!             where sigma^2 = RSS / (n*T - k)
*!   robust  - HC1 heteroskedasticity-consistent VCE (alias: hc1)
*!             Implements White's (1980) robust standard errors
*!   hac     - HC3 panel-robust VCE following Arellano (1987)
*!             Accounts for within-individual correlation structure
*!   cluster - CR0 cluster-robust VCE without small-sample adjustment (alias: cl)
*!             Standard clustered errors without finite-sample correction
*!
*! VCE Types (Native Stata Implementation via _equitrends_vce_native):
*!   hc2     - HC2 leverage-adjusted heteroskedasticity-robust VCE
*!             V_HC2 = (X'X)^{-1} * sum[e_i^2/(1-h_ii) * x_i*x_i'] * (X'X)^{-1}
*!   hc3     - HC3 more conservative leverage adjustment (Davidson & MacKinnon 1993)
*!             V_HC3 = (X'X)^{-1} * sum[e_i^2/(1-h_ii)^2 * x_i*x_i'] * (X'X)^{-1}
*!   cr0     - CR0 cluster-robust without small-sample adjustment
*!             V_CR0 = V_CR1 * (G-1)/G
*!   cr1     - CR1 cluster-robust with DF adjustment (Stata default)
*!             V_CR1 = [G/(G-1)] * (X'X)^{-1} * sum_g[X_g'e_g e_g'X_g] * (X'X)^{-1}
*!
*! Note: The cluster option implements CR0 (unadjusted) clustering to maintain
*! consistency with the R package implementation and theoretical framework.
*! For CR1 (Stata default), use vce(cr1).

program define _equitrends_vcov, rclass
    version 16.0
    
    syntax, x_ddm(name) y_ddm(name) beta(name) id(name) ///
            [vce(string) cluster(name) n(integer 0) n_t(integer 0) ///
             no_placebos(integer 0) vcov_result(name) se_result(name)]
    
    // -------------------------------------------------------------------------
    // Input Validation and Normalization
    // -------------------------------------------------------------------------
    
    // OLS VCE is used by default for homoskedastic panel data
    if "`vce'" == "" {
        local vce "ols"
    }
    
    // VCE type is normalized to lowercase for consistent matching
    local vce = lower("`vce'")
    
    // VCE type must be one of the supported estimators
    // Custom Mata: ols, robust/hc1, hac, cluster/cl
    // Native Stata: hc2, hc3, cr0, cr1
    if !inlist("`vce'", "ols", "robust", "hc1", "hac", "cluster", "cl", "hc2", "hc3", "cr0", "cr1") {
        di as error "Invalid VCE type: `vce'"
        di as error "Valid options: ols, robust/hc1, hac, cluster/cl, hc2, hc3, cr0, cr1"
        exit 198
    }
    
    // Aliases are mapped to canonical names for unified handling
    // Note: hc1 maps to robust (custom), hc3 is now distinct from hac (Arellano)
    if "`vce'" == "hc1" local vce "robust"
    if "`vce'" == "cl"  local vce "cluster"
    
    // Required input matrices are verified to exist before computation
    capture confirm matrix `x_ddm'
    if _rc {
        di as error "Matrix `x_ddm' not found"
        exit 111
    }
    
    capture confirm matrix `y_ddm'
    if _rc {
        di as error "Matrix `y_ddm' not found"
        exit 111
    }
    
    capture confirm matrix `beta'
    if _rc {
        di as error "Matrix `beta' not found"
        exit 111
    }
    
    capture confirm matrix `id'
    if _rc {
        di as error "Matrix `id' not found"
        exit 111
    }
    
    // Panel dimensions are inferred from the data if not explicitly specified
    if `n' == 0 {
        tempname id_mata
        mata: st_matrix("`id_mata'", uniqrows(st_matrix("`id'")))
        local n = rowsof(`id_mata')
    }
    
    if `n_t' == 0 {
        // Placeholder value; calling function provides actual T+1 period count
        local n_t = 1
    }
    
    // Clustering is performed at the individual level by default
    if "`cluster'" == "" {
        local cluster "`id'"
    }
    
    // -------------------------------------------------------------------------
    // VCE Computation
    // -------------------------------------------------------------------------
    
    // The appropriate VCE estimator is selected based on the specified type
    tempname vcov_mat se_vec
    
    if "`vce'" == "ols" {
        // Homoskedastic VCE: Sigma = sigma^2 * (X'X)^{-1}
        mata: st_matrix("`vcov_mat'", _eqt_vcov_ols( ///
            st_matrix("`x_ddm'"), ///
            st_matrix("`y_ddm'")[., 1], ///
            `n', `n_t'))
    }
    else if "`vce'" == "robust" {
        // HC1 heteroskedasticity-consistent VCE
        mata: st_matrix("`vcov_mat'", _eqt_vcov_hc1( ///
            st_matrix("`x_ddm'"), ///
            st_matrix("`y_ddm'")[., 1], ///
            st_matrix("`beta'")[., 1]))
    }
    else if "`vce'" == "hac" {
        // HC3 panel-robust VCE accounting for within-individual correlation
        mata: st_matrix("`vcov_mat'", _eqt_vcov_hc3( ///
            st_matrix("`x_ddm'"), ///
            st_matrix("`y_ddm'")[., 1], ///
            st_matrix("`beta'")[., 1], ///
            st_matrix("`id'")[., 1]))
    }
    else if "`vce'" == "cluster" {
        // CR0 cluster-robust VCE without finite-sample adjustment
        mata: st_matrix("`vcov_mat'", _eqt_vcov_cr0( ///
            st_matrix("`x_ddm'"), ///
            st_matrix("`y_ddm'")[., 1], ///
            st_matrix("`beta'")[., 1], ///
            st_matrix("`cluster'")[., 1]))
    }
    
    // Standard errors are computed from diagonal of VCE matrix
    mata: st_matrix("`se_vec'", _eqt_compute_se(st_matrix("`vcov_mat'")))
    
    // Placebo-only submatrix is extracted when testing hypotheses (3.1)-(3.3)
    if `no_placebos' > 0 {
        tempname vcov_placebo
        mata: st_matrix("`vcov_placebo'", ///
            _eqt_extract_vcov_placebo(st_matrix("`vcov_mat'"), `no_placebos'))
        return matrix vcov_placebo = `vcov_placebo'
    }
    
    // -------------------------------------------------------------------------
    // Results Storage and Return
    // -------------------------------------------------------------------------
    
    // Results are optionally stored in caller-specified matrices
    if "`vcov_result'" != "" {
        matrix `vcov_result' = `vcov_mat'
    }
    
    if "`se_result'" != "" {
        matrix `se_result' = `se_vec'
    }
    
    // Matrices are returned via r() for flexible downstream usage
    return matrix vcov = `vcov_mat'
    return matrix se = `se_vec'
    return local vce_type "`vce'"
    
end
