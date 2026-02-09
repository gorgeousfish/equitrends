*! _equitrends_vce_native.ado - Native Stata VCE computation using regress command
*!
*! This module provides variance-covariance matrix estimation using Stata's
*! built-in vce() options with the regress command. It operates on double-demeaned
*! data (after removing individual and time fixed effects) and uses noconstant
*! since the constant is absorbed by the fixed effects transformation.
*!
*! Syntax:
*!   _equitrends_vce_native varlist, VCE(string) [CLuster(varname)]
*!
*! Arguments:
*!   varlist   - First variable is Y_ddm (double-demeaned outcome), rest are X_ddm
*!   vce       - VCE type: hc2, hc3, cr0, cr1
*!   cluster   - Cluster variable (required for cr0, cr1)
*!
*! VCE Types:
*!   hc2   - HC2 heteroskedasticity-robust with leverage adjustment
*!           V_HC2 = (X'X)^{-1} * sum[e_i^2/(1-h_ii) * x_i*x_i'] * (X'X)^{-1}
*!   hc3   - HC3 more conservative leverage adjustment (Davidson & MacKinnon 1993)
*!           V_HC3 = (X'X)^{-1} * sum[e_i^2/(1-h_ii)^2 * x_i*x_i'] * (X'X)^{-1}
*!   cr0   - CR0 cluster-robust without small-sample adjustment
*!           V_CR0 = V_CR1 * (G-1)/G
*!   cr1   - CR1 cluster-robust with DF adjustment (Stata default)
*!           V_CR1 = [G/(G-1)] * (X'X)^{-1} * sum_g[X_g'e_g e_g'X_g] * (X'X)^{-1}
*!
*! Returns:
*!   r(b)   - Coefficient vector (1 x k)
*!   r(V)   - Variance-covariance matrix (k x k)
*!   r(N)   - Number of observations
*!   r(G)   - Number of clusters (for cr0/cr1 only)
*!
*! References:
*!   MacKinnon, J.G., and White, H. (1985). Some heteroskedasticity-consistent
*!     covariance matrix estimators with improved finite sample properties.
*!     Journal of Econometrics, 29(3):305-325.
*!   Cameron, A.C., and Miller, D.L. (2015). A practitioner's guide to
*!     cluster-robust inference. Journal of Human Resources, 50(2):317-372.

program define _equitrends_vce_native, rclass
    version 16.0
    
    // -------------------------------------------------------------------------
    // Syntax parsing
    // -------------------------------------------------------------------------
    syntax varlist(min=2 numeric), VCE(string) [CLuster(varname numeric)]
    
    // Parse dependent and independent variables
    gettoken y_ddm x_ddm : varlist
    
    // Normalize VCE type to lowercase
    local vce = lower("`vce'")
    
    // -------------------------------------------------------------------------
    // Input validation
    // -------------------------------------------------------------------------
    
    // Validate VCE type
    if !inlist("`vce'", "hc2", "hc3", "cr0", "cr1") {
        di as error "Invalid vce(`vce') for native VCE computation"
        di as error "Valid options: hc2, hc3, cr0, cr1"
        exit 198
    }
    
    // Cluster variable is required for cluster-robust VCE
    if inlist("`vce'", "cr0", "cr1") & "`cluster'" == "" {
        di as error "vce(`vce') requires cluster(varname) option"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // VCE computation using native Stata regression
    // -------------------------------------------------------------------------
    
    if "`vce'" == "hc2" {
        // HC2: Leverage-adjusted heteroskedasticity-robust
        // Uses Stata's native vce(hc2) option
        quietly regress `y_ddm' `x_ddm', vce(hc2) noconstant
        
        // Store results
        tempname b_mat V_mat
        matrix `b_mat' = e(b)
        matrix `V_mat' = e(V)
        local N_obs = e(N)
        
        return matrix b = `b_mat'
        return matrix V = `V_mat'
        return scalar N = `N_obs'
    }
    else if "`vce'" == "hc3" {
        // HC3: More conservative leverage adjustment
        // Uses Stata's native vce(hc3) option
        quietly regress `y_ddm' `x_ddm', vce(hc3) noconstant
        
        // Store results
        tempname b_mat V_mat
        matrix `b_mat' = e(b)
        matrix `V_mat' = e(V)
        local N_obs = e(N)
        
        return matrix b = `b_mat'
        return matrix V = `V_mat'
        return scalar N = `N_obs'
    }
    else if "`vce'" == "cr1" {
        // CR1: Cluster-robust with small-sample adjustment (Stata default)
        // Uses Stata's native vce(cluster) option
        quietly regress `y_ddm' `x_ddm', vce(cluster `cluster') noconstant
        
        // Store results
        tempname b_mat V_mat
        matrix `b_mat' = e(b)
        matrix `V_mat' = e(V)
        local N_obs = e(N)
        local N_clust = e(N_clust)
        
        return matrix b = `b_mat'
        return matrix V = `V_mat'
        return scalar N = `N_obs'
        return scalar G = `N_clust'
    }
    else if "`vce'" == "cr0" {
        // CR0: Cluster-robust WITHOUT small-sample adjustment
        // Stata's vce(cluster) uses CR1 with G/(G-1) adjustment
        // To get CR0, we multiply V_CR1 by (G-1)/G
        
        quietly regress `y_ddm' `x_ddm', vce(cluster `cluster') noconstant
        
        // Get number of clusters
        local G = e(N_clust)
        
        // Boundary check: need at least 2 clusters
        if `G' < 2 {
            di as error "CR0 VCE requires at least 2 clusters (got G=`G')"
            exit 198
        }
        
        // Apply CR0 adjustment: V_CR0 = V_CR1 * (G-1)/G
        // This undoes the small-sample correction that Stata applies
        tempname b_mat V_cr1 V_cr0
        matrix `b_mat' = e(b)
        matrix `V_cr1' = e(V)
        
        // Calculate adjustment factor
        local adj_factor = (`G' - 1) / `G'
        matrix `V_cr0' = `V_cr1' * `adj_factor'
        
        local N_obs = e(N)
        
        return matrix b = `b_mat'
        return matrix V = `V_cr0'
        return scalar N = `N_obs'
        return scalar G = `G'
    }
    
end
