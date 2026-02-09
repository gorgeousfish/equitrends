*! equitrends_vcov.mata - Variance-covariance matrix estimation for TWFE models
*!
*! Implements robust variance-covariance estimators for placebo coefficient
*! inference in the two-way fixed effects (TWFE) regression framework.
*!
*! Core VCE Functions:
*!   - _eqt_vcov_ols()           : Default OLS VCE (homoskedastic)
*!   - _eqt_vcov_hc1()           : HC1 heteroskedasticity-robust VCE
*!   - _eqt_vcov_hc1_cluster()   : HC1 cluster-robust VCE with small sample adjustment
*!   - _eqt_vcov_hc3()           : HC3 leverage-adjusted VCE (Arellano-type)
*!   - _eqt_vcov_cr0()           : CR0 cluster-robust VCE without adjustment
*!
*! Helper Functions:
*!   - _eqt_hat_diagonal()       : Hat matrix diagonal elements for HC3
*!   - _eqt_extract_vcov_placebo(): Extract placebo coefficient submatrix
*!   - _eqt_compute_se()         : Compute standard errors from VCE
*!
*! Legacy Functions (for backward compatibility):
*!   - _equitrends_extract_submatrix() : Extract upper-left submatrix
*!   - _equitrends_cluster_indices()   : Get cluster row indices
*!   - _equitrends_meat_matrix()       : Compute meat matrix
*!
*! Note: CR0 does NOT apply the (G/(G-1)) small sample adjustment

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// CORE VCE FUNCTIONS
// ============================================================================

// ============================================================================
// _eqt_vcov_ols()
// Compute default OLS variance-covariance matrix
//
// Formula: V = sigma^2 * (X'X)^{-1}
// where sigma^2 = sum(u^2) / df
//       df = N - p - n - T + 1
//
// Arguments:
//   X_ddm : real matrix - double-demeaned design matrix (N x p)
//   Y_ddm : real colvector - double-demeaned dependent variable (N x 1)
//   n     : real scalar - number of individuals
//   n_t   : real scalar - number of time periods
//
// Returns:
//   real matrix - OLS VCE matrix (p x p)
// ============================================================================
real matrix _eqt_vcov_ols(real matrix X_ddm, real colvector Y_ddm, 
                          real scalar n, real scalar n_t)
{
    real matrix XtX, XtX_inv, vcov
    real colvector beta, residuals
    real scalar N, p, df, sigma_sq
    
    // Dimensions
    N = rows(X_ddm)
    p = cols(X_ddm)
    
    // Compute (X'X)^{-1}
    XtX = cross(X_ddm, X_ddm)
    XtX_inv = invsym(XtX)
    
    // OLS coefficients
    beta = XtX_inv * cross(X_ddm, Y_ddm)
    
    // Residuals
    residuals = Y_ddm - X_ddm * beta
    
    // Degrees of freedom: df = N - p - n - T + 1
    // Accounts for individual fixed effects (n) and time fixed effects (T)
    df = N - p - n - n_t + 1
    
    // Check for valid degrees of freedom
    if (df <= 0) {
        errprintf("Error: Degrees of freedom is non-positive (%g)\n", df)
        return(J(p, p, .))
    }
    
    // Variance estimate
    sigma_sq = cross(residuals, residuals) / df
    
    // OLS VCE
    vcov = sigma_sq * XtX_inv
    
    return(vcov)
}

// ============================================================================
// _eqt_hat_diagonal()
// Compute diagonal elements of the hat matrix H = X(X'X)^{-1}X'
//
// Reference: Used for HC3 VCE computation
//
// Formula: h_i = x_i' (X'X)^{-1} x_i
// Optimized: H_diag = rowsum(X * (X'X)^{-1} :* X)
//
// Arguments:
//   X_ddm : real matrix - double-demeaned design matrix (N x p)
//
// Returns:
//   real colvector - hat matrix diagonal elements (N x 1)
//                    all elements are in [0, 1]
// ============================================================================
real colvector _eqt_hat_diagonal(real matrix X_ddm)
{
    real matrix XtX_inv, X_XtXinv
    real colvector h_diag
    
    // Compute (X'X)^{-1}
    XtX_inv = invsym(cross(X_ddm, X_ddm))
    
    // Optimized computation: h_i = sum_j (X[i,j] * (X * XtX_inv)[i,j])
    X_XtXinv = X_ddm * XtX_inv
    h_diag = rowsum(X_XtXinv :* X_ddm)
    
    return(h_diag)
}

// ============================================================================
// _eqt_vcov_hc1()
// Compute HC1 (White) heteroskedasticity-robust VCE
//
// Formula: V = (N/(N-K)) * (X'X)^{-1} * M * (X'X)^{-1}
// where M = sum_i (u_i^2 * x_i * x_i') is the meat matrix
//
// Arguments:
//   X_ddm : real matrix - double-demeaned design matrix (N x p)
//   Y_ddm : real colvector - double-demeaned dependent variable (N x 1)
//   beta  : real colvector - OLS coefficient estimates (p x 1)
//
// Returns:
//   real matrix - HC1 VCE matrix (p x p)
// ============================================================================
real matrix _eqt_vcov_hc1(real matrix X_ddm, real colvector Y_ddm, 
                          real colvector beta)
{
    real matrix XtX_inv, meat, vcov
    real colvector residuals, u_sq
    real scalar N, K, adjustment
    
    // Dimensions
    N = rows(X_ddm)
    K = cols(X_ddm)
    
    // BOUNDARY CHECK: N must be greater than K for HC1 small sample adjustment
    if (N <= K) {
        errprintf("Error: HC1 VCE requires N > K (got N=%g, K=%g)\n", N, K)
        return(J(K, K, .))
    }
    
    // Compute (X'X)^{-1}
    XtX_inv = invsym(cross(X_ddm, X_ddm))
    
    // Residuals
    residuals = Y_ddm - X_ddm * beta
    
    // Squared residuals
    u_sq = residuals :^ 2
    
    // HC1 meat matrix: sum_i (u_i^2 * x_i * x_i')
    // Efficient computation: X' * diag(u^2) * X = (X :* u)' * (X :* u)
    // But we need X' * diag(u^2) * X which is cross(X, u_sq, X)
    meat = cross(X_ddm, u_sq, X_ddm)
    
    // HC1 small sample adjustment: N / (N - K)
    adjustment = N / (N - K)
    meat = meat * adjustment
    
    // Sandwich estimator
    vcov = XtX_inv * meat * XtX_inv
    
    return(vcov)
}

// ============================================================================
// _eqt_vcov_hc1_cluster()
// Compute HC1 cluster-robust VCE
//
// Formula: V = adj * (X'X)^{-1} * M * (X'X)^{-1}
// where M = sum_g (score_g * score_g')
//       score_g = sum_{i in g} (u_i * x_i)
//       adj = (G/(G-1)) * (N-1)/(N-K)
//
// CRITICAL: This is different from _eqt_vcov_hc1 (White HC1) which does NOT cluster
//           and from _eqt_vcov_cr0 which uses different adjustment
//
// Arguments:
//   X_ddm       : real matrix - double-demeaned design matrix (N x p)
//   Y_ddm       : real colvector - double-demeaned dependent variable (N x 1)
//   beta        : real colvector - OLS coefficient estimates (p x 1)
//   cluster_var : real colvector - cluster variable (N x 1)
//
// Returns:
//   real matrix - HC1 cluster VCE matrix (p x p)
// ============================================================================
real matrix _eqt_vcov_hc1_cluster(real matrix X_ddm, real colvector Y_ddm, 
                                   real colvector beta, real colvector cluster_var)
{
    real matrix XtX_inv, meat, vcov, X_g
    real colvector residuals, unique_clusters, idx, u_g, score_g
    real scalar G, K, N, g, adjustment
    
    // Dimensions
    N = rows(X_ddm)
    K = cols(X_ddm)
    
    // BOUNDARY CHECK: N must be greater than K
    if (N <= K) {
        errprintf("Error: HC1 cluster VCE requires N > K (got N=%g, K=%g)\n", N, K)
        return(J(K, K, .))
    }
    
    // Compute (X'X)^{-1}
    XtX_inv = invsym(cross(X_ddm, X_ddm))
    
    // Residuals
    residuals = Y_ddm - X_ddm * beta
    
    // Get unique clusters
    unique_clusters = uniqrows(cluster_var)
    G = rows(unique_clusters)
    
    // BOUNDARY CHECK: G must be at least 2 for cluster adjustment
    if (G < 2) {
        errprintf("Error: HC1 cluster VCE requires at least 2 clusters (got G=%g)\n", G)
        return(J(K, K, .))
    }
    
    // Initialize meat matrix
    meat = J(K, K, 0)
    
    // Compute cluster meat matrix: sum_g (score_g * score_g')
    // where score_g = X_g' * u_g = sum_{i in g} (u_i * x_i)
    for (g = 1; g <= G; g++) {
        // Get indices for this cluster
        idx = selectindex(cluster_var :== unique_clusters[g])
        
        // Extract cluster data
        X_g = X_ddm[idx, .]
        u_g = residuals[idx]
        
        // Cluster score: X_g' * u_g
        score_g = cross(X_g, u_g)
        
        // Add outer product to meat matrix
        meat = meat + score_g * score_g'
    }
    
    // HC1 cluster small sample adjustment: (G/(G-1)) * (N-1)/(N-K)
    adjustment = (G / (G - 1)) * ((N - 1) / (N - K))
    meat = meat * adjustment
    
    // Sandwich estimator
    vcov = XtX_inv * meat * XtX_inv
    
    return(vcov)
}

// ============================================================================
// _eqt_vcov_hc3()
// Compute HC3 (Arellano) panel HAC VCE
//
// Formula: V = (X'X)^{-1} * M * (X'X)^{-1}
// where M = sum_i (X_i' * Omega_i * X_i)
//       Omega_i = diag(u_it^2 / (1 - h_it)^2)
//
// Arguments:
//   X_ddm : real matrix - double-demeaned design matrix (N x p)
//   Y_ddm : real colvector - double-demeaned dependent variable (N x 1)
//   beta  : real colvector - OLS coefficient estimates (p x 1)
//   id    : real colvector - individual identifiers (N x 1)
//
// Returns:
//   real matrix - HC3 VCE matrix (p x p)
// ============================================================================
real matrix _eqt_vcov_hc3(real matrix X_ddm, real colvector Y_ddm, 
                          real colvector beta, real colvector id)
{
    real matrix XtX_inv, meat, vcov, X_i
    real colvector residuals, h_diag, unique_ids, idx, u_i, h_i, u_adj, score_i
    real scalar n, K, i
    
    // Dimensions
    K = cols(X_ddm)
    
    // Compute (X'X)^{-1}
    XtX_inv = invsym(cross(X_ddm, X_ddm))
    
    // Residuals
    residuals = Y_ddm - X_ddm * beta
    
    // Hat matrix diagonal elements
    h_diag = _eqt_hat_diagonal(X_ddm)
    
    // Get unique individuals
    unique_ids = uniqrows(id)
    n = rows(unique_ids)
    
    // Initialize meat matrix
    meat = J(K, K, 0)
    
    // Compute HC3 meat matrix by individual
    for (i = 1; i <= n; i++) {
        // Get indices for this individual
        idx = selectindex(id :== unique_ids[i])
        
        // Extract individual data
        X_i = X_ddm[idx, .]
        u_i = residuals[idx]
        h_i = h_diag[idx]
        
        // HC3 adjustment: u_adj = u / (1 - h)
        // STABILITY CHECK: cap h_i to avoid numerical instability when h_i approaches 1
        if (max(h_i) > 0.99) {
            printf("{txt}Warning: High-leverage observations detected (max h_i = %g). HC3 VCE may be unstable.\n", max(h_i))
        }
        // Cap denominator to prevent division by near-zero using element-wise maximum
        u_adj = u_i :/ rowmax((1 :- h_i, J(rows(h_i), 1, 1e-10)))
        
        // Individual score: X_i' * u_adj
        score_i = cross(X_i, u_adj)
        
        // Add outer product to meat matrix
        meat = meat + score_i * score_i'
    }
    
    // Sandwich estimator (no additional adjustment for HC3)
    vcov = XtX_inv * meat * XtX_inv
    
    return(vcov)
}

// ============================================================================
// _eqt_vcov_cr0()
// Compute CR0 cluster-robust VCE WITHOUT small sample adjustment
//
// CRITICAL: CR0 does NOT apply the G/(G-1) adjustment that Stata's
//           vce(cluster) uses by default (which is CR1)
//
// Formula: V = (X'X)^{-1} * M * (X'X)^{-1}
// where M = sum_g (X_g' * u_g * u_g' * X_g)
//
// CR0 vs CR1:
//   CR0: V = (X'X)^{-1} * M * (X'X)^{-1}           (no adjustment)
//   CR1: V = (G/(G-1)) * (X'X)^{-1} * M * (X'X)^{-1}  (Stata default)
//   Therefore: CR0 = CR1 * (G-1)/G
//
// Arguments:
//   X_ddm       : real matrix - double-demeaned design matrix (N x p)
//   Y_ddm       : real colvector - double-demeaned dependent variable (N x 1)
//   beta        : real colvector - OLS coefficient estimates (p x 1)
//   cluster_var : real colvector - cluster variable (N x 1)
//                 if empty/missing, defaults to individual ID
//
// Returns:
//   real matrix - CR0 VCE matrix (p x p)
// ============================================================================
real matrix _eqt_vcov_cr0(real matrix X_ddm, real colvector Y_ddm, 
                          real colvector beta, real colvector cluster_var)
{
    real matrix XtX_inv, meat, vcov, X_g
    real colvector residuals, unique_clusters, idx, u_g, score_g
    real scalar G, K, g
    
    // Dimensions
    K = cols(X_ddm)
    
    // Compute (X'X)^{-1}
    XtX_inv = invsym(cross(X_ddm, X_ddm))
    
    // Residuals
    residuals = Y_ddm - X_ddm * beta
    
    // Get unique clusters
    unique_clusters = uniqrows(cluster_var)
    G = rows(unique_clusters)
    
    // BOUNDARY CHECK: G must be at least 2 for cluster-robust variance estimation
    // With only 1 cluster, cluster-robust VCE is undefined (no between-cluster variation)
    if (G < 2) {
        errprintf("Error: CR0 VCE requires at least 2 clusters (got G=%g)\n", G)
        return(J(K, K, .))
    }
    
    // Initialize meat matrix
    meat = J(K, K, 0)
    
    // Compute CR0 meat matrix by cluster
    for (g = 1; g <= G; g++) {
        // Get indices for this cluster
        idx = selectindex(cluster_var :== unique_clusters[g])
        
        // Extract cluster data
        X_g = X_ddm[idx, .]
        u_g = residuals[idx]
        
        // Cluster score: X_g' * u_g
        score_g = cross(X_g, u_g)
        
        // Add outer product to meat matrix
        meat = meat + score_g * score_g'
    }
    
    // CRITICAL: CR0 does NOT apply small sample adjustment
    // Do NOT multiply by G/(G-1)
    
    // Sandwich estimator
    vcov = XtX_inv * meat * XtX_inv
    
    return(vcov)
}

// ============================================================================
// _eqt_extract_vcov_placebo()
// Extract the placebo coefficient submatrix from full VCE matrix
//
// Arguments:
//   vcov        : real matrix - full VCE matrix (p x p)
//   no_placebos : real scalar - number of placebo coefficients
//
// Returns:
//   real matrix - placebo VCE submatrix (no_placebos x no_placebos)
//                 returns empty matrix if no_placebos invalid
// ============================================================================
real matrix _eqt_extract_vcov_placebo(real matrix vcov, real scalar no_placebos)
{
    if (no_placebos <= 0 | no_placebos > rows(vcov) | no_placebos > cols(vcov)) {
        return(J(0, 0, .))
    }
    return(vcov[1..no_placebos, 1..no_placebos])
}

// ============================================================================
// _eqt_compute_se()
// Compute standard errors from VCE matrix
//
// Formula: SE = sqrt(diag(V))
//
// Arguments:
//   vcov : real matrix - VCE matrix (k x k)
//
// Returns:
//   real colvector - standard errors (k x 1)
// ============================================================================
real colvector _eqt_compute_se(real matrix vcov)
{
    if (rows(vcov) == 0 | cols(vcov) == 0) {
        return(J(0, 1, .))
    }
    return(sqrt(diagonal(vcov)))
}

// ============================================================================
// LEGACY FUNCTIONS (for backward compatibility)
// ============================================================================

// ============================================================================
// _equitrends_extract_submatrix()
// Extract upper-left submatrix of specified dimension
//
// Arguments:
//   M    - real matrix, the source matrix
//   dim  - real scalar, dimension of submatrix to extract
//
// Returns:
//   real matrix - upper-left (dim x dim) submatrix
//
// Notes:
//   - Used to extract vcov_placebo from full vcov matrix
//   - Returns empty matrix if dim > rows(M) or dim > cols(M)
// ============================================================================
real matrix _equitrends_extract_submatrix(real matrix M, real scalar dim)
{
    if (dim <= 0 | dim > rows(M) | dim > cols(M)) {
        return(J(0, 0, .))
    }
    return(M[1..dim, 1..dim])
}

// ============================================================================
// _equitrends_cluster_indices()
// Get row indices for each cluster (individual)
//
// Arguments:
//   id         - real colvector, individual identifiers
//   unique_ids - real colvector, unique individual identifiers
//
// Returns:
//   pointer vector - pointers to index vectors for each cluster
//
// Notes:
//   - Used for clustered variance-covariance computation
//   - Each pointer points to a colvector of row indices
// ============================================================================
pointer(real colvector) colvector _equitrends_cluster_indices(
    real colvector id, 
    real colvector unique_ids)
{
    pointer(real colvector) colvector indices
    real scalar n_clusters, i
    
    n_clusters = rows(unique_ids)
    indices = J(n_clusters, 1, NULL)
    
    for (i = 1; i <= n_clusters; i++) {
        indices[i] = &selectindex(id :== unique_ids[i])
    }
    
    return(indices)
}

// ============================================================================
// _equitrends_meat_matrix()
// Compute the "meat" of the sandwich estimator for clustered standard errors
//
// Arguments:
//   X_ddm      - real matrix, double-demeaned design matrix (N x p)
//   residuals  - real colvector, OLS residuals (N x 1)
//   id         - real colvector, individual identifiers (N x 1)
//   unique_ids - real colvector, unique individual identifiers (n x 1)
//
// Returns:
//   real matrix - meat matrix (p x p)
//
// Formula:
//   M = sum_i (X_i' * u_i) * (X_i' * u_i)'
//
// Notes:
//   - Does NOT include small sample adjustment (n/(n-1))
//   - Caller should apply adjustment after calling this function
// ============================================================================
real matrix _equitrends_meat_matrix(
    real matrix X_ddm,
    real colvector residuals,
    real colvector id,
    real colvector unique_ids)
{
    real matrix meat, X_i
    real colvector u_i, score_i, idx_i
    real scalar n_clusters, p, i
    
    n_clusters = rows(unique_ids)
    p = cols(X_ddm)
    meat = J(p, p, 0)
    
    for (i = 1; i <= n_clusters; i++) {
        // Get indices for this cluster
        idx_i = selectindex(id :== unique_ids[i])
        
        // Extract cluster data
        X_i = X_ddm[idx_i, .]
        u_i = residuals[idx_i]
        
        // Compute score for this cluster: X_i' * u_i (p x 1)
        score_i = cross(X_i, u_i)
        
        // Add outer product to meat matrix
        meat = meat + score_i * score_i'
    }
    
    return(meat)
}

// ============================================================================
// NATIVE STATA VCE FUNCTIONS
// These functions use Stata's built-in vce() options via the regress command
// ============================================================================

// ============================================================================
// _eqt_vcov_native()
// Compute VCE using Stata's native regression command
//
// This function creates temporary Stata variables from Mata matrices,
// calls the _equitrends_vce_native ado program, and retrieves the results.
//
// Arguments:
//   X_ddm       : real matrix - double-demeaned design matrix (N x p)
//   Y_ddm       : real colvector - double-demeaned dependent variable (N x 1)
//   vce_type    : string scalar - VCE type: "hc2", "hc3", "cr0", "cr1"
//   cluster_var : real colvector - cluster variable (N x 1), required for cr0/cr1
//
// Returns:
//   real matrix - VCE matrix (p x p)
//
// Notes:
//   - Creates temporary variables __eqt_y_ddm, __eqt_x_ddm_*, __eqt_cluster
//   - These are dropped after VCE computation
//   - Requires _equitrends_vce_native.ado to be available
// ============================================================================
real matrix _eqt_vcov_native(real matrix X_ddm, real colvector Y_ddm, 
                              string scalar vce_type, real colvector cluster_var)
{
    real scalar N, p, i, rc
    real matrix vcov
    string scalar varlist, cmd, x_varname
    string rowvector x_varnames
    
    N = rows(X_ddm)
    p = cols(X_ddm)
    
    // Validate inputs
    if (N == 0 | p == 0) {
        errprintf("Error: Empty data matrices provided to _eqt_vcov_native\n")
        return(J(p, p, .))
    }
    
    if (rows(Y_ddm) != N) {
        errprintf("Error: Y_ddm and X_ddm have different row counts\n")
        return(J(p, p, .))
    }
    
    // Preserve current data and create temporary dataset
    stata("preserve")
    stata("quietly clear")
    stata(sprintf("quietly set obs %g", N))
    
    // Create temporary Y variable
    (void) st_addvar("double", "__eqt_y_ddm")
    st_store(., "__eqt_y_ddm", Y_ddm)
    
    // Create temporary X variables
    x_varnames = J(1, p, "")
    for (i = 1; i <= p; i++) {
        x_varname = sprintf("__eqt_x_ddm_%g", i)
        x_varnames[i] = x_varname
        (void) st_addvar("double", x_varname)
        st_store(., x_varname, X_ddm[., i])
    }
    
    // Create cluster variable if needed
    if (vce_type == "cr0" | vce_type == "cr1") {
        if (rows(cluster_var) != N) {
            errprintf("Error: cluster_var has wrong length for native VCE\n")
            stata("restore")
            return(J(p, p, .))
        }
        (void) st_addvar("double", "__eqt_cluster")
        st_store(., "__eqt_cluster", cluster_var)
    }
    
    // Build varlist for _equitrends_vce_native
    varlist = "__eqt_y_ddm"
    for (i = 1; i <= p; i++) {
        varlist = varlist + " " + x_varnames[i]
    }
    
    // Build and execute command
    if (vce_type == "cr0" | vce_type == "cr1") {
        cmd = sprintf("_equitrends_vce_native %s, vce(%s) cluster(__eqt_cluster)", 
                      varlist, vce_type)
    }
    else {
        cmd = sprintf("_equitrends_vce_native %s, vce(%s)", varlist, vce_type)
    }
    
    // Execute native VCE computation
    rc = _stata(cmd, 1)
    
    if (rc != 0) {
        errprintf("Error: _equitrends_vce_native failed with return code %g\n", rc)
        stata("restore")
        return(J(p, p, .))
    }
    
    // Retrieve VCE matrix from r(V)
    vcov = st_matrix("r(V)")
    
    // Restore original data
    stata("restore")
    
    // Validate result dimensions
    if (rows(vcov) != p | cols(vcov) != p) {
        errprintf("Error: Native VCE returned wrong dimensions (%g x %g, expected %g x %g)\n",
                  rows(vcov), cols(vcov), p, p)
        return(J(p, p, .))
    }
    
    return(vcov)
}

// ============================================================================
// _eqt_vcov_hc2_native()
// Convenience wrapper for HC2 native VCE
// ============================================================================
real matrix _eqt_vcov_hc2_native(real matrix X_ddm, real colvector Y_ddm)
{
    return(_eqt_vcov_native(X_ddm, Y_ddm, "hc2", J(0, 1, .)))
}

// ============================================================================
// _eqt_vcov_hc3_native()
// Convenience wrapper for HC3 native VCE
// ============================================================================
real matrix _eqt_vcov_hc3_native(real matrix X_ddm, real colvector Y_ddm)
{
    return(_eqt_vcov_native(X_ddm, Y_ddm, "hc3", J(0, 1, .)))
}

// ============================================================================
// _eqt_vcov_cr0_native()
// Convenience wrapper for CR0 native VCE
// ============================================================================
real matrix _eqt_vcov_cr0_native(real matrix X_ddm, real colvector Y_ddm, 
                                  real colvector cluster_var)
{
    return(_eqt_vcov_native(X_ddm, Y_ddm, "cr0", cluster_var))
}

// ============================================================================
// _eqt_vcov_cr1_native()
// Convenience wrapper for CR1 native VCE
// ============================================================================
real matrix _eqt_vcov_cr1_native(real matrix X_ddm, real colvector Y_ddm, 
                                  real colvector cluster_var)
{
    return(_eqt_vcov_native(X_ddm, Y_ddm, "cr1", cluster_var))
}

end
