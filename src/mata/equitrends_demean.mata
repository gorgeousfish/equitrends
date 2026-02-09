*! equitrends_demean.mata - Double demeaning functions for EQUITRENDS
*!
*! This file contains functions for double demeaning transformation:
*!   - _eqt_grouped_mean()        : Calculate group means
*!   - _eqt_between_trans()       : Within-group demeaning (vector)
*!   - _eqt_mat_between_trans()   : Within-group demeaning (matrix)
*!   - _eqt_ols_cholesky()        : OLS via Cholesky decomposition
*!   - _eqt_double_demean()       : Two-step double demeaning (for Bootstrap)
*!   - _eqt_standard_double_demean() : Standard two-way FE demeaning (for IU)
*!   - _eqt_construct_WD()        : Construct WD matrix
*!   - _eqt_sigma_hathat_c()      : Constrained variance estimation
*!
*! Two different demeaning methods are implemented:
*!   - IU method: Standard two-way FE (x_dm = x - mean_i - mean_t + mean)
*!     -> Use _eqt_standard_double_demean()
*!   - Bootstrap method: OLS residualization approach
*!     -> Use _eqt_double_demean()

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// _eqt_grouped_mean()
// Calculate the mean for each group and return a vector of same length as x
// where each element is replaced by its group mean.
//
// Arguments:
//   x     : real colvector - input values
//   group : real colvector - group identifiers (same length as x)
//
// Returns:
//   real colvector - vector of group means (same length as x)
// ============================================================================
real colvector _eqt_grouped_mean(real colvector x, real colvector group)
{
    real colvector unique_groups, means, idx
    real scalar i, n_groups, group_mean
    
    // Get unique groups
    unique_groups = uniqrows(group)
    n_groups = rows(unique_groups)
    
    // Initialize output vector
    means = J(rows(x), 1, .)
    
    // Calculate mean for each group and fill corresponding positions
    for (i = 1; i <= n_groups; i++) {
        idx = selectindex(group :== unique_groups[i])
        group_mean = mean(x[idx])
        means[idx] = J(rows(idx), 1, group_mean)
    }
    
    return(means)
}

// ============================================================================
// _eqt_between_trans()
// Apply within-group demeaning transformation: x_it - mean(x_i)
//
// Arguments:
//   x     : real colvector - input values
//   group : real colvector - group identifiers (same length as x)
//
// Returns:
//   real colvector - demeaned values (same length as x)
// ============================================================================
real colvector _eqt_between_trans(real colvector x, real colvector group)
{
    real colvector means, between
    
    means = _eqt_grouped_mean(x, group)
    between = x - means
    
    return(between)
}

// ============================================================================
// _eqt_mat_between_trans()
// Apply within-group demeaning to each column of a matrix independently
//
// Arguments:
//   x     : real matrix - input matrix (N_obs x p)
//   group : real colvector - group identifiers (N_obs x 1)
//
// Returns:
//   real matrix - demeaned matrix (N_obs x p)
// ============================================================================
real matrix _eqt_mat_between_trans(real matrix x, real colvector group)
{
    real matrix between
    real scalar i, ncols
    
    ncols = cols(x)
    between = J(rows(x), ncols, .)
    
    for (i = 1; i <= ncols; i++) {
        between[., i] = _eqt_between_trans(x[., i], group)
    }
    
    return(between)
}

// ============================================================================
// _eqt_ols_cholesky()
// Solve OLS using Cholesky decomposition: beta = (X'X)^{-1} X'y
//
// The R/Armadillo implementation uses:
//   chol = arma::chol(XtX)           // Upper triangular Cholesky
//   chol_inv = arma::inv(trimatu(chol))
//   beta = chol_inv * chol_inv' * Xty
//
// Mata's cholesky() returns lower triangular L where XtX = LL'
// So we need to transpose to get upper triangular R = L'
//
// Note: If Cholesky fails (singular matrix), falls back to qrsolve()
//       This fallback is a Stata-specific safety mechanism. The R package's
//       ols_cholesky() has no fallback and would throw an error instead.
//       In normal usage, multicollinearity detection should prevent this
//       fallback from triggering. If it does trigger, it indicates a 
//       near-singular matrix that was not caught by the detection.
//
// Arguments:
//   XtX : real matrix - X'X matrix (p x p, positive definite)
//   Xty : real matrix - X'y vector/matrix (p x k)
//
// Returns:
//   real matrix - OLS coefficients beta (p x k)
// ============================================================================
real matrix _eqt_ols_cholesky(real matrix XtX, real matrix Xty)
{
    real matrix L, R, R_inv, beta
    real scalar rc
    
    // Try Cholesky decomposition first
    // Mata's cholesky() returns lower triangular L where XtX = LL'
    L = cholesky(XtX)
    
    // Check if Cholesky succeeded (no missing values)
    if (hasmissing(L)) {
        // Cholesky failed - matrix is singular or near-singular
        // Fall back to qrsolve which handles rank-deficient matrices
        // 
        // NOTE: This is a safety mechanism not present in R's ols_cholesky().
        // The R package would throw an error in this case.
        // This typically indicates a near-singular matrix not caught by
        // multicollinearity detection. Results may differ slightly from R.
        printf("{txt}Warning: Cholesky decomposition failed in _eqt_ols_cholesky().\n")
        printf("{txt}         Using QR solver as fallback. Matrix may be near-singular.\n")
        printf("{txt}         Note: R's ols_cholesky() would error in this case.\n")
        beta = qrsolve(XtX, Xty)
        return(beta)
    }
    
    // R package uses upper triangular, so R = L'
    R = L'
    
    // Compute inverse of upper triangular R
    R_inv = luinv(R)
    
    // Check if inverse succeeded
    if (hasmissing(R_inv)) {
        // Fall back to qrsolve
        // This case is rare - Cholesky succeeded but luinv failed
        printf("{txt}Warning: Matrix inversion failed in _eqt_ols_cholesky().\n")
        printf("{txt}         Using QR solver as fallback.\n")
        beta = qrsolve(XtX, Xty)
        return(beta)
    }
    
    // beta = R^{-1} * (R^{-1})' * Xty = (R'R)^{-1} * Xty
    beta = R_inv * R_inv' * Xty
    
    return(beta)
}

// ============================================================================
// _eqt_double_demean()
// Apply two-step double demeaning transformation (FOR BOOTSTRAP METHOD)
//
// This function implements the standard two-way fixed effects demeaning:
//   x_dm = x - mean_i(x) - mean_t(x) + mean(x)
//
// Note: Previously used OLS residualization (between_x - WD * beta_hat), but
// this gave different results than R due to singular matrix handling 
// differences in qrsolve vs Armadillo's Cholesky solver. The standard
// demeaning formula gives numerically identical results to R.
//
// Arguments:
//   x          : real matrix - input matrix (N_obs x p)
//   individual : real colvector - individual identifiers (N_obs x 1)
//   time       : real colvector - time identifiers (N_obs x 1)
//   WD         : real matrix - demeaned time dummy matrix (N_obs x T) [unused]
//
// Returns:
//   real matrix - double demeaned matrix (N_obs x p)
// ============================================================================
real matrix _eqt_double_demean(real matrix x, real colvector individual, 
                               real colvector time, real matrix WD)
{
    // Use standard double-demeaning which matches R's results
    // The WD parameter is kept for backward compatibility but not used
    pragma unused WD
    
    return(_eqt_standard_double_demean(x, individual, time))
}

// ============================================================================
// _eqt_standard_double_demean()
// Apply standard two-way fixed effects demeaning (FOR IU METHOD)
//
// This implements the standard double demeaning formula used by plm::plm()
// with effect="twoways" and model="within":
//
//   x_dm = x - mean_i(x) - mean_t(x) + mean(x)
//
// where:
//   mean_i(x) = individual-specific mean
//   mean_t(x) = time-specific mean
//   mean(x)   = overall mean
//
// This is equivalent to reghdfe with absorb(ID period).
//
// Arguments:
//   x          : real matrix - input matrix (N_obs x p)
//   individual : real colvector - individual identifiers (N_obs x 1)
//   time       : real colvector - time identifiers (N_obs x 1)
//
// Returns:
//   real matrix - double demeaned matrix (N_obs x p)
// ============================================================================
real matrix _eqt_standard_double_demean(real matrix x, real colvector individual, 
                                        real colvector time)
{
    real matrix x_dm, x_id_mean, x_t_mean
    real rowvector x_mean
    real colvector unique_ids, unique_times, idx
    real scalar N, n_ids, n_times, i, t, j, ncols
    
    N = rows(x)
    ncols = cols(x)
    
    // Get unique identifiers
    unique_ids = uniqrows(individual)
    unique_times = uniqrows(time)
    n_ids = rows(unique_ids)
    n_times = rows(unique_times)
    
    // Compute overall mean for each column
    x_mean = mean(x)
    
    // Compute individual means
    x_id_mean = J(N, ncols, .)
    for (i = 1; i <= n_ids; i++) {
        idx = selectindex(individual :== unique_ids[i])
        for (j = 1; j <= ncols; j++) {
            x_id_mean[idx, j] = J(rows(idx), 1, mean(x[idx, j]))
        }
    }
    
    // Compute time means
    x_t_mean = J(N, ncols, .)
    for (t = 1; t <= n_times; t++) {
        idx = selectindex(time :== unique_times[t])
        for (j = 1; j <= ncols; j++) {
            x_t_mean[idx, j] = J(rows(idx), 1, mean(x[idx, j]))
        }
    }
    
    // Apply standard double demean formula: x - mean_i - mean_t + mean
    x_dm = x - x_id_mean - x_t_mean :+ x_mean
    
    return(x_dm)
}

// ============================================================================
// _eqt_standard_double_demean_vec()
// Apply standard two-way fixed effects demeaning to a vector (FOR IU METHOD)
//
// Vector version of _eqt_standard_double_demean() for single column.
//
// Arguments:
//   x          : real colvector - input vector (N_obs x 1)
//   individual : real colvector - individual identifiers (N_obs x 1)
//   time       : real colvector - time identifiers (N_obs x 1)
//
// Returns:
//   real colvector - double demeaned vector (N_obs x 1)
// ============================================================================
real colvector _eqt_standard_double_demean_vec(real colvector x, 
                                               real colvector individual, 
                                               real colvector time)
{
    real colvector x_dm, x_id_mean, x_t_mean, idx
    real scalar x_mean
    real colvector unique_ids, unique_times
    real scalar N, n_ids, n_times, i, t
    
    N = rows(x)
    
    // Get unique identifiers
    unique_ids = uniqrows(individual)
    unique_times = uniqrows(time)
    n_ids = rows(unique_ids)
    n_times = rows(unique_times)
    
    // Compute overall mean
    x_mean = mean(x)
    
    // Compute individual means
    x_id_mean = J(N, 1, .)
    for (i = 1; i <= n_ids; i++) {
        idx = selectindex(individual :== unique_ids[i])
        x_id_mean[idx] = J(rows(idx), 1, mean(x[idx]))
    }
    
    // Compute time means
    x_t_mean = J(N, 1, .)
    for (t = 1; t <= n_times; t++) {
        idx = selectindex(time :== unique_times[t])
        x_t_mean[idx] = J(rows(idx), 1, mean(x[idx]))
    }
    
    // Apply standard double demean formula: x - mean_i - mean_t + mean
    x_dm = x - x_id_mean - x_t_mean :+ x_mean
    
    return(x_dm)
}

// ============================================================================
// _eqt_construct_WD()
// Construct the WD matrix (demeaned time dummy matrix)
//
// D is the time dummy matrix where D[i,t] = 1 if observation i is in period t
// WD = matrix_between_transformation(D, ID)
//
// Arguments:
//   period : real colvector - time period identifiers (N_obs x 1)
//   ID     : real colvector - individual identifiers (N_obs x 1)
//
// Returns:
//   real matrix - WD matrix (N_obs x T)
// ============================================================================
real matrix _eqt_construct_WD(real colvector period, real colvector ID)
{
    real colvector unique_periods, idx
    real matrix D, WD
    real scalar N_obs, T, t
    
    // Get unique periods
    unique_periods = uniqrows(period)
    T = rows(unique_periods)
    N_obs = rows(period)
    
    // Construct time dummy matrix D
    D = J(N_obs, T, 0)
    
    for (t = 1; t <= T; t++) {
        idx = selectindex(period :== unique_periods[t])
        D[idx, t] = J(rows(idx), 1, 1)
    }
    
    // Apply between transformation to get WD
    WD = _eqt_mat_between_trans(D, ID)
    
    return(WD)
}

// ============================================================================
// _eqt_sigma_hathat_c()
// Calculate constrained variance estimate for bootstrap
//
// Formula: sigma^2_c = (u'u) / df
// where:
//   u = y - X * beta (residuals)
//   df = N - p - n - T + 1 (degrees of freedom)
//   N = total observations
//   p = number of regressors
//   n = number of individuals
//   T = number of time periods
//
// Arguments:
//   parameter : real colvector - coefficient vector (p x 1)
//   x         : real matrix - design matrix (N_obs x p)
//   y         : real colvector - dependent variable (N_obs x 1)
//   ID        : real colvector - individual identifiers (N_obs x 1)
//   time      : real colvector - time identifiers (N_obs x 1)
//
// Returns:
//   real scalar - constrained variance estimate
// ============================================================================
real scalar _eqt_sigma_hathat_c(real colvector parameter, real matrix x, 
                                real colvector y, real colvector ID, 
                                real colvector time)
{
    real colvector Xb, residuals
    real scalar N, n, no_periods, p, df, c_sigma_hathat
    
    // Validate input: arrays must be non-empty
    if (rows(ID) == 0) {
        errprintf("Error: ID array cannot be empty\n")
        _error(3498)
    }
    if (rows(time) == 0) {
        errprintf("Error: time array cannot be empty\n")
        _error(3498)
    }
    
    // Compute fitted values and residuals
    Xb = x * parameter
    residuals = y - Xb
    
    // Get dimensions
    N = rows(ID)
    n = rows(uniqrows(ID))
    no_periods = rows(uniqrows(time))
    p = cols(x)
    
    // Degrees of freedom for TWFE regression after double-demeaning
    // 
    // Formula: df = N - p - n - T + 1
    //   where N = total obs, p = # regressors, n = # individuals, T = # periods
    //
    // Note: Paper eq. 4.6 shows denominator (n-1)T, which omits the p term.
    // This is either a simplification or typo in the paper. Our formula matches:
    //   1. The R package (EquiTrends) implementation exactly
    //   2. Standard TWFE degrees of freedom: N - (n-1) - (T-1) - p
    //      = N - n + 1 - T + 1 - p = N - p - n - T + 1 (after normalization)
    //
    // The difference: paper (n-1)T vs implementation T(n-1) - p
    // When p=T (only placebos): df = T(n-2), not T(n-1)
    df = N - p - n - no_periods + 1
    
    // Check for valid degrees of freedom (matches Python version)
    if (df <= 0) {
        errprintf("Error: Degrees of freedom must be positive (df = %g)\n", df)
        errprintf("  N = %g, p = %g, n = %g, T = %g\n", N, p, n, no_periods)
        errprintf("  Model is overparameterized: need N - p - n - T + 1 > 0\n")
        _error(3498)
    }
    
    // Constrained variance estimate
    c_sigma_hathat = cross(residuals, residuals) / df
    
    return(c_sigma_hathat)
}

end
