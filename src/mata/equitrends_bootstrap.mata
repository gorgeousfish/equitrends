*! equitrends_bootstrap.mata - Bootstrap functions for EQUITRENDS
*!
*! This file contains bootstrap-related functions for the EQUITRENDS package:
*!   - _eqt_bootstrap_result        : Result structure for bootstrap tests
*!   - _eqt_rademacher()            : Generate Rademacher random variables
*!   - _eqt_id_to_index()           : Map IDs to consecutive indices
*!   - _eqt_quantile()              : Compute quantile of a vector
*!   - _eqt_boot_objective()        : MSE objective function for optimization
*!   - _eqt_boot_constraint()       : Constraint function for optimization
*!   - _eqt_boot_optimization()     : Constrained OLS optimization
*!   - _eqt_spherical_bootstrap()   : Spherical error bootstrap
*!   - _eqt_wild_bootstrap()        : Wild bootstrap with Rademacher weights
*!   - _eqt_bootstrap_critical_value() : Compute bootstrap critical value
*!   - _eqt_bootstrap_test()        : Main bootstrap test function

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// STRUCTURE DEFINITION: _eqt_bootstrap_result
// ============================================================================

// Result structure for bootstrap equivalence tests
struct _eqt_bootstrap_result {
    // Coefficient estimates
    real colvector placebo_coefs      // Placebo coefficients (unconstrained)
    real colvector constrained_coefs  // Constrained coefficients
    real scalar    max_abs_coef       // Maximum absolute placebo coefficient
    
    // Bootstrap results
    real colvector bootstrap_maxcoefs // Bootstrap max absolute coefficients (B x 1)
    real scalar    critical_value     // Critical value c*_alpha
    real scalar    reject_H0          // Whether to reject H0 (1=reject, 0=not reject)
    
    // Variance estimate
    real scalar    sigma_hathat_c     // Constrained variance estimate
    
    // Metadata
    real scalar    B                  // Number of bootstrap replications
    real scalar    alpha              // Significance level
    real scalar    delta              // Equivalence threshold
    string scalar  type               // "Boot" or "Wild"
    real scalar    n_individuals      // Number of individuals
    real scalar    n_periods          // Number of time periods
    real scalar    n_obs              // Number of observations
    
    // Diagnostic information
    real scalar    optimization_converged  // Whether optimization converged
    real scalar    optimization_iterations // Number of optimization iterations
}


// ============================================================================
// _eqt_matrix_exists()
// Helper function to check if a Stata matrix exists
// Must be defined before _eqt_boot_optimization_python which uses it
// ============================================================================
real scalar _eqt_matrix_exists(string scalar name)
{
    real scalar rc
    stata("capture confirm matrix " + name)
    stata("scalar __eqt_rc = _rc")
    rc = st_numscalar("__eqt_rc")
    return(rc == 0 ? 1 : 0)
}

// ============================================================================
// _eqt_scalar_exists()
// Helper function to check if a Stata scalar exists
// ============================================================================
real scalar _eqt_scalar_exists(string scalar name)
{
    real scalar rc
    stata("capture confirm scalar " + name)
    stata("scalar __eqt_rc = _rc")
    rc = st_numscalar("__eqt_rc")
    return(rc == 0 ? 1 : 0)
}


// ============================================================================
// _eqt_rademacher()
// Generate Rademacher random variables: P(R=1) = P(R=-1) = 0.5
//
// Arguments:
//   n : real scalar - number of random variables to generate
//
// Returns:
//   real colvector - vector of Rademacher variables (n x 1)
// ============================================================================
real colvector _eqt_rademacher(real scalar n)
{
    real colvector R
    
    // Generate Bernoulli(0.5) and transform to {-1, +1}
    // rbinomial(n, 1, 1, 0.5) gives 0 or 1 with equal probability
    // 2*x - 1 transforms to -1 or +1
    R = 2 * rbinomial(n, 1, 1, 0.5) :- 1
    
    return(R)
}

// ============================================================================
// _eqt_id_to_index()
// Create mapping from ID values to consecutive indices (1, 2, ..., n)
// This handles non-consecutive IDs (e.g., IDs like 5, 10, 15, ...)
//
// Arguments:
//   ID : real colvector - individual identifiers (N_obs x 1)
//
// Returns:
//   real colvector - index for each observation (N_obs x 1)
//                    where index[i] is the position of ID[i] in unique IDs
// ============================================================================
real colvector _eqt_id_to_index(real colvector ID)
{
    real colvector unique_IDs, indices
    real scalar N, n, i, j
    
    N = rows(ID)
    unique_IDs = uniqrows(ID)
    n = rows(unique_IDs)
    
    // Create index vector
    indices = J(N, 1, .)
    
    // Map each ID to its index in unique_IDs
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= n; j++) {
            if (ID[i] == unique_IDs[j]) {
                indices[i] = j
                break
            }
        }
    }
    
    return(indices)
}

// ============================================================================
// _eqt_quantile()
// Compute the alpha quantile of a vector
// Matches R's quantile() function with type=7 (default)
//
// Arguments:
//   x     : real colvector - input vector
//   alpha : real scalar    - quantile level (0 < alpha < 1)
//
// Returns:
//   real scalar - the alpha quantile of x
// ============================================================================
real scalar _eqt_quantile(real colvector x, real scalar alpha)
{
    real colvector sorted_x
    real scalar n, h, q_low, q_high, quantile_val
    
    n = rows(x)
    sorted_x = sort(x, 1)
    
    // R's type=7 quantile: linear interpolation
    // h = (n-1)*p + 1, then interpolate between floor(h) and ceil(h)
    h = (n - 1) * alpha + 1
    
    if (h <= 1) {
        quantile_val = sorted_x[1]
    }
    else if (h >= n) {
        quantile_val = sorted_x[n]
    }
    else {
        q_low = sorted_x[floor(h)]
        q_high = sorted_x[ceil(h)]
        quantile_val = q_low + (h - floor(h)) * (q_high - q_low)
    }
    
    return(quantile_val)
}


// ============================================================================
// _eqt_boot_objective()
// MSE objective function for constrained optimization
//
// Arguments:
//   beta        : real colvector - coefficient vector (p x 1)
//   X           : real matrix    - design matrix (N x p)
//   Y           : real colvector - dependent variable (N x 1)
//   no_placebos : real scalar    - number of placebo coefficients
//   delta       : real scalar    - equivalence threshold
//
// Returns:
//   real scalar - MSE = (1/N) * sum((Y - X*beta)^2)
// ============================================================================
real scalar _eqt_boot_objective(real colvector beta, real matrix X, 
                                real colvector Y, real scalar no_placebos,
                                real scalar delta)
{
    real colvector residuals
    real scalar N, MSE
    
    N = rows(Y)
    residuals = Y - X * beta
    MSE = cross(residuals, residuals) / N
    
    return(MSE)
}

// ============================================================================
// _eqt_boot_constraint()
// Constraint function for constrained optimization
// Constraint: max|beta[1:no_placebos]| >= delta
// Transformed to: -max|beta[1:no_placebos]| + delta <= 0
//
// Arguments:
//   beta        : real colvector - coefficient vector (p x 1)
//   X           : real matrix    - design matrix (N x p) [unused, for interface]
//   Y           : real colvector - dependent variable (N x 1) [unused]
//   no_placebos : real scalar    - number of placebo coefficients
//   delta       : real scalar    - equivalence threshold
//
// Returns:
//   real scalar - constraint value (-max|beta_placebo| + delta)
//                 constraint satisfied when return value <= 0
// ============================================================================
real scalar _eqt_boot_constraint(real colvector beta, real matrix X,
                                 real colvector Y, real scalar no_placebos,
                                 real scalar delta)
{
    real colvector placebo_coefs
    real scalar max_abs, constraint_val
    
    // Extract placebo coefficients (first no_placebos elements)
    placebo_coefs = beta[1::no_placebos]
    
    // Compute max absolute value
    max_abs = max(abs(placebo_coefs))
    
    // Constraint: -max|beta| + delta <= 0
    constraint_val = -max_abs + delta
    
    return(constraint_val)
}


// ============================================================================
// _eqt_boot_penalized_objective()
// Penalized MSE objective function for Mata optimize()
// Combines MSE with penalty for constraint violation
//
// Arguments (in order required by optimize()):
//   todo        : real scalar    - what to compute (0=value only)
//   beta        : real rowvector - coefficient vector (1 x p)
//   X           : real matrix    - design matrix (N x p) [extra arg 1]
//   Y           : real colvector - dependent variable (N x 1) [extra arg 2]
//   no_placebos : real scalar    - number of placebo coefficients [extra arg 3]
//   delta       : real scalar    - equivalence threshold [extra arg 4]
//   lambda      : real scalar    - penalty coefficient [extra arg 5]
//   v           : real scalar    - objective value (output)
//   g           : real rowvector - gradient (output, unused for d0)
//   H           : real matrix    - Hessian (output, unused for d0)
// ============================================================================
void _eqt_boot_penalized_objective(real scalar todo, real rowvector beta,
                                   real matrix X, real colvector Y,
                                   real scalar no_placebos, real scalar delta,
                                   real scalar lambda,
                                   real scalar v, real rowvector g, real matrix H)
{
    real colvector residuals, placebo_coefs
    real scalar N, MSE, max_abs, penalty
    
    N = rows(Y)
    
    // Compute MSE
    residuals = Y - X * beta'
    MSE = cross(residuals, residuals) / N
    
    // Compute constraint violation penalty
    // Constraint: max|beta[1:no_placebos]| >= delta
    // Penalty when violated: lambda * (delta - max|beta|)^2
    placebo_coefs = beta[1::no_placebos]'
    max_abs = max(abs(placebo_coefs))
    
    if (max_abs < delta) {
        penalty = lambda * (delta - max_abs)^2
    }
    else {
        penalty = 0
    }
    
    // Total objective (minimize)
    v = MSE + penalty
}

// ============================================================================
// _eqt_boot_optimization()
// Constrained OLS optimization using Python COBYLA (matches R)
// Minimizes MSE subject to max|beta_placebo| >= delta
//
// Arguments:
//   X           : real matrix    - double-demeaned design matrix (N x p)
//   Y           : real colvector - double-demeaned dependent variable (N x 1)
//   no_placebos : real scalar    - number of placebo coefficients
//   delta       : real scalar    - equivalence threshold
//   start_val   : real colvector - starting values (usually unconstrained OLS)
//
// Returns:
//   real colvector - constrained coefficients (p x 1)
// ============================================================================
real colvector _eqt_boot_optimization(real matrix X, real colvector Y,
                                      real scalar no_placebos, real scalar delta,
                                      real colvector start_val)
{
    real colvector result
    real scalar max_abs

    // Check if unconstrained solution already satisfies constraint
    max_abs = max(abs(start_val[1::no_placebos]))
    if (max_abs >= delta) {
        return(start_val)
    }

    // Use Python COBYLA implementation for constrained optimization
    result = _eqt_boot_optimization_python(X, Y, no_placebos, delta, start_val)
    if (rows(result) > 0) {
        return(result)
    }

    errprintf("Error: Python constrained OLS unavailable or failed.\n")
    errprintf("Please configure Stata Python and ensure numpy/scipy are installed.\n")
    exit(198)
}

// ============================================================================
// _eqt_boot_optimization_python()
// Python COBYLA wrapper via Stata-Python interface
// Uses _equitrends_constrained_ols.ado which handles dynamic path discovery
// Returns empty vector on failure
// ============================================================================
real colvector _eqt_boot_optimization_python(real matrix X, real colvector Y,
                                             real scalar no_placebos, real scalar delta,
                                             real colvector start_val)
{
    real colvector result
    string scalar x_name, y_name, sv_name, res_name, cmd

    x_name = "__eqt_py_X"
    y_name = "__eqt_py_Y"
    sv_name = "__eqt_py_start"
    res_name = "__eqt_py_result"

    // Store data in Stata matrices for Python access
    st_matrix(x_name, X)
    st_matrix(y_name, Y)
    st_matrix(sv_name, start_val)
    
    // Clear any previous result
    stata("capture matrix drop " + res_name)

    // Build command string for _equitrends_constrained_ols
    // This ADO file handles dynamic Python path discovery via _equitrends_load_python
    cmd = "capture noisily _equitrends_constrained_ols, " +
          "y(" + y_name + ") " +
          "x(" + x_name + ") " +
          "delta(" + strofreal(delta, "%21.15g") + ") " +
          "no_placebos(" + strofreal(no_placebos) + ") " +
          "startval(" + sv_name + ") " +
          "result(" + res_name + ")"
    
    stata(cmd)
    stata("scalar __eqt_rc = _rc")
    
    if (st_numscalar("__eqt_rc") != 0) {
        return(J(0, 1, .))
    }

    // Retrieve result - check if matrix exists
    if (_eqt_matrix_exists(res_name) == 0) {
        return(J(0, 1, .))
    }
    
    result = st_matrix(res_name)
    if (rows(result) == 0) {
        return(J(0, 1, .))
    }

    return(result)
}


// ============================================================================
// _eqt_spherical_bootstrap()
// Perform one spherical bootstrap iteration
// Generates N(0, σ²_c) errors and computes bootstrap max|β̂|
//
// Algorithm:
//   1. Generate u ~ N(0, σ²_c) for each observation
//   2. Construct Y* = Xβ̂_c + u
//   3. Double-demean Y*
//   4. OLS estimate β̂* from demeaned data
//   5. Return max|β̂*[1:no_placebos]|
//
// Arguments:
//   Xb          : real colvector - fitted values Xβ̂_c (N x 1)
//   X           : real matrix    - double-demeaned design matrix (N x p)
//   ID          : real colvector - individual identifiers (N x 1)
//   period      : real colvector - time period identifiers (N x 1)
//   WD          : real matrix    - demeaned time dummy matrix (N x T)
//   variance    : real scalar    - constrained variance σ²_c
//   no_placebos : real scalar    - number of placebo coefficients
//
// Returns:
//   real scalar - max|β̂*[1:no_placebos]|
// ============================================================================
real scalar _eqt_spherical_bootstrap(real colvector Xb, real matrix X,
                                     real colvector ID, real colvector period,
                                     real matrix WD, real scalar variance,
                                     real scalar no_placebos)
{
    real colvector boot_u, boot_y, boot_y_demeaned, boot_coef, placebo_coefs
    real scalar max_abs
    real scalar N
    
    N = rows(X)
    
    // Step 1: Generate N(0, σ²_c) error terms
    // rnormal(n, 1, 0, sqrt(variance)) generates N(0, variance)
    boot_u = rnormal(N, 1, 0, sqrt(variance))
    
    // Step 2: Construct bootstrap Y = Xβ̂_c + u
    boot_y = Xb + boot_u
    
    // Step 3: Double-demean the bootstrap Y
    boot_y_demeaned = _eqt_double_demean(boot_y, ID, period, WD)
    
    // Step 4: OLS estimate bootstrap coefficients
    boot_coef = _eqt_ols_cholesky(cross(X, X), cross(X, boot_y_demeaned))
    
    // Step 5: Extract placebo coefficients and compute max absolute value
    placebo_coefs = boot_coef[1::no_placebos]
    max_abs = max(abs(placebo_coefs))
    
    return(max_abs)
}


// ============================================================================
// _eqt_wild_bootstrap()
// Perform one wild bootstrap iteration with Rademacher weights
// Uses individual-clustered Rademacher variables
//
// The optional cached_id_indices parameter improves efficiency by avoiding
// redundant ID-to-index mapping across bootstrap iterations.
//
// Algorithm:
//   1. Generate Rademacher R_i for each individual (n values)
//   2. Expand R to observation level: R_full[obs] = R[individual[obs]]
//   3. Construct u* = R_full ⊙ û (element-wise product)
//   4. Construct Y* = Xβ̂_c + u*
//   5. Double-demean Y*
//   6. OLS estimate β̂* from demeaned data
//   7. Return max|β̂*[1:no_placebos]|
//
// Arguments:
//   Xb              : real colvector - fitted values Xβ̂_c (N x 1)
//   X               : real matrix    - double-demeaned design matrix (N x p)
//   u_ddot          : real colvector - constrained residuals û (N x 1)
//   ID              : real colvector - individual identifiers (N x 1)
//   period          : real colvector - time period identifiers (N x 1)
//   WD              : real matrix    - demeaned time dummy matrix (N x T)
//   no_placebos     : real scalar    - number of placebo coefficients
//   cached_id_indices: real colvector (optional) - pre-computed ID->index mapping
//
// Returns:
//   real scalar - max|β̂*[1:no_placebos]|
// ============================================================================
real scalar _eqt_wild_bootstrap(real colvector Xb, real matrix X,
                                real colvector u_ddot, real colvector ID,
                                real colvector period, real matrix WD,
                                real scalar no_placebos,
                                | real colvector cached_id_indices)
{
    real colvector unique_ID, R_vars, R_vars_full, id_indices
    real colvector boot_u, boot_y, boot_y_demeaned, boot_coef, placebo_coefs
    real scalar n, N_obs, max_abs, i
    
    N_obs = rows(X)
    unique_ID = uniqrows(ID)
    n = rows(unique_ID)
    
    // Step 1: Generate Rademacher variables for each individual
    R_vars = _eqt_rademacher(n)
    
    // Step 2: Map individual-level R to observation level
    // Use cached mapping if provided, otherwise compute it
    if (args() >= 8 && rows(cached_id_indices) == N_obs) {
        id_indices = cached_id_indices
    }
    else {
        id_indices = _eqt_id_to_index(ID)
    }
    
    // Expand R_vars to full observation level
    R_vars_full = J(N_obs, 1, .)
    for (i = 1; i <= N_obs; i++) {
        R_vars_full[i] = R_vars[id_indices[i]]
    }
    
    // Step 3: Construct bootstrap errors u* = R ⊙ û
    boot_u = R_vars_full :* u_ddot
    
    // Step 4: Construct bootstrap Y = Xβ̂_c + u*
    boot_y = Xb + boot_u
    
    // Step 5: Double-demean the bootstrap Y
    boot_y_demeaned = _eqt_double_demean(boot_y, ID, period, WD)
    
    // Step 6: OLS estimate bootstrap coefficients
    boot_coef = _eqt_ols_cholesky(cross(X, X), cross(X, boot_y_demeaned))
    
    // Step 7: Extract placebo coefficients and compute max absolute value
    placebo_coefs = boot_coef[1::no_placebos]
    max_abs = max(abs(placebo_coefs))
    
    return(max_abs)
}


// ============================================================================
// _eqt_bootstrap_critical_value()
// Compute bootstrap critical value as α quantile of bootstrap distribution
//
// Arguments:
//   bootstrap_maxcoefs : real colvector - bootstrap max|β̂*| values (B x 1)
//   alpha              : real scalar    - significance level (e.g., 0.05)
//
// Returns:
//   real scalar - critical value c*_α
// ============================================================================
real scalar _eqt_bootstrap_critical_value(real colvector bootstrap_maxcoefs,
                                          real scalar alpha)
{
    real scalar critical_value
    
    // Critical value is the α quantile (not 1-α)
    // This is because we reject H0 when max|β̂| < c*_α
    critical_value = _eqt_quantile(bootstrap_maxcoefs, alpha)
    
    return(critical_value)
}



// ============================================================================
// _eqt_bootstrap_test()
// Main bootstrap equivalence test function
// Implements both spherical ("Boot") and wild ("Wild") bootstrap methods
//
// Arguments:
//   Y           : real colvector - dependent variable (N x 1)
//   X           : real matrix    - design matrix (N x p), first no_placebos cols are placebo
//   ID          : real colvector - individual identifiers (N x 1)
//   period      : real colvector - time period identifiers (N x 1)
//   no_placebos : real scalar    - number of placebo coefficients
//   delta       : real scalar    - equivalence threshold
//   alpha       : real scalar    - significance level
//   B           : real scalar    - number of bootstrap replications
//   type        : string scalar  - "Boot" (spherical) or "Wild"
//   seed        : real scalar    - random seed (0 = no seed setting)
//   x_varnames  : string rowvector - (optional) original X variable names for warnings
//
// Returns:
//   struct _eqt_bootstrap_result - complete test results
// ============================================================================
struct _eqt_bootstrap_result scalar _eqt_bootstrap_test(
    real colvector Y, real matrix X, real colvector ID, real colvector period,
    real scalar no_placebos, real scalar delta, real scalar alpha,
    real scalar B, string scalar type, real scalar seed,
    | string rowvector x_varnames)
{
    struct _eqt_bootstrap_result scalar result
    real matrix WD, X_demeaned
    real colvector Y_demeaned, unconstrained_coef, constrained_coef
    real colvector Xb, u_ddot, bootstrap_maxcoefs
    real scalar sigma2_c, max_abs_unconstrained, critical_value
    real scalar b, n, T_periods, N_obs
    real rowvector problematic_wd, problematic_x
    real scalar effective_no_placebos, i, idx
    real rowvector removed_placebo_idx, removed_control_idx
    string scalar period_list
    real colvector unique_periods_ordered, placebo_period_values
    real scalar base_period_val
    real scalar use_external_draws, use_external_bootmax
    real matrix XtX, WtW
    real matrix boot_draws, wild_draws
    real colvector boot_u, boot_y, boot_y_between, boot_y_demeaned, boot_coef, placebo_coefs
    real colvector unique_ID, id_indices, R_vars, R_vars_full
    real scalar max_abs, N_obs_wild, denom, max_rel_diff
    
    // Set random seed if specified
    if (seed > 0) {
        rseed(seed)
    }
    
    // Get dimensions
    N_obs = rows(Y)
    n = rows(uniqrows(ID))
    T_periods = rows(uniqrows(period))
    
    // Step 1: Construct WD matrix (demeaned time dummies)
    WD = _eqt_construct_WD(period, ID)
    
    // Step 1.1: Check for multicollinearity in WD matrix (matches R package)
    // R: if(qr(WD)$rank < ncol(WD)){ WD <- remove_multicollinearity(WD, asmatrix = TRUE)$df }
    // Use QR method to match R's qr() with column pivoting behavior
    problematic_wd = _eqt_detect_multicol_qr(WD)
    if (cols(problematic_wd) > 0) {
        WD = _eqt_remove_multicol(WD, problematic_wd)
    }
    
    // Step 2: Double-demean Y and X
    Y_demeaned = _eqt_double_demean(Y, ID, period, WD)
    X_demeaned = _eqt_double_demean(X, ID, period, WD)
    
    // Map column indices to actual period values for warning messages
    unique_periods_ordered = _equitrends_uniq_preserve_order(period)
    base_period_val = max(period)
    placebo_period_values = select(unique_periods_ordered, 
                                   abs(unique_periods_ordered :- base_period_val) :>= 1e-10)
    
    // Step 2.1: Check for multicollinearity in X_demeaned (matches R package)
    // R: if(qr(X)$rank < ncol(X)){ new_X <- remove_multicollinearity(X) ... }
    // Use QR method to match R's qr() with column pivoting behavior
    effective_no_placebos = no_placebos
    problematic_x = _eqt_detect_multicol_qr(X_demeaned)
    if (cols(problematic_x) > 0) {
        // Separate placebo and control removals
        removed_placebo_idx = J(1, 0, .)
        removed_control_idx = J(1, 0, .)
        
        for (i = 1; i <= cols(problematic_x); i++) {
            idx = problematic_x[i]
            if (idx <= no_placebos) {
                removed_placebo_idx = removed_placebo_idx, idx
            } else {
                removed_control_idx = removed_control_idx, idx
            }
        }
        
        // Display warning for removed placebos (matches R package format)
        // Uses "multicolinearity" spelling to match R package exactly
        if (cols(removed_placebo_idx) > 0) {
            period_list = strofreal(placebo_period_values[removed_placebo_idx[1]])
            for (i = 2; i <= cols(removed_placebo_idx); i++) {
                period_list = period_list + ", " + strofreal(placebo_period_values[removed_placebo_idx[i]])
            }
            printf("{txt}Warning: The placebo corresponding to period(s) %s removed due to multicolinearity.\n", period_list)
        }
        
        // Display warning for removed controls (matches R package format)
        // Display actual variable names if provided
        if (cols(removed_control_idx) > 0) {
            real scalar ctrl_idx
            string scalar var_name_list
            ctrl_idx = removed_control_idx[1] - no_placebos
            if (args() >= 11 & cols(x_varnames) > 0 & ctrl_idx >= 1 & ctrl_idx <= cols(x_varnames)) {
                var_name_list = x_varnames[ctrl_idx]
            } else {
                var_name_list = strofreal(ctrl_idx)
            }
            for (i = 2; i <= cols(removed_control_idx); i++) {
                ctrl_idx = removed_control_idx[i] - no_placebos
                if (args() >= 11 & cols(x_varnames) > 0 & ctrl_idx >= 1 & ctrl_idx <= cols(x_varnames)) {
                    var_name_list = var_name_list + ", " + x_varnames[ctrl_idx]
                } else {
                    var_name_list = var_name_list + ", " + strofreal(ctrl_idx)
                }
            }
            printf("{txt}Warning: The following control variables were removed due to multicolinearity: %s\n", var_name_list)
        }
        
        // Remove problematic columns
        X_demeaned = _eqt_remove_multicol(X_demeaned, problematic_x)
        
        // Update effective number of placebos
        effective_no_placebos = no_placebos - cols(removed_placebo_idx)
    }
    
    // Step 2.2: Check if all placebos were removed due to multicollinearity
    // This is a fatal condition - equivalence test requires at least one placebo
    if (effective_no_placebos <= 0) {
        errprintf("Error: All placebo coefficients were removed due to multicollinearity.\n")
        errprintf("The bootstrap equivalence test cannot be performed without placebo coefficients.\n")
        errprintf("This indicates severe multicollinearity in your data.\n")
        exit(459)
    }
    
    // Step 3: Unconstrained OLS
    unconstrained_coef = _eqt_ols_cholesky(cross(X_demeaned, X_demeaned), 
                                           cross(X_demeaned, Y_demeaned))
    
    // Step 4: Check if constraint is binding
    // Use effective_no_placebos after multicollinearity removal
    max_abs_unconstrained = max(abs(unconstrained_coef[1::effective_no_placebos]))
    
    // Step 5: Constrained OLS (if needed)
    constrained_coef = _eqt_boot_optimization(X_demeaned, Y_demeaned, 
                                              effective_no_placebos, delta, 
                                              unconstrained_coef)
    
    // Step 6: Compute constrained variance estimate
    sigma2_c = _eqt_sigma_hathat_c(constrained_coef, X_demeaned, Y_demeaned, 
                                   ID, period)
    
    // Step 7: Compute fitted values and residuals under constrained model
    Xb = X_demeaned * constrained_coef
    u_ddot = Y_demeaned - Xb
    
    // Step 8: Bootstrap loop
    bootstrap_maxcoefs = J(B, 1, .)
    XtX = cross(X_demeaned, X_demeaned)
    WtW = cross(WD, WD)
    use_external_draws = 0
    if (_eqt_scalar_exists("__eqt_use_external_draws")) {
        use_external_draws = (st_numscalar("__eqt_use_external_draws") == 1)
    }
    use_external_bootmax = 0
    if (_eqt_scalar_exists("__eqt_use_external_bootmax")) {
        use_external_bootmax = (st_numscalar("__eqt_use_external_bootmax") == 1)
    }
    
    // Check if progress display is enabled (controlled by Stata scalar)
    real scalar show_dots
    show_dots = 1  // Default: show progress
    if (_eqt_scalar_exists("__eqt_show_dots")) {
        show_dots = (st_numscalar("__eqt_show_dots") == 1)
    }
    
    // Initialize progress bar
    struct _eqt_progress_bar scalar pb
    _eqt_progress_init(pb, B, "Bootstrap replications (" + type + ")", show_dots)
    
    if (type == "Boot") {
        // Spherical bootstrap
        if (use_external_draws == 1) {
            if (_eqt_matrix_exists("__eqt_boot_draws") == 0) {
                errprintf("External bootstrap draws requested but __eqt_boot_draws not found.\n")
                exit(198)
            }
            boot_draws = st_matrix("__eqt_boot_draws")
            if (rows(boot_draws) != rows(X_demeaned)) {
                errprintf("External bootstrap draws row mismatch: %g vs %g.\n", rows(boot_draws), rows(X_demeaned))
                exit(198)
            }
            if (cols(boot_draws) < B) {
                errprintf("External bootstrap draws column mismatch: %g < %g.\n", cols(boot_draws), B)
                exit(198)
            }
            
            for (b = 1; b <= B; b++) {
                boot_u = boot_draws[., b] * sqrt(sigma2_c)
                boot_y = Xb + boot_u
                // Double-demeaning ensures consistency with the two-way FE specification
                boot_y_demeaned = _eqt_double_demean(boot_y, ID, period, WD)
                boot_coef = _eqt_ols_cholesky(XtX, cross(X_demeaned, boot_y_demeaned))
                placebo_coefs = boot_coef[1::effective_no_placebos]
                max_abs = max(abs(placebo_coefs))
                bootstrap_maxcoefs[b] = max_abs
                
                // Update progress bar
                _eqt_progress_update(pb, b)
            }
        }
        else {
            for (b = 1; b <= B; b++) {
                bootstrap_maxcoefs[b] = _eqt_spherical_bootstrap(
                    Xb, X_demeaned, ID, period, WD, sigma2_c, effective_no_placebos)
                
                // Update progress bar
                _eqt_progress_update(pb, b)
            }
        }
    }
    else if (type == "Wild") {
        // Wild bootstrap
        if (use_external_draws == 1) {
            if (_eqt_matrix_exists("__eqt_wild_draws") == 0) {
                errprintf("External wild draws requested but __eqt_wild_draws not found.\n")
                exit(198)
            }
            wild_draws = st_matrix("__eqt_wild_draws")
            unique_ID = uniqrows(ID)
            if (rows(wild_draws) != rows(unique_ID)) {
                errprintf("External wild draws row mismatch: %g vs %g.\n", rows(wild_draws), rows(unique_ID))
                exit(198)
            }
            if (cols(wild_draws) < B) {
                errprintf("External wild draws column mismatch: %g < %g.\n", cols(wild_draws), B)
                exit(198)
            }
            
            id_indices = _eqt_id_to_index(ID)
            N_obs_wild = rows(X_demeaned)
            
            for (b = 1; b <= B; b++) {
                R_vars = wild_draws[., b]
                R_vars_full = J(N_obs_wild, 1, .)
                for (i = 1; i <= N_obs_wild; i++) {
                    R_vars_full[i] = R_vars[id_indices[i]]
                }
                
                boot_u = R_vars_full :* u_ddot
                boot_y = Xb + boot_u
                // Double-demeaning ensures consistency with the two-way FE specification
                boot_y_demeaned = _eqt_double_demean(boot_y, ID, period, WD)
                boot_coef = _eqt_ols_cholesky(XtX, cross(X_demeaned, boot_y_demeaned))
                placebo_coefs = boot_coef[1::effective_no_placebos]
                max_abs = max(abs(placebo_coefs))
                bootstrap_maxcoefs[b] = max_abs
                
                // Update progress bar
                _eqt_progress_update(pb, b)
            }
        }
        else {
            // Pre-compute ID-to-index mapping for efficiency
            id_indices = _eqt_id_to_index(ID)
            N_obs_wild = rows(X_demeaned)
            unique_ID = uniqrows(ID)
            
            for (b = 1; b <= B; b++) {
                // Generate Rademacher variables for each individual - O(n)
                R_vars = _eqt_rademacher(rows(unique_ID))
                
                // Expand to observation level using precomputed id_indices - O(N)
                R_vars_full = J(N_obs_wild, 1, .)
                for (i = 1; i <= N_obs_wild; i++) {
                    R_vars_full[i] = R_vars[id_indices[i]]
                }
                
                // Bootstrap iteration (same as _eqt_wild_bootstrap)
                boot_u = R_vars_full :* u_ddot
                boot_y = Xb + boot_u
                boot_y_demeaned = _eqt_double_demean(boot_y, ID, period, WD)
                boot_coef = _eqt_ols_cholesky(XtX, cross(X_demeaned, boot_y_demeaned))
                placebo_coefs = boot_coef[1::effective_no_placebos]
                bootstrap_maxcoefs[b] = max(abs(placebo_coefs))
                
                // Update progress bar
                _eqt_progress_update(pb, b)
            }
        }
    }
    else {
        // Invalid type - should not happen if validated upstream
        errprintf("Invalid bootstrap type: %s. Use 'Boot' or 'Wild'.\n", type)
        exit(198)
    }
    
    // Finalize progress bar display
    _eqt_progress_finish(pb)
    
    // Step 8.5: Optional override with external bootstrap maxcoefs
    if (use_external_bootmax == 1) {
        if (type == "Boot") {
            if (_eqt_matrix_exists("__eqt_boot_maxcoefs") == 0) {
                errprintf("External boot maxcoefs requested but __eqt_boot_maxcoefs not found.\n")
                exit(198)
            }
            boot_draws = st_matrix("__eqt_boot_maxcoefs")
            if (rows(boot_draws) != B) {
                errprintf("External boot maxcoefs length mismatch: %g vs %g.\n", rows(boot_draws), B)
                exit(198)
            }
            denom = max(abs(boot_draws))
            if (denom < 1e-12) denom = 1e-12
            max_rel_diff = max(abs(bootstrap_maxcoefs - boot_draws)) / denom
            st_numscalar("__eqt_boot_maxcoefs_rel_diff", max_rel_diff)
            bootstrap_maxcoefs = boot_draws
        }
        else if (type == "Wild") {
            if (_eqt_matrix_exists("__eqt_wild_maxcoefs") == 0) {
                errprintf("External wild maxcoefs requested but __eqt_wild_maxcoefs not found.\n")
                exit(198)
            }
            wild_draws = st_matrix("__eqt_wild_maxcoefs")
            if (rows(wild_draws) != B) {
                errprintf("External wild maxcoefs length mismatch: %g vs %g.\n", rows(wild_draws), B)
                exit(198)
            }
            denom = max(abs(wild_draws))
            if (denom < 1e-12) denom = 1e-12
            max_rel_diff = max(abs(bootstrap_maxcoefs - wild_draws)) / denom
            st_numscalar("__eqt_wild_maxcoefs_rel_diff", max_rel_diff)
            bootstrap_maxcoefs = wild_draws
        }
    }
    
    // Step 9: Compute critical value
    critical_value = _eqt_bootstrap_critical_value(bootstrap_maxcoefs, alpha)
    
    // Step 10: Populate result structure
    // Use effective_no_placebos for extracting placebo coefficients
    result.placebo_coefs = unconstrained_coef[1::effective_no_placebos]
    result.constrained_coefs = constrained_coef
    result.max_abs_coef = max_abs_unconstrained
    result.bootstrap_maxcoefs = bootstrap_maxcoefs
    result.critical_value = critical_value
    result.reject_H0 = (max_abs_unconstrained < critical_value)
    result.sigma_hathat_c = sigma2_c
    result.B = B
    result.alpha = alpha
    result.delta = delta
    result.type = type
    result.n_individuals = n
    result.n_periods = T_periods
    result.n_obs = N_obs
    result.optimization_converged = 1  // Assume converged for now
    result.optimization_iterations = .  // Not tracked in current implementation
    
    return(result)
}


// ============================================================================
// ADO INTERFACE HELPER FUNCTIONS
// ============================================================================

// ============================================================================
// _eqt_ado_constrained_ols()
// Helper function for ado interface - constrained OLS optimization
//
// Arguments:
//   x_name       : string scalar - name of Stata matrix with X
//   y_name       : string scalar - name of Stata matrix with Y
//   no_placebos  : real scalar   - number of placebo coefficients
//   delta        : real scalar   - equivalence threshold
//   startval_name: string scalar - name of Stata matrix with starting values
//   result_name  : string scalar - name for result matrix
//
// Returns: (via Stata matrix)
//   Constrained coefficient estimates
// ============================================================================
void _eqt_ado_constrained_ols(string scalar x_name, string scalar y_name,
                              real scalar no_placebos, real scalar delta,
                              string scalar startval_name, string scalar result_name)
{
    real matrix X
    real colvector Y, start_val, result
    
    // Get matrices from Stata
    X = st_matrix(x_name)
    Y = st_matrix(y_name)
    start_val = st_matrix(startval_name)
    
    // Ensure column vectors
    if (cols(Y) > rows(Y)) Y = Y'
    if (cols(start_val) > rows(start_val)) start_val = start_val'
    
    // Call optimization function
    result = _eqt_boot_optimization(X, Y, no_placebos, delta, start_val)
    
    // Store result in Stata
    st_matrix(result_name, result)
}

end
