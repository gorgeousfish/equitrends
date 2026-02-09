*! equitrends_min_delta.mata - Bootstrap minimum equivalence threshold search
*!
*! Implements numerical optimization to find the minimum equivalence threshold
*! delta* for which the null hypothesis H0: ||beta||_inf >= delta can be rejected.
*! Uses golden section search with a continuous transformation of the discrete
*! hypothesis test outcome to enable standard optimization algorithms.

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// Structure: _eqt_min_delta_result
// Container for minimum threshold search results
// ============================================================================
struct _eqt_min_delta_result {
    real scalar    min_delta           // Minimum equivalence threshold delta*
    real scalar    search_lower        // Search interval lower bound
    real scalar    search_upper        // Search interval upper bound
    real scalar    max_abs_coef        // Maximum absolute placebo coefficient
    real scalar    max_sd              // Maximum standard error
    real colvector placebo_coefs       // Placebo coefficients (unconstrained)
    real scalar    n_iterations        // Optimization iterations performed
    real scalar    converged           // Convergence indicator (1=yes, 0=no)
    real scalar    alpha               // Significance level
    real scalar    B                   // Number of bootstrap replications
    string scalar  type                // Bootstrap type ("Boot" or "Wild")
    real scalar    seed                // Random seed used
}


// ============================================================================
// _eqt_max_sd()
// Compute maximum standard error of placebo coefficients
//
// The search interval upper bound is determined as max|beta| + c * max(se),
// where max(se) is computed from the OLS variance-covariance matrix
// V = sigma^2 * (X'X)^{-1}, extracting the placebo coefficient submatrix.
//
// Arguments:
//   X           : real matrix    - double-demeaned design matrix (N x p)
//   Y           : real colvector - double-demeaned dependent variable (N x 1)
//   ID          : real colvector - individual identifiers (N x 1)
//   period      : real colvector - time period identifiers (N x 1)
//   no_placebos : real scalar    - number of placebo coefficients T
//
// Returns:
//   real scalar - max_{l in {1,...,T}} sqrt(V_ll)
// ============================================================================
real scalar _eqt_max_sd(real matrix X, real colvector Y,
                        real colvector ID, real colvector period,
                        real scalar no_placebos)
{
    real colvector unconstrained_coefs
    real scalar placebo_variance
    real matrix XtX_inv, vcov_mat, vcov_mat_placebos
    real colvector placebo_se
    real scalar max_sd
    
    // Compute OLS coefficients: beta = (X'X)^{-1}X'Y
    unconstrained_coefs = _eqt_ols_cholesky(cross(X, X), cross(X, Y))
    
    // Compute residual variance sigma^2 = RSS / df
    placebo_variance = _eqt_sigma_hathat_c(unconstrained_coefs, X, Y, ID, period)
    
    // Compute variance-covariance matrix: V = sigma^2 * (X'X)^{-1}
    XtX_inv = invsym(cross(X, X))
    vcov_mat = XtX_inv * placebo_variance
    
    // Extract placebo coefficient standard errors and return maximum
    vcov_mat_placebos = vcov_mat[1::no_placebos, 1::no_placebos]
    placebo_se = sqrt(diagonal(vcov_mat_placebos))
    max_sd = max(placebo_se)
    
    return(max_sd)
}


// ============================================================================
// _eqt_min_delta_wrapper()
// Transform hypothesis test outcome to continuous optimization objective
//
// Converts the binary reject/not-reject decision into a continuous value:
//   f(delta) = -exp(-delta)  if H0 is rejected at delta
//   f(delta) =  exp(-delta)  if H0 is not rejected at delta
//
// This ensures the minimum of f corresponds to the smallest delta where
// rejection occurs, enabling standard optimization algorithms.
//
// Arguments:
//   delta       : real scalar    - equivalence threshold being evaluated
//   Y           : real colvector - dependent variable (N x 1)
//   X           : real matrix    - design matrix (N x p)
//   ID          : real colvector - individual identifiers (N x 1)
//   period      : real colvector - time period identifiers (N x 1)
//   no_placebos : real scalar    - number of placebo coefficients
//   alpha       : real scalar    - significance level
//   B           : real scalar    - number of bootstrap replications
//   type        : string scalar  - "Boot" (parametric) or "Wild" (cluster)
//   seed        : real scalar    - random seed (0 = use current state)
//
// Returns:
//   real scalar - transformed objective value
// ============================================================================
real scalar _eqt_min_delta_wrapper(real scalar delta,
                                   real colvector Y, real matrix X,
                                   real colvector ID, real colvector period,
                                   real scalar no_placebos, real scalar alpha,
                                   real scalar B, string scalar type,
                                   real scalar seed)
{
    struct _eqt_bootstrap_result scalar boot_result
    real scalar wrapper_value
    
    // Evaluate bootstrap test at current threshold
    boot_result = _eqt_bootstrap_test(Y, X, ID, period, no_placebos, 
                                       delta, alpha, B, type, seed)
    
    // Transform to continuous objective for optimization
    if (boot_result.reject_H0) {
        wrapper_value = -exp(-delta)
    }
    else {
        wrapper_value = exp(-delta)
    }
    
    return(wrapper_value)
}


// ============================================================================
// _eqt_golden_search()
// One-dimensional optimization via golden section search
//
// Iteratively narrows the search interval [a,b] using the golden ratio
// phi = (sqrt(5) - 1)/2, guaranteeing convergence to a local minimum.
//
// Arguments:
//   a, b        : real scalars   - search interval bounds
//   tol         : real scalar    - convergence tolerance
//   max_iter    : real scalar    - maximum iterations
//   Y           : real colvector - dependent variable (N x 1)
//   X           : real matrix    - design matrix (N x p)
//   ID          : real colvector - individual identifiers
//   period      : real colvector - time period identifiers
//   no_placebos : real scalar    - number of placebo coefficients
//   alpha       : real scalar    - significance level
//   B           : real scalar    - bootstrap replications
//   type        : string scalar  - bootstrap type
//   seed        : real scalar    - random seed
//   n_iter      : real scalar    - (output) iterations performed
//   converged   : real scalar    - (output) convergence indicator
//
// Returns:
//   real scalar - minimizer delta*
// ============================================================================
real scalar _eqt_golden_search(real scalar a, real scalar b,
                               real scalar tol, real scalar max_iter,
                               real colvector Y, real matrix X,
                               real colvector ID, real colvector period,
                               real scalar no_placebos, real scalar alpha,
                               real scalar B, string scalar type,
                               real scalar seed,
                               real scalar n_iter, real scalar converged)
{
    real scalar phi, x1, x2, f1, f2
    real scalar iter
    
    // Golden ratio phi = (sqrt(5) - 1) / 2
    phi = (sqrt(5) - 1) / 2
    
    // Initialize interior evaluation points
    x1 = a + (1 - phi) * (b - a)
    x2 = a + phi * (b - a)
    f1 = _eqt_min_delta_wrapper(x1, Y, X, ID, period, no_placebos, alpha, B, type, seed)
    f2 = _eqt_min_delta_wrapper(x2, Y, X, ID, period, no_placebos, alpha, B, type, seed)
    
    // Iteratively narrow the search interval
    converged = 0
    
    // Check if progress display is enabled
    real scalar show_search_progress
    show_search_progress = 1
    if (_eqt_scalar_exists("__eqt_show_search_progress")) {
        show_search_progress = (st_numscalar("__eqt_show_search_progress") == 1)
    }
    
    for (iter = 1; iter <= max_iter; iter++) {
        if (abs(b - a) < tol) {
            converged = 1
            break
        }
        
        // Update search progress display
        _eqt_progress_search_update(iter, max_iter, abs(b - a), show_search_progress)
        
        if (f1 < f2) {
            b = x2
            x2 = x1
            f2 = f1
            x1 = a + (1 - phi) * (b - a)
            f1 = _eqt_min_delta_wrapper(x1, Y, X, ID, period, no_placebos, alpha, B, type, seed)
        }
        else {
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + phi * (b - a)
            f2 = _eqt_min_delta_wrapper(x2, Y, X, ID, period, no_placebos, alpha, B, type, seed)
        }
    }
    
    // Finalize search progress display
    _eqt_progress_search_finish(iter, converged, show_search_progress)
    
    n_iter = iter
    return((a + b) / 2)
}



// ============================================================================
// _eqt_min_delta_bootstrap()
// Find minimum equivalence threshold delta* via bootstrap hypothesis inversion
//
// Computes the smallest delta such that H0: ||beta||_inf >= delta can be
// rejected. The search interval is [max|beta|, max|beta| + 15*max(se)].
//
// Arguments:
//   Y           : real colvector - dependent variable (N x 1)
//   X           : real matrix    - design matrix (N x p)
//   ID          : real colvector - individual identifiers (N x 1)
//   period      : real colvector - time period identifiers (N x 1)
//   no_placebos : real scalar    - number of placebo coefficients T
//   alpha       : real scalar    - significance level
//   B           : real scalar    - number of bootstrap replications
//   type        : string scalar  - "Boot" (parametric) or "Wild" (cluster)
//   seed        : real scalar    - random seed (0 = use current state)
//
// Returns:
//   struct _eqt_min_delta_result - search results including delta*
// ============================================================================
struct _eqt_min_delta_result scalar _eqt_min_delta_bootstrap(
    real colvector Y, real matrix X, real colvector ID, real colvector period,
    real scalar no_placebos, real scalar alpha, real scalar B,
    string scalar type, real scalar seed)
{
    struct _eqt_min_delta_result scalar result
    real matrix WD, X_demeaned
    real colvector Y_demeaned, unconstrained_coefs, placebo_coefs
    real scalar max_abs_coef, max_sd, lower, upper
    real scalar min_delta, n_iter, converged
    real scalar tol, max_iter
    
    // Optimization parameters
    tol = 1e-6
    max_iter = 100
    
    // Construct WD matrix and double-demean data
    WD = _eqt_construct_WD(period, ID)
    Y_demeaned = _eqt_double_demean(Y, ID, period, WD)
    X_demeaned = _eqt_double_demean(X, ID, period, WD)
    
    // Handle multicollinearity in design matrices
    real rowvector problematic_wd, problematic_x
    real scalar effective_no_placebos, i, idx
    real rowvector removed_placebo_idx
    
    // Check WD matrix for multicollinearity
    problematic_wd = _eqt_detect_multicol_qr(WD)
    if (cols(problematic_wd) > 0) {
        WD = _eqt_remove_multicol(WD, problematic_wd)
        Y_demeaned = _eqt_double_demean(Y, ID, period, WD)
        X_demeaned = _eqt_double_demean(X, ID, period, WD)
    }
    
    // Check demeaned design matrix for multicollinearity
    effective_no_placebos = no_placebos
    problematic_x = _eqt_detect_multicol_qr(X_demeaned)
    if (cols(problematic_x) > 0) {
        removed_placebo_idx = J(1, 0, .)
        for (i = 1; i <= cols(problematic_x); i++) {
            idx = problematic_x[i]
            if (idx <= no_placebos) {
                removed_placebo_idx = removed_placebo_idx, idx
            }
        }
        X_demeaned = _eqt_remove_multicol(X_demeaned, problematic_x)
        X = _eqt_remove_multicol(X, problematic_x)
        effective_no_placebos = no_placebos - cols(removed_placebo_idx)
    }
    
    // Compute unconstrained OLS coefficients
    unconstrained_coefs = _eqt_ols_cholesky(cross(X_demeaned, X_demeaned), 
                                            cross(X_demeaned, Y_demeaned))
    placebo_coefs = unconstrained_coefs[1::effective_no_placebos]
    
    // Determine search interval bounds
    max_abs_coef = max(abs(placebo_coefs))
    max_sd = _eqt_max_sd(X_demeaned, Y_demeaned, ID, period, effective_no_placebos)
    lower = max_abs_coef
    upper = max_abs_coef + 15 * max_sd
    
    // Ensure valid search interval
    if (lower < 0) lower = 0
    if (upper <= lower) upper = lower + 1
    
    // Initialize random number generator if seed specified
    if (seed > 0) {
        rseed(seed)
    }
    
    // Disable progress display during golden section search
    // The search involves multiple bootstrap calls; displaying progress for each
    // would clutter the output. Instead, display a single message.
    real scalar saved_show_dots, show_dots_existed
    show_dots_existed = _eqt_scalar_exists("__eqt_show_dots")
    if (show_dots_existed) {
        saved_show_dots = st_numscalar("__eqt_show_dots")
    }
    else {
        saved_show_dots = 1
    }
    
    // Display search progress message if progress display is enabled
    if (saved_show_dots == 1) {
        _eqt_progress_search_init(max_iter, 
            "Searching for minimum equivalence threshold", saved_show_dots)
    }
    
    // Disable progress display for internal bootstrap calls
    st_numscalar("__eqt_show_dots", 0)
    
    // Enable search-level progress display
    st_numscalar("__eqt_show_search_progress", saved_show_dots)
    
    // Perform golden section search
    n_iter = 0
    converged = 0
    min_delta = _eqt_golden_search(lower, upper, tol, max_iter,
                                   Y, X, ID, period,
                                   effective_no_placebos, alpha, B, type, 0,
                                   n_iter, converged)
    
    // Restore original progress display setting
    if (show_dots_existed) {
        st_numscalar("__eqt_show_dots", saved_show_dots)
    }
    else {
        stata("capture scalar drop __eqt_show_dots")
    }
    
    // Clean up search progress scalar
    stata("capture scalar drop __eqt_show_search_progress")
    
    // Populate result structure
    result.min_delta = min_delta
    result.search_lower = lower
    result.search_upper = upper
    result.max_abs_coef = max_abs_coef
    result.max_sd = max_sd
    result.placebo_coefs = placebo_coefs
    result.n_iterations = n_iter
    result.converged = converged
    result.alpha = alpha
    result.B = B
    result.type = type
    result.seed = seed
    
    return(result)
}


// ============================================================================
// _eqt_min_delta_bootstrap_simple()
// Simplified interface returning only the minimum threshold delta*
// ============================================================================
real scalar _eqt_min_delta_bootstrap_simple(
    real colvector Y, real matrix X, real colvector ID, real colvector period,
    real scalar no_placebos, real scalar alpha, real scalar B,
    string scalar type, real scalar seed)
{
    struct _eqt_min_delta_result scalar result
    
    result = _eqt_min_delta_bootstrap(Y, X, ID, period, no_placebos, 
                                       alpha, B, type, seed)
    
    return(result.min_delta)
}


// ============================================================================
// ADO Interface Functions
// ============================================================================

// _eqt_ado_max_sd() - Stata interface for maximum standard error calculation
void _eqt_ado_max_sd(string scalar X_name, string scalar Y_name,
                     string scalar ID_name, string scalar period_name,
                     real scalar no_placebos)
{
    real matrix X
    real colvector Y, ID, period
    real scalar max_sd
    
    // Get data from Stata
    X = st_matrix(X_name)
    Y = st_matrix(Y_name)
    if (cols(Y) > rows(Y)) Y = Y'
    
    ID = st_data(., ID_name)
    period = st_data(., period_name)
    
    // Calculate max_sd
    max_sd = _eqt_max_sd(X, Y, ID, period, no_placebos)
    
    // Return to Stata
    st_numscalar("r(max_sd)", max_sd)
}


// ============================================================================
// _eqt_ado_min_delta_bootstrap()
// Stata interface for bootstrap minimum threshold search
// ============================================================================
void _eqt_ado_min_delta_bootstrap(string scalar Y_name, string scalar X_names,
                                   string scalar ID_name, string scalar period_name,
                                   real scalar no_placebos, real scalar alpha,
                                   real scalar B, string scalar type,
                                   real scalar seed)
{
    struct _eqt_min_delta_result scalar result
    real matrix X
    real colvector Y, ID, period
    string rowvector X_varnames
    real scalar i
    
    // Parse X variable names
    X_varnames = tokens(X_names)
    
    // Get data from Stata
    Y = st_data(., Y_name)
    X = st_data(., X_varnames)
    ID = st_data(., ID_name)
    period = st_data(., period_name)
    
    // Run minimum threshold search
    result = _eqt_min_delta_bootstrap(Y, X, ID, period, no_placebos,
                                       alpha, B, type, seed)
    
    // Return results to Stata
    st_numscalar("r(min_delta)", result.min_delta)
    st_numscalar("r(max_abs_coef)", result.max_abs_coef)
    st_numscalar("r(max_sd)", result.max_sd)
    st_numscalar("r(search_lower)", result.search_lower)
    st_numscalar("r(search_upper)", result.search_upper)
    st_numscalar("r(n_iterations)", result.n_iterations)
    st_numscalar("r(converged)", result.converged)
    st_numscalar("r(alpha)", result.alpha)
    st_numscalar("r(B)", result.B)
    
    // Return placebo coefficients as matrix
    st_matrix("r(placebo_coefs)", result.placebo_coefs)
    st_matrixrowstripe("r(placebo_coefs)", (J(no_placebos, 1, ""), strofreal(1::no_placebos)))
    st_matrixcolstripe("r(placebo_coefs)", ("", "coef"))
}

end
