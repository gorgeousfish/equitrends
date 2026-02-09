*! equitrends_multicollin.mata - Multicollinearity detection and handling
*!
*! This file contains functions for multicollinearity detection and removal:
*!   - _eqt_detect_multicol()      : Detect collinear columns using SVD
*!   - _eqt_detect_multicol_qr()   : Detect using QR (alternative)
*!   - _eqt_remove_multicol()      : Remove collinear columns from matrix
*!   - _eqt_warn_placebo_rm()      : Display warning for removed placebos
*!   - _eqt_warn_control_rm()      : Display warning for removed controls
*!   - _eqt_process_multicol()     : Complete multicollinearity workflow

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// _eqt_detect_multicol()
// Detect multicollinearity in a matrix using SVD-based rank detection
//
// Uses SVD for rank detection and identifies columns that contribute least
// to the column space.
//
// Arguments:
//   X   : real matrix - input matrix (N_obs x p)
//   tol : real scalar - tolerance for rank detection (default 1e-10)
//
// Returns:
//   real rowvector - indices of problematic columns (1-based), empty if no collinearity
// ============================================================================
real rowvector _eqt_detect_multicol(real matrix X, | real scalar tol)
{
    real matrix U, Vt, V
    real colvector s
    real scalar r, ncols, i
    real rowvector problematic_vars
    real colvector col_contributions, sorted_idx
    real scalar n_remove
    
    // Default tolerance
    if (args() < 2) {
        tol = 1e-10
    }
    
    ncols = cols(X)
    
    // Handle edge cases
    if (ncols == 0) {
        return(J(1, 0, .))
    }
    if (ncols == 1) {
        // Single column - check if all zeros
        if (max(abs(X)) < tol) {
            return(1)
        }
        return(J(1, 0, .))
    }
    
    // Compute SVD
    svd(X, U, s, Vt)
    V = Vt'
    
    // Determine numerical rank
    r = sum(s :> (max(s) * tol * max((rows(X), ncols))))
    
    // If full rank, no multicollinearity
    if (r >= ncols) {
        return(J(1, 0, .))
    }
    
    // Number of columns to remove
    n_remove = ncols - r
    
    // Identify which columns to remove
    // We look at the contribution of each column to the null space
    // Columns with largest contribution to null space singular vectors are problematic
    
    // Compute column contributions to null space
    col_contributions = J(ncols, 1, 0)
    for (i = r + 1; i <= ncols; i++) {
        if (i <= rows(s)) {
            col_contributions = col_contributions + abs(V[., i]):^2
        }
    }
    
    // If all contributions are zero (shouldn't happen), use last columns
    if (max(col_contributions) < tol) {
        problematic_vars = (ncols - n_remove + 1)..ncols
        return(problematic_vars)
    }
    
    // Sort columns by contribution (descending) and take top n_remove
    sorted_idx = order(-col_contributions, 1)
    // sorted_idx is a column vector, extract first n_remove elements and convert to rowvector
    problematic_vars = sort(sorted_idx[1..n_remove], 1)'
    
    return(problematic_vars)
}

// ============================================================================
// _eqt_detect_multicol_qr()
// Implementation using Modified Gram-Schmidt to match R's qr() LINPACK behavior
//
// R's default qr() (LINPACK) processes columns in original order without
// aggressive maximum-norm pivoting. It identifies linearly dependent columns
// as those whose norm becomes negligible after orthogonalization.
//
// Algorithm (Modified Gram-Schmidt WITHOUT column pivoting):
// 1. Process columns in original order (k = 1, 2, ..., ncol)
// 2. For each column k:
//    a. Orthogonalize against all previous columns (1..k-1)
//    b. Check if orthogonalized norm is below tolerance
//    c. If yes, mark column k as linearly dependent
//    d. If no, add to orthonormal basis
// 3. Return indices of all linearly dependent columns
//
// This matches R's behavior where qr(X)$pivot for matrix with col3=col1+col2
// returns [1,2,3] (no reordering), rank=2, and problematic=pivot[3]=3.
//
// Arguments:
//   X   : real matrix - input matrix (N_obs x p)
//   tol : real scalar - tolerance for rank detection (default 1e-10)
//
// Returns:
//   real rowvector - indices of problematic columns (1-based), empty if no collinearity
// ============================================================================
real rowvector _eqt_detect_multicol_qr(real matrix X, | real scalar tol)
{
    real matrix Q
    real colvector x_orth
    real rowvector problematic_vars
    real scalar ncols, nrows, k, j, n_basis
    real scalar col_norm, max_norm, rel_tol, proj
    
    // Default tolerance (matches R's behavior)
    if (args() < 2) {
        tol = 1e-10
    }
    
    ncols = cols(X)
    nrows = rows(X)
    
    // Handle edge cases
    if (ncols == 0) {
        return(J(1, 0, .))
    }
    if (ncols == 1) {
        // Single column - check if all zeros
        if (max(abs(X)) < tol) {
            return(1)
        }
        return(J(1, 0, .))
    }
    
    // Compute maximum column norm for relative tolerance
    max_norm = 0
    for (j = 1; j <= ncols; j++) {
        col_norm = sqrt(cross(X[., j], X[., j]))
        if (col_norm > max_norm) {
            max_norm = col_norm
        }
    }
    
    if (max_norm < tol) {
        // All columns are essentially zero
        return(1..ncols)
    }
    
    // Relative tolerance (matches R's LINPACK behavior)
    rel_tol = tol * max_norm * max((nrows, ncols))
    
    // Initialize orthonormal basis storage
    // Q will store the orthonormal basis vectors (only for independent columns)
    Q = J(nrows, ncols, 0)
    n_basis = 0
    
    // Track problematic (linearly dependent) columns
    problematic_vars = J(1, 0, .)
    
    // Process columns in ORIGINAL order (matching R's LINPACK behavior)
    for (k = 1; k <= ncols; k++) {
        // Start with original column k
        x_orth = X[., k]
        
        // Orthogonalize against all existing basis vectors
        for (j = 1; j <= n_basis; j++) {
            proj = cross(Q[., j], x_orth)
            x_orth = x_orth - proj * Q[., j]
        }
        
        // Compute norm of orthogonalized column
        col_norm = sqrt(cross(x_orth, x_orth))
        
        // Check if column is linearly dependent (norm near zero after orthogonalization)
        if (col_norm < rel_tol) {
            // Column k is linearly dependent on previous columns
            problematic_vars = problematic_vars, k
        }
        else {
            // Column k is linearly independent - add to basis
            n_basis = n_basis + 1
            Q[., n_basis] = x_orth / col_norm
        }
    }
    
    return(problematic_vars)
}

// ============================================================================
// _eqt_remove_multicol()
// Remove multicollinear columns from a matrix
//
// Arguments:
//   X                : real matrix - input matrix (N_obs x p)
//   problematic_vars : real rowvector - indices of columns to remove (1-based)
//
// Returns:
//   real matrix - matrix with problematic columns removed
// ============================================================================
real matrix _eqt_remove_multicol(real matrix X, real rowvector problematic_vars)
{
    real rowvector keep_cols
    real scalar j, ncols, n_prob
    real scalar is_problematic
    
    n_prob = cols(problematic_vars)
    
    // If no problematic columns, return original matrix
    if (n_prob == 0) {
        return(X)
    }
    
    ncols = cols(X)
    keep_cols = J(1, 0, .)
    
    // Build list of columns to keep
    for (j = 1; j <= ncols; j++) {
        is_problematic = 0
        if (anyof(problematic_vars, j)) {
            is_problematic = 1
        }
        if (!is_problematic) {
            keep_cols = keep_cols, j
        }
    }
    
    // Return matrix with only kept columns
    if (cols(keep_cols) == 0) {
        return(J(rows(X), 0, .))
    }
    
    return(X[., keep_cols])
}

// ============================================================================
// _eqt_warn_placebo_rm()
// Display warning message when placebo variables are removed due to multicollinearity
//
// Note: Uses "multicolinearity" spelling to match R package exactly
//
// Arguments:
//   periods : string vector - period names of removed placebos
// ============================================================================
void _eqt_warn_placebo_rm(string vector periods)
{
    string scalar period_list, msg
    real scalar i
    
    if (cols(periods) == 0) {
        return
    }
    
    // Build comma-separated list of periods
    period_list = periods[1]
    for (i = 2; i <= cols(periods); i++) {
        period_list = period_list + ", " + periods[i]
    }
    
    // Display warning (matches R package format exactly)
    msg = "The placebo corresponding to period(s) " + period_list + " removed due to multicolinearity."
    printf("{txt}Warning: %s\n", msg)
}

// ============================================================================
// _eqt_warn_control_rm()
// Display warning message when control variables are removed due to multicollinearity
//
// Note: Uses "multicolinearity" spelling to match R package exactly
//
// Arguments:
//   varnames : string vector - names of removed control variables
// ============================================================================
void _eqt_warn_control_rm(string vector varnames)
{
    string scalar var_list, msg
    real scalar i
    
    if (cols(varnames) == 0) {
        return
    }
    
    // Build comma-separated list of variable names
    var_list = varnames[1]
    for (i = 2; i <= cols(varnames); i++) {
        var_list = var_list + ", " + varnames[i]
    }
    
    // Display warning (matches R package format exactly)
    msg = "The following control variables were removed due to multicolinearity: " + var_list
    printf("{txt}Warning: %s\n", msg)
}

// ============================================================================
// _eqt_process_multicol()
// Complete multicollinearity processing workflow
//
// This function:
// 1. Detects multicollinearity in the design matrix
// 2. Identifies which columns are problematic
// 3. Separates placebo and control variable removals
// 4. Displays appropriate warnings
// 5. Returns the cleaned matrix
//
// Arguments:
//   X             : real matrix - design matrix (N_obs x p)
//   no_placebos   : real scalar - number of placebo columns (first no_placebos cols)
//   placebo_names : string vector - names of placebo variables
//   control_names : string vector - names of control variables
//
// Returns:
//   struct containing:
//     - X_clean : cleaned matrix
//     - removed_placebo_idx : indices of removed placebo columns
//     - removed_control_idx : indices of removed control columns
//     - new_no_placebos : updated number of placebos
// ============================================================================
struct _eqt_multicol_result {
    real matrix X_clean
    real rowvector removed_placebo_idx
    real rowvector removed_control_idx
    real scalar new_no_placebos
    string vector new_placebo_names
    string vector new_control_names
}

struct _eqt_multicol_result scalar _eqt_process_multicol(
    real matrix X,
    real scalar no_placebos,
    string vector placebo_names,
    string vector control_names
)
{
    struct _eqt_multicol_result scalar result
    real rowvector problematic_vars
    real rowvector removed_placebo_idx, removed_control_idx
    string vector removed_placebo_names, removed_control_names
    real scalar i, idx, n_prob
    real rowvector keep_placebo_idx, keep_control_idx
    
    // Initialize result
    result.removed_placebo_idx = J(1, 0, .)
    result.removed_control_idx = J(1, 0, .)
    
    // Detect multicollinearity
    problematic_vars = _eqt_detect_multicol(X)
    n_prob = cols(problematic_vars)
    
    // If no multicollinearity, return original
    if (n_prob == 0) {
        result.X_clean = X
        result.new_no_placebos = no_placebos
        result.new_placebo_names = placebo_names
        result.new_control_names = control_names
        return(result)
    }
    
    // Separate placebo and control removals
    removed_placebo_idx = J(1, 0, .)
    removed_control_idx = J(1, 0, .)
    removed_placebo_names = J(1, 0, "")
    removed_control_names = J(1, 0, "")
    
    for (i = 1; i <= n_prob; i++) {
        idx = problematic_vars[i]
        if (idx <= no_placebos) {
            removed_placebo_idx = removed_placebo_idx, idx
            if (idx <= cols(placebo_names)) {
                removed_placebo_names = removed_placebo_names, placebo_names[idx]
            }
        } else {
            removed_control_idx = removed_control_idx, idx
            if (idx - no_placebos <= cols(control_names)) {
                removed_control_names = removed_control_names, control_names[idx - no_placebos]
            }
        }
    }
    
    // Display warnings
    if (cols(removed_placebo_names) > 0) {
        // Extract period numbers from placebo names (remove "placebo_" prefix)
        string vector period_nums
        period_nums = J(1, cols(removed_placebo_names), "")
        for (i = 1; i <= cols(removed_placebo_names); i++) {
            period_nums[i] = subinstr(removed_placebo_names[i], "placebo_", "")
        }
        _eqt_warn_placebo_rm(period_nums)
    }
    
    if (cols(removed_control_names) > 0) {
        _eqt_warn_control_rm(removed_control_names)
    }
    
    // Remove problematic columns
    result.X_clean = _eqt_remove_multicol(X, problematic_vars)
    
    // Update placebo count
    result.new_no_placebos = no_placebos - cols(removed_placebo_idx)
    
    // Update names
    keep_placebo_idx = J(1, 0, .)
    for (i = 1; i <= no_placebos; i++) {
        if (!anyof(removed_placebo_idx, i)) {
            keep_placebo_idx = keep_placebo_idx, i
        }
    }
    if (cols(keep_placebo_idx) > 0 & cols(placebo_names) > 0) {
        result.new_placebo_names = placebo_names[keep_placebo_idx]
    } else {
        result.new_placebo_names = J(1, 0, "")
    }
    
    keep_control_idx = J(1, 0, .)
    for (i = no_placebos + 1; i <= cols(X); i++) {
        if (!anyof(removed_control_idx, i)) {
            keep_control_idx = keep_control_idx, i - no_placebos
        }
    }
    if (cols(keep_control_idx) > 0 & cols(control_names) > 0) {
        result.new_control_names = control_names[keep_control_idx]
    } else {
        result.new_control_names = J(1, 0, "")
    }
    
    // Store removed indices
    result.removed_placebo_idx = removed_placebo_idx
    result.removed_control_idx = removed_control_idx
    
    return(result)
}

end
