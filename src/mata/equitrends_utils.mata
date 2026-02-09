*! equitrends_utils.mata - Utility functions for EQUITRENDS
*!
*! This file contains utility functions used by the EQUITRENDS package:
*!   - _equitrends_norm_inf()        : Infinity norm (max absolute value)
*!   - _equitrends_rms()             : Root mean square
*!   - _equitrends_mean()            : Mean of vector
*!   - _equitrends_validate_threshold() : Validate equivalence threshold
*!   - _equitrends_validate_alpha()  : Validate significance level
*!   - _equitrends_W_critical_value(): W distribution critical value lookup
*!   - _equitrends_validate_rms_alpha(): Validate alpha for RMS test

version 16.0

mata:
mata set matastrict on

// ============================================================================
// _equitrends_norm_inf()
// Compute the infinity norm (maximum absolute value) of a vector
//
// Arguments:
//   x - real colvector
//
// Returns:
//   real scalar - max(|x_i|)
//
// Notes:
//   - Returns 0 for empty vectors
//   - Used for Max equivalence test statistic
// ============================================================================
real scalar _equitrends_norm_inf(real colvector x)
{
    if (rows(x) == 0) {
        return(0)
    }
    return(max(abs(x)))
}

// ============================================================================
// _equitrends_rms()
// Compute the root mean square of a vector
//
// Arguments:
//   x - real colvector with T elements
//
// Returns:
//   real scalar - sqrt(sum(x.^2)/T)
//
// Notes:
//   - Returns 0 for empty vectors
//   - Used for RMS equivalence test statistic
// ============================================================================
real scalar _equitrends_rms(real colvector x)
{
    real scalar T
    
    T = rows(x)
    if (T == 0) {
        return(0)
    }
    return(sqrt(sum(x:^2) / T))
}

// ============================================================================
// _equitrends_mean()
// Compute the mean of a vector
//
// Arguments:
//   x - real colvector
//
// Returns:
//   real scalar - mean(x)
//
// Notes:
//   - Wrapper around Mata's built-in mean() function
//   - Returns missing for empty vectors
// ============================================================================
real scalar _equitrends_mean(real colvector x)
{
    return(mean(x))
}

// ============================================================================
// _equitrends_validate_threshold()
// Validate that equivalence threshold is positive
//
// Arguments:
//   delta     - real scalar, the equivalence threshold
//   test_type - string scalar, name of the test (for error message)
//
// Returns:
//   void (exits with error if invalid)
//
// Notes:
//   - Threshold must be strictly positive
//   - Calls _error(3498) on failure
// ============================================================================
void _equitrends_validate_threshold(real scalar delta, string scalar test_type)
{
    if (delta <= 0) {
        errprintf("Error in %s test: Equivalence threshold must be positive\n", test_type)
        errprintf("  Got: %g\n", delta)
        _error(3498)
    }
}

// ============================================================================
// _equitrends_validate_alpha()
// Validate that significance level is in (0, 1)
//
// Arguments:
//   alpha - real scalar, the significance level
//
// Returns:
//   void (exits with error if invalid)
//
// Notes:
//   - Alpha must be strictly between 0 and 1
//   - Calls _error(3498) on failure
// ============================================================================
void _equitrends_validate_alpha(real scalar alpha)
{
    if (alpha <= 0 | alpha >= 1) {
        errprintf("Error: Significance level alpha must be in (0, 1)\n")
        errprintf("  Got: %g\n", alpha)
        _error(3498)
    }
}

// ============================================================================
// _equitrends_W_critical_value()
// Look up critical value from W distribution
//
// Arguments:
//   alpha - real scalar, the significance level
//
// Returns:
//   real scalar - W distribution quantile at alpha
//
// Supported alpha values:
//   0.005, 0.01, 0.025, 0.05, 0.1, 0.2 (left tail, for RMS test)
//   0.8, 0.9, 0.95, 0.975, 0.99, 0.995 (right tail, for confidence intervals)
//
// Notes:
//   - Uses tolerance 1e-10 for floating point comparison
//   - Calls _error(3498) for unsupported alpha values
// ============================================================================
real scalar _equitrends_W_critical_value(real scalar alpha)
{
    real matrix crit_table
    real scalar i, tol
    
    // Tolerance for floating point comparison
    tol = 1e-10
    
    // W distribution critical values obtained via Monte Carlo simulation
    // Format: (alpha, Q_W(alpha))
    // Note: 0.005 and 0.995 are extrapolated from the W distribution tail behavior
    crit_table = (
        0.005, -5.3000000 \
        0.010, -4.2329959 \
        0.025, -2.9047318 \
        0.050, -2.1431720 \
        0.100, -1.4601327 \
        0.200, -0.8561188 \
        0.800,  0.8344549 \
        0.900,  1.4928501 \
        0.950,  2.2003307 \
        0.975,  3.0065821 \
        0.990,  4.1150156 \
        0.995,  5.2000000
    )
    
    // Search for matching alpha
    for (i = 1; i <= rows(crit_table); i++) {
        if (abs(crit_table[i, 1] - alpha) < tol) {
            return(crit_table[i, 2])
        }
    }
    
    // Alpha not found - report error
    errprintf("Error: alpha must be one of: 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, ")
    errprintf("0.8, 0.9, 0.95, 0.975, 0.99, 0.995\n")
    errprintf("  Got: %g\n", alpha)
    _error(3498)
    
    // This line is never reached, but needed for compiler
    return(.)
}

// ============================================================================
// _equitrends_validate_rms_alpha()
// Validate that alpha is supported for RMS equivalence test
//
// Arguments:
//   alpha - real scalar, the significance level
//
// Returns:
//   real scalar - 1 if valid, 0 if invalid
//
// Valid alpha values for RMS test:
//   0.01, 0.025, 0.05, 0.1, 0.2
//
// Notes:
//   - RMS test only supports specific alpha values due to W distribution
//   - Displays error message if invalid (but does not call _error)
// ============================================================================
real scalar _equitrends_validate_rms_alpha(real scalar alpha)
{
    real rowvector valid_alphas
    real scalar tol, i
    
    tol = 1e-10
    valid_alphas = (0.01, 0.025, 0.05, 0.1, 0.2)
    
    for (i = 1; i <= cols(valid_alphas); i++) {
        if (abs(valid_alphas[i] - alpha) < tol) {
            return(1)
        }
    }
    
    // Invalid alpha - display error message
    errprintf("Error: For RMS equivalence test, alpha must be one of: ")
    errprintf("0.01, 0.025, 0.05, 0.1, 0.2\n")
    errprintf("  Got: %g\n", alpha)
    
    return(0)
}


// ============================================================================
// _eqt_omit_na()
// Omit rows with missing values from data vectors/matrices
//
// Arguments:
//   Y       : real colvector - outcome variable (required)
//   ID      : real colvector - individual ID (required)
//   G       : real colvector - treatment group indicator (required)
//   period  : real colvector - time period (required)
//   X       : real matrix    - control variables (optional, can be empty)
//   cluster : real colvector - cluster variable (optional, can be empty)
//
// Returns:
//   real colvector - indices of valid rows (1-indexed)
//
// Notes:
//   - Returns indices of rows where all specified variables are non-missing
//   - For matrix X, a row is valid if no element in that row is missing
//   - Empty optional arguments (0 rows/cols) are ignored
// ============================================================================
real colvector _eqt_omit_na(real colvector Y, real colvector ID,
                             real colvector G, real colvector period,
                             | real matrix X, real colvector cluster)
{
    real colvector valid_rows
    real scalar n, has_X, has_cluster
    
    n = rows(Y)
    if (n == 0) {
        return(J(0, 1, .))
    }
    
    // Initialize all rows as valid
    valid_rows = J(n, 1, 1)
    
    // Check required variables for missing values
    valid_rows = valid_rows :& !rowmissing(Y)
    valid_rows = valid_rows :& !rowmissing(ID)
    valid_rows = valid_rows :& !rowmissing(G)
    valid_rows = valid_rows :& !rowmissing(period)
    
    // Check optional X matrix (if provided and non-empty)
    has_X = (args() >= 5) & (cols(X) > 0) & (rows(X) > 0)
    if (has_X) {
        // rowmissing(X) returns the count of missing values in each row
        // A row is valid if this count is 0
        valid_rows = valid_rows :& (rowmissing(X) :== 0)
    }
    
    // Check optional cluster variable (if provided and non-empty)
    has_cluster = (args() >= 6) & (rows(cluster) > 0)
    if (has_cluster) {
        valid_rows = valid_rows :& !rowmissing(cluster)
    }
    
    // Return indices of valid rows
    return(selectindex(valid_rows))
}


// ============================================================================
// _eqt_apply_na_filter()
// Apply missing value filter to data vectors and matrices
// Helper function to apply the filter indices from _eqt_omit_na()
//
// Arguments:
//   valid_idx : real colvector - indices of valid rows from _eqt_omit_na()
//   Y         : real colvector - outcome variable
//   ID        : real colvector - individual ID
//   G         : real colvector - treatment group indicator
//   period    : real colvector - time period
//   X         : real matrix    - control variables (optional)
//   cluster   : real colvector - cluster variable (optional)
//
// Returns:
//   real scalar - number of rows omitted due to missing values
//
// Notes:
//   - Returns the count of omitted observations for warning display
//   - Actual filtering is performed by the caller
// ============================================================================
real scalar _eqt_apply_na_filter(real colvector valid_idx,
                                  real colvector Y, real colvector ID,
                                  real colvector G, real colvector period,
                                  | real matrix X, real colvector cluster)
{
    real scalar n_orig, n_valid, na_omitted
    real scalar has_X, has_cluster
    
    n_orig = rows(Y)
    n_valid = rows(valid_idx)
    na_omitted = n_orig - n_valid
    
    // If no rows omitted, nothing to do
    if (na_omitted == 0) {
        return(0)
    }
    
    // This function cannot modify in place since Mata passes by value
    // The caller must handle the filtering
    // This function just returns the count for the warning message
    
    return(na_omitted)
}

end
