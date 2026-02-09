*! equitrends_validate.mata - Input validation functions for EQUITRENDS
*!
*! This file contains validation functions for the EQUITRENDS package:
*!   - _equitrends_is_colvector()         : Check if matrix is a column vector
*!   - _equitrends_check_lengths()        : Check all vectors have equal length
*!   - _equitrends_is_binary()            : Check if vector contains only 0/1
*!   - _equitrends_is_subset()            : Check if A is subset of B
*!   - _equitrends_base_in_pretreatment() : Check if base_period is in pretreatment
*!   - _equitrends_count_unique()         : Count unique values in vector
*!   - _equitrends_val_equiv_thresh() : Validate equiv_threshold >= 0
*!   - _equitrends_validate_alpha_range() : Validate 0 < alpha < 1
*!   - _equitrends_validate_type()        : Validate type parameter
*!   - _equitrends_validate_method()      : Validate method parameter
*!   - _equitrends_validate_B()           : Validate bootstrap B parameter

version 16.0

mata:
mata set matastrict on

// ============================================================================
// _equitrends_is_colvector()
// Check if a matrix is a column vector (n x 1)
//
// Arguments:
//   x - transmorphic matrix to check
//
// Returns:
//   real scalar - 1 if column vector, 0 otherwise
//
// Notes:
//   - Empty matrix (0 x 0) returns 0
//   - Row vector (1 x n where n > 1) returns 0
//   - Scalar (1 x 1) returns 1 (treated as column vector)
// ============================================================================
real scalar _equitrends_is_colvector(transmorphic matrix x)
{
    // Check if it's a column vector: cols == 1 and rows >= 1
    if (cols(x) == 1 & rows(x) >= 1) {
        return(1)
    }
    return(0)
}

// ============================================================================
// _equitrends_check_lengths()
// Check that all provided vectors/matrices have equal length (rows)
//
// Arguments:
//   Y       - real colvector, dependent variable
//   ID      - real colvector, individual identifier
//   G       - real colvector, treatment indicator
//   period  - real colvector, time period
//   X       - real matrix, control variables (optional, can be empty)
//   cluster - real colvector, cluster variable (optional, can be empty)
//
// Returns:
//   real scalar - 1 if all lengths match, 0 otherwise
//
// Notes:
//   - X is checked by rows (number of observations)
//   - Empty X (0 columns) is ignored
//   - Empty cluster (0 rows) is ignored
// ============================================================================
real scalar _equitrends_check_lengths(
    real colvector Y,
    real colvector ID,
    real colvector G,
    real colvector period,
    | real matrix X,
    real colvector cluster
)
{
    real scalar n_Y, n_ID, n_G, n_period
    
    n_Y = rows(Y)
    n_ID = rows(ID)
    n_G = rows(G)
    n_period = rows(period)
    
    // Check basic vectors
    if (n_Y != n_ID | n_Y != n_G | n_Y != n_period) {
        return(0)
    }
    
    // Check X if provided and non-empty
    if (args() >= 5 & cols(X) > 0) {
        if (rows(X) != n_Y) {
            return(0)
        }
    }
    
    // Check cluster if provided and non-empty
    if (args() >= 6 & rows(cluster) > 0) {
        if (rows(cluster) != n_Y) {
            return(0)
        }
    }
    
    return(1)
}

// ============================================================================
// _equitrends_is_binary()
// Check if a vector contains only binary values (0 or 1)
//
// Arguments:
//   G - real colvector to check
//
// Returns:
//   real scalar - 1 if all values are 0 or 1, 0 otherwise
//
// Notes:
//   - Returns 0 if any missing values are present
//   - Empty vector returns 1 (vacuously true)
// ============================================================================
real scalar _equitrends_is_binary(real colvector G)
{
    real scalar i, n, val
    
    n = rows(G)
    
    // Empty vector is vacuously binary
    if (n == 0) {
        return(1)
    }
    
    // Check each element
    for (i = 1; i <= n; i++) {
        val = G[i]
        
        // Check for missing
        if (missing(val)) {
            return(0)
        }
        
        // Check if 0 or 1
        if (val != 0 & val != 1) {
            return(0)
        }
    }
    
    return(1)
}

// ============================================================================
// _equitrends_is_subset()
// Check if all elements of A are contained in B
//
// Arguments:
//   A - real colvector, the subset to check
//   B - real colvector, the superset
//
// Returns:
//   real scalar - 1 if A ⊆ B, 0 otherwise
//
// Notes:
//   - Empty A is always a subset (vacuously true)
//   - Tolerance of 1e-10 is used for floating point comparison
// ============================================================================
real scalar _equitrends_is_subset(real colvector A, real colvector B)
{
    real scalar i, j, n_A, n_B, found, tol
    
    tol = 1e-10
    n_A = rows(A)
    n_B = rows(B)
    
    // Empty A is always a subset
    if (n_A == 0) {
        return(1)
    }
    
    // Empty B cannot contain non-empty A
    if (n_B == 0) {
        return(0)
    }
    
    // Check each element of A
    for (i = 1; i <= n_A; i++) {
        found = 0
        for (j = 1; j <= n_B; j++) {
            if (abs(A[i] - B[j]) < tol) {
                found = 1
                break
            }
        }
        if (!found) {
            return(0)
        }
    }
    
    return(1)
}

// ============================================================================
// _equitrends_base_in_pretreatment()
// Check if base_period is an element of pretreatment_period
//
// Arguments:
//   base_period         - real scalar, the base period
//   pretreatment_period - real colvector, the pretreatment periods
//
// Returns:
//   real scalar - 1 if base_period ∈ pretreatment_period, 0 otherwise
//
// Notes:
//   - Uses tolerance 1e-10 for floating point comparison
//   - Empty pretreatment_period returns 0
// ============================================================================
real scalar _equitrends_base_in_pretreatment(
    real scalar base_period,
    real colvector pretreatment_period
)
{
    real scalar i, n, tol
    
    tol = 1e-10
    n = rows(pretreatment_period)
    
    // Empty pretreatment_period cannot contain base_period
    if (n == 0) {
        return(0)
    }
    
    // Search for base_period
    for (i = 1; i <= n; i++) {
        if (abs(pretreatment_period[i] - base_period) < tol) {
            return(1)
        }
    }
    
    return(0)
}

// ============================================================================
// _equitrends_count_unique()
// Count the number of unique values in a vector
//
// Arguments:
//   x - real colvector
//
// Returns:
//   real scalar - number of unique values
//
// Notes:
//   - Uses uniqrows() for counting
//   - Empty vector returns 0
// ============================================================================
real scalar _equitrends_count_unique(real colvector x)
{
    if (rows(x) == 0) {
        return(0)
    }
    return(rows(uniqrows(x)))
}

// ============================================================================
// _equitrends_val_equiv_thresh()
// Validate that equivalence threshold is strictly positive
//
// Arguments:
//   equiv_threshold - real scalar
//
// Returns:
//   real scalar - 1 if valid (> 0), 0 if invalid (<= 0)
//
// Notes:
//   - The equivalence threshold delta must be strictly positive
//   - Mathematical requirement: if delta = 0, H0: ||beta||_inf >= 0 is
//     always satisfied since norms are non-negative, rendering the test
//     vacuous and incapable of rejecting the null hypothesis
//   - Error message: "equiv_threshold must be strictly positive (>0)."
// ============================================================================
real scalar _equitrends_val_equiv_thresh(real scalar equiv_threshold)
{
    if (equiv_threshold <= 0) {
        return(0)
    }
    return(1)
}

// ============================================================================
// _equitrends_validate_alpha_range()
// Validate that alpha is strictly between 0 and 1
//
// Arguments:
//   alpha - real scalar
//
// Returns:
//   real scalar - 1 if valid (0 < alpha < 1), 0 otherwise
//
// Notes:
//   - The significance level must satisfy 0 < alpha < 1
// ============================================================================
real scalar _equitrends_validate_alpha_range(real scalar alpha)
{
    if (alpha <= 0 | alpha >= 1) {
        return(0)
    }
    return(1)
}

// ============================================================================
// _equitrends_validate_type()
// Validate that type parameter is one of: max, mean, rms
//
// Arguments:
//   type - string scalar
//
// Returns:
//   real scalar - 1 if valid, 0 otherwise
//
// Notes:
//   - Case-sensitive comparison
//   - Error message: "type must be max, mean, or rms"
// ============================================================================
real scalar _equitrends_validate_type(string scalar type)
{
    if (type == "max" | type == "mean" | type == "rms") {
        return(1)
    }
    return(0)
}

// ============================================================================
// _equitrends_validate_method()
// Validate that method parameter is one of: IU, Boot, Wild
//
// Arguments:
//   method - string scalar
//
// Returns:
//   real scalar - 1 if valid, 0 otherwise
//
// Notes:
//   - Case-sensitive comparison
//   - Only applicable when type == "max"
//   - Error message: "method must be IU, Boot, or Wild"
// ============================================================================
real scalar _equitrends_validate_method(string scalar method)
{
    if (method == "IU" | method == "Boot" | method == "Wild") {
        return(1)
    }
    return(0)
}

// ============================================================================
// _equitrends_validate_B()
// Validate that B is a strictly positive integer
//
// Arguments:
//   B - real scalar
//
// Returns:
//   real scalar - 1 if valid (positive integer), 0 otherwise
//
// Notes:
//   - B must be > 0 and an integer
//   - Error message: "B must be a strictly positive integer scalar"
// ============================================================================
real scalar _equitrends_validate_B(real scalar B)
{
    // Check if positive
    if (B <= 0) {
        return(0)
    }
    
    // Check if integer (no fractional part)
    if (B != floor(B)) {
        return(0)
    }
    
    return(1)
}

// ============================================================================
// _equitrends_is_scalar()
// Check if a value is a scalar (1 x 1 matrix)
//
// Arguments:
//   x - transmorphic matrix
//
// Returns:
//   real scalar - 1 if scalar, 0 otherwise
// ============================================================================
real scalar _equitrends_is_scalar(transmorphic matrix x)
{
    if (rows(x) == 1 & cols(x) == 1) {
        return(1)
    }
    return(0)
}

// ============================================================================
// _equitrends_has_X_columns()
// Check if X matrix has at least one column
//
// Arguments:
//   X - real matrix
//
// Returns:
//   real scalar - 1 if cols(X) >= 1, 0 otherwise
//
// Notes:
//   - Error message: "X must have at least one column."
// ============================================================================
real scalar _equitrends_has_X_columns(real matrix X)
{
    if (cols(X) >= 1) {
        return(1)
    }
    return(0)
}

end
