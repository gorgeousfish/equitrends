*! _equitrends_check_deps.ado - Dependency verification for equitrends package
*!
*! Verifies availability of required dependencies. Python with
*! numpy and scipy is required for constrained optimization in the bootstrap
*! test for hypothesis (3.1).

version 16.0

// ============================================================================
// _equitrends_check_python
// Verifies availability of Python environment and required packages
//
// Python with numpy and scipy is required for constrained least squares
// optimization under the constraint ||beta||_inf = delta in the bootstrap
// test for hypothesis (3.1).
//
// Arguments:
//   quiet : option - Suppresses display output when specified
//
// Returns:
//   r(python_available) : scalar - 1 if Python is available, 0 otherwise
//   r(numpy_available)  : scalar - 1 if numpy is available, 0 otherwise
//   r(scipy_available)  : scalar - 1 if scipy is available, 0 otherwise
// ============================================================================
capture program drop _equitrends_check_python
program define _equitrends_check_python, rclass
    version 16.0
    syntax [, Quiet]
    
    local is_quiet = ("`quiet'" != "")
    
    capture python query
    if _rc != 0 {
        if !`is_quiet' {
            display as error "Error: Python is not available in Stata."
            display as text "Configure with: python set exec <path>"
        }
        return scalar python_available = 0
        return scalar numpy_available = 0
        return scalar scipy_available = 0
        exit
    }
    
    local numpy_ok = 1
    capture python: import numpy
    if _rc != 0 {
        local numpy_ok = 0
        if !`is_quiet' {
            display as error "Error: numpy is not installed."
        }
    }
    
    local scipy_ok = 1
    capture python: import scipy
    if _rc != 0 {
        local scipy_ok = 0
        if !`is_quiet' {
            display as error "Error: scipy is not installed."
        }
    }
    
    if !`is_quiet' {
        display as text "  python: " as result "OK"
        if `numpy_ok' {
            display as text "  numpy: " as result "OK"
        }
        else {
            display as text "  numpy: " as error "Missing"
        }
        if `scipy_ok' {
            display as text "  scipy: " as result "OK"
        }
        else {
            display as text "  scipy: " as error "Missing"
        }
    }
    
    return scalar python_available = 1
    return scalar numpy_available = `numpy_ok'
    return scalar scipy_available = `scipy_ok'
end

// ============================================================================
// _equitrends_check_deps
// Main entry point for complete dependency verification
//
// Performs verification of all required dependencies.
// Required dependencies (Python, numpy, scipy) must be satisfied for
// constrained optimization in bootstrap tests.
//
// Arguments:
//   quiet : option - Suppresses display output when specified
//
// Returns:
//   r(python_available)      : scalar - 1 if Python available
//   r(numpy_available)       : scalar - 1 if numpy available
//   r(scipy_available)       : scalar - 1 if scipy available
//   r(all_dependencies_ok)   : scalar - 1 if all required deps satisfied
// ============================================================================
capture program drop _equitrends_check_deps
program define _equitrends_check_deps, rclass
    version 16.0
    syntax [, Quiet]
    
    local all_ok = 1
    local is_quiet = ("`quiet'" != "")
    
    if !`is_quiet' {
        display as text ""
        display as text "============================================"
        display as text "EQUITRENDS Dependency Check"
        display as text "============================================"
        display as text ""
        display as text "Required Dependencies (for bootstrap methods):"
    }
    
    // -------------------------------------------------------------------------
    // Required dependencies: Python environment with numpy and scipy
    // -------------------------------------------------------------------------
    capture noisily _equitrends_check_python, `quiet'
    if _rc != 0 {
        local python_avail = 0
        local numpy_avail = 0
        local scipy_avail = 0
    }
    else {
        local python_avail = r(python_available)
        local numpy_avail = r(numpy_available)
        local scipy_avail = r(scipy_available)
    }
    if `python_avail' == 0 | `numpy_avail' == 0 | `scipy_avail' == 0 {
        local all_ok = 0
    }
    
    return scalar python_available = `python_avail'
    return scalar numpy_available = `numpy_avail'
    return scalar scipy_available = `scipy_avail'
    return scalar all_dependencies_ok = `all_ok'
    
    if !`is_quiet' {
        display as text ""
        display as text "============================================"
    }
    
    if `all_ok' == 0 {
        if !`is_quiet' {
            display as error "Status: Some required dependencies are missing."
            display as text "============================================"
            display as text ""
        }
        exit 198
    }
    else {
        if !`is_quiet' {
            display as result "Status: All dependencies satisfied."
            display as text "============================================"
            display as text ""
        }
    }
end
