*! _equitrends_constrained_ols.ado - Constrained OLS with infinity-norm constraint
*!
*! Computes constrained OLS estimates subject to ||beta||_inf = delta for the
*! bootstrap test of hypothesis H0: ||beta||_inf >= delta (maximum placebo effect).
*! The constrained estimator is used to generate bootstrap samples under the null.

capture program drop _equitrends_constrained_ols
program define _equitrends_constrained_ols, rclass
    version 16.0
    
    // =========================================================================
    // _equitrends_constrained_ols
    // Constrained OLS estimation with infinity-norm constraint
    //
    // The sum of squared residuals is minimized subject to:
    //   ||beta||_inf = max_{l=1,...,T} |beta_l| = delta
    //
    // The constrained estimator beta_c satisfies the null hypothesis boundary.
    // Bootstrap samples are then generated from:
    //   Y_it^{(b)} = W_it' * beta_c + u_it^{(b)}
    //
    // Optimization is performed via COBYLA (derivative-free constrained
    // optimization algorithm) through the Python/SciPy interface.
    //
    // Arguments:
    //   y(name)          : N x 1 response vector (double-demeaned)
    //   x(name)          : N x p design matrix (double-demeaned)
    //   delta(real)      : equivalence threshold (positive)
    //   no_placebos(int) : number of placebo coefficients (T)
    //   startval(name)   : p x 1 initial values for optimization
    //   result(name)     : name for output matrix (p x 1)
    //
    // Returns:
    //   r(result)                       : constrained coefficient estimates
    //   eqt_constrained_ols_converged   : convergence flag (scalar, 1 if converged)
    // =========================================================================
    
    // -------------------------------------------------------------------------
    // Parse syntax and validate inputs
    // -------------------------------------------------------------------------
    syntax , Y(name) X(name) Delta(real) No_placebos(integer) ///
             Startval(name) Result(name)
    
    if `delta' <= 0 {
        display as error "Error: delta must be strictly positive, got `delta'"
        exit 198
    }
    
    if `no_placebos' < 1 {
        display as error "Error: no_placebos must be at least 1, got `no_placebos'"
        exit 198
    }
    
    capture confirm matrix `y'
    if _rc != 0 {
        display as error "Error: matrix `y' not found"
        exit 198
    }
    
    capture confirm matrix `x'
    if _rc != 0 {
        display as error "Error: matrix `x' not found"
        exit 198
    }
    
    capture confirm matrix `startval'
    if _rc != 0 {
        display as error "Error: matrix `startval' not found"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Validate matrix dimensions
    // -------------------------------------------------------------------------
    local y_rows = rowsof(`y')
    local y_cols = colsof(`y')
    local x_rows = rowsof(`x')
    local x_cols = colsof(`x')
    local sv_rows = rowsof(`startval')
    local sv_cols = colsof(`startval')
    
    if `y_cols' != 1 {
        display as error "Error: Y must be a column vector (N x 1), got `y_rows' x `y_cols'"
        exit 198
    }
    
    if `y_rows' != `x_rows' {
        display as error "Error: Y and X must have same number of rows"
        display as error "       Y has `y_rows' rows, X has `x_rows' rows"
        exit 198
    }
    
    // Starting values are accepted as p x 1 or 1 x p (flattened internally)
    local sv_len = max(`sv_rows', `sv_cols')
    if `sv_len' != `x_cols' {
        display as error "Error: startval length must equal X columns"
        display as error "       startval has `sv_len' elements, X has `x_cols' columns"
        exit 198
    }
    
    if `no_placebos' > `x_cols' {
        display as error "Error: no_placebos cannot exceed number of coefficients"
        display as error "       no_placebos=`no_placebos', X has `x_cols' columns"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Solve constrained optimization via COBYLA (Python/SciPy)
    // -------------------------------------------------------------------------
    // Python environment is initialized and dependencies are imported
    capture _equitrends_load_python, quiet
    if _rc != 0 {
        display as error "Error: Python module unavailable for constrained OLS."
        exit 198
    }
    
    capture python: from sfi import Matrix, Scalar
    if _rc != 0 {
        display as error "Error: Python SFI bridge unavailable."
        exit 603
    }

    capture python: import numpy as np
    if _rc != 0 {
        display as error "Error: numpy not available in Python."
        exit 603
    }

    capture python: import equitrends_python as etp
    if _rc != 0 {
        display as error "Error: equitrends_python not importable."
        exit 603
    }

    // Data matrices are transferred from Stata to Python
    capture python: Y_np = np.array(Matrix.get("`y'")).flatten()
    if _rc != 0 {
        display as error "Error: Failed to load Y matrix into Python."
        exit 603
    }

    capture python: X_np = np.array(Matrix.get("`x'"))
    if _rc != 0 {
        display as error "Error: Failed to load X matrix into Python."
        exit 603
    }

    capture python: sv_np = np.array(Matrix.get("`startval'")).flatten()
    if _rc != 0 {
        display as error "Error: Failed to load start values into Python."
        exit 603
    }

    // RSS is minimized subject to ||beta[1:no_placebos]||_inf = delta
    capture python: beta_py, conv = etp.constrained_ols(Y_np, X_np, `delta', `no_placebos', sv_np)
    if _rc != 0 {
        display as error "Error: Python constrained OLS failed."
        exit 603
    }

    // Results are transferred back to Stata
    capture python: Matrix.store("`result'", beta_py.reshape(-1, 1))
    if _rc != 0 {
        display as error "Error: Failed to store constrained coefficients."
        exit 603
    }

    capture python: Scalar.setValue("eqt_constrained_ols_converged", 1 if conv else 0)
    if _rc != 0 {
        display as error "Error: Failed to store convergence flag."
        exit 603
    }
    
    // -------------------------------------------------------------------------
    // Store and return constrained estimates
    // -------------------------------------------------------------------------
    capture confirm matrix `result'
    if _rc == 0 {
        tempname result_copy
        matrix `result_copy' = `result'
        return matrix `result' = `result_copy'
    }
    else {
        display as error "Error: Result matrix not created"
        exit 603
    }
end
