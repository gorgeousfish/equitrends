*! _equitrends_demean.ado - Double demeaning transformations for TWFE estimation
*!
*! Implements the double demeaning transformation for two-way fixed effects
*! (TWFE) model estimation. The transformation removes individual and time
*! fixed effects from panel data:
*!
*!   Y_ddot_{i,t} = Y_{i,t} - Y_bar_{i,.} - Y_bar_{.,t} + Y_bar_{..}
*!
*! where Y_bar_{i,.} is the individual mean, Y_bar_{.,t} is the time mean,
*! and Y_bar_{..} is the grand mean. The transformed variables are used for
*! pooled OLS estimation of placebo treatment effects beta.

version 16.0

// ============================================================================
// Main Program - Subcommand dispatcher
// ============================================================================
program define _equitrends_demean
    
    gettoken subcmd 0 : 0
    
    if "`subcmd'" == "double_demean" {
        _equitrends_demean_dd `0'
    }
    else if "`subcmd'" == "construct_wd" {
        _equitrends_demean_wd `0'
    }
    else if "`subcmd'" == "sigma_c" {
        _equitrends_demean_sigma `0'
    }
    else if "`subcmd'" == "grouped_mean" {
        _equitrends_demean_gm `0'
    }
    else if "`subcmd'" == "between_transform" {
        _equitrends_demean_bt `0'
    }
    else {
        di as error "Unknown subcommand: `subcmd'"
        di as error "Valid subcommands: double_demean, construct_wd, sigma_c, grouped_mean, between_transform"
        exit 198
    }
end

// ============================================================================
// _equitrends_demean_dd
// Apply double demeaning transformation for TWFE estimation
//
// The two-way fixed effects transformation removes individual and time
// fixed effects: X_ddot = X - X_bar_i - X_bar_t + X_bar
//
// Arguments:
//   x      : string - Mata matrix name to be demeaned (N x p)
//   id     : string - Mata vector of individual identifiers (N x 1)
//   time   : string - Mata vector of time period indicators (N x 1)
//   wd     : string - Mata matrix name for W_ddot (interface compatibility)
//   result : string - Output Mata matrix name for demeaned data (N x p)
// ============================================================================
program define _equitrends_demean_dd
    version 16.0
    syntax, x(string) id(string) time(string) wd(string) result(string)
    
    // Validate input matrix existence
    mata: st_local("x_exists", strofreal(st_matrix("`x'") != J(0,0,.)))
    if "`x_exists'" == "0" {
        di as error "Matrix `x' not found"
        exit 198
    }
    
    // Apply double demeaning transformation
    mata: st_matrix("`result'", _eqt_double_demean( ///
        st_matrix("`x'"), ///
        st_matrix("`id'"), ///
        st_matrix("`time'"), ///
        st_matrix("`wd'") ///
    ))
end

// ============================================================================
// _equitrends_demean_wd
// Construct demeaned treatment indicator matrix W_ddot
//
// Constructs the matrix W_ddot where W_{i,t,l} = G_i * D_l(t) represents the
// interaction between treatment group membership G_i and time period dummies
// D_l(t). Individual means are removed: W_ddot_{i,t} = W_{i,t} - W_bar_i
//
// Arguments:
//   period : string - Mata vector of time period indicators (N x 1)
//   id     : string - Mata vector of individual identifiers (N x 1)
//   result : string - Output Mata matrix name for W_ddot (N x T)
// ============================================================================
program define _equitrends_demean_wd
    version 16.0
    syntax, period(string) id(string) result(string)
    
    mata: st_matrix("`result'", _eqt_construct_WD( ///
        st_matrix("`period'"), ///
        st_matrix("`id'") ///
    ))
end

// ============================================================================
// _equitrends_demean_sigma
// Compute constrained residual variance for bootstrap inference
//
// Computes the constrained variance estimate sigma_c^2 for the parametric
// bootstrap test of the maximum equivalence hypothesis. The variance is
// estimated from residuals under the constraint ||beta||_inf = delta.
//
// Formula: sigma_c^2 = sum((Y_ddot - W_ddot' * beta_c)^2) / ((n-1) * T)
//
// Arguments:
//   param  : string - Mata vector of constrained beta_c (p x 1)
//   x      : string - Mata design matrix W_ddot (N x p)
//   y      : string - Mata outcome vector Y_ddot (N x 1)
//   id     : string - Mata vector of individual identifiers (N x 1)
//   time   : string - Mata vector of time period indicators (N x 1)
//   result : string - Output scalar name for sigma_c^2
// ============================================================================
program define _equitrends_demean_sigma
    version 16.0
    syntax, param(string) x(string) y(string) id(string) time(string) result(string)
    
    mata: st_numscalar("`result'", _eqt_sigma_hathat_c( ///
        st_matrix("`param'"), ///
        st_matrix("`x'"), ///
        st_matrix("`y'"), ///
        st_matrix("`id'"), ///
        st_matrix("`time'") ///
    ))
end

// ============================================================================
// _equitrends_demean_gm
// Compute group means for panel transformations
//
// Computes group means and expands them to a vector of the same length as
// the input, where each element is replaced by its corresponding group mean.
//
// Arguments:
//   x      : string - Mata vector of values (N x 1)
//   group  : string - Mata vector of group identifiers (N x 1)
//   result : string - Output Mata vector for expanded group means (N x 1)
// ============================================================================
program define _equitrends_demean_gm
    version 16.0
    syntax, x(string) group(string) result(string)
    
    mata: st_matrix("`result'", _eqt_grouped_mean( ///
        st_matrix("`x'"), ///
        st_matrix("`group'") ///
    ))
end

// ============================================================================
// _equitrends_demean_bt
// Apply within-group demeaning transformation
//
// Removes group-specific means: x_demeaned = x - x_bar_group
//
// This transformation removes group-level fixed effects and serves as a
// component of the double demeaning procedure.
//
// Arguments:
//   x      : string - Mata vector of values (N x 1)
//   group  : string - Mata vector of group identifiers (N x 1)
//   result : string - Output Mata vector for demeaned values (N x 1)
// ============================================================================
program define _equitrends_demean_bt
    version 16.0
    syntax, x(string) group(string) result(string)
    
    mata: st_matrix("`result'", _eqt_between_trans( ///
        st_matrix("`x'"), ///
        st_matrix("`group'") ///
    ))
end
