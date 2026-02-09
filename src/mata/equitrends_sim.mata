*! equitrends_sim.mata - Panel data simulation for Monte Carlo studies
*!
*! Generates panel data following the data generating process in Dette and 
*! Schumann (2024). Supports AR(p) error processes, heteroskedasticity, 
*! Ashenfelter's dip, and linear trend violations of the parallel trends 
*! assumption.

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// FUNCTION: _eqt_sim_ar_process()
// Generate AR(p) process with burn-in period
//
// Arguments:
//   n_periods  - Number of periods to generate
//   phi        - AR coefficients (1 x p row vector)
//   sigma      - Innovation standard deviation
//   burnins    - Number of burn-in periods to discard
//
// Returns:
//   Column vector of length n_periods containing AR(p) process values
//
// Algorithm:
//   1. Initialize first p values with N(0, sigma) random draws
//   2. Generate (burnins + n_periods) values using AR recursion:
//      u_t = sum_{k=1}^{p} phi_k * u_{t-k} + eta_t, where eta_t ~ N(0, sigma^2)
//   3. Discard first burnins values, return last n_periods values
// ============================================================================
real colvector _eqt_sim_ar_process(
    real scalar n_periods,
    real rowvector phi,
    real scalar sigma,
    real scalar burnins
)
{
    real scalar p, total_periods, t, k
    real colvector u, eta
    real scalar ar_sum
    
    // Get AR order from phi length
    p = cols(phi)
    
    // Handle special case: phi = 0 (white noise)
    if (p == 1 && phi[1] == 0) {
        return(rnormal(n_periods, 1, 0, sigma))
    }
    
    // Total periods including burn-in
    total_periods = burnins + n_periods
    
    // Initialize AR process vector
    u = J(total_periods, 1, 0)
    
    // Generate innovations
    eta = rnormal(total_periods, 1, 0, sigma)
    
    // Initialize first p values with innovations
    for (t = 1; t <= p; t++) {
        u[t] = eta[t]
    }
    
    // Generate AR(p) process
    for (t = p + 1; t <= total_periods; t++) {
        ar_sum = 0
        for (k = 1; k <= p; k++) {
            ar_sum = ar_sum + phi[k] * u[t - k]
        }
        u[t] = ar_sum + eta[t]
    }
    
    // Return last n_periods values (discard burn-in)
    return(u[(burnins + 1)..total_periods])
}


// ============================================================================
// FUNCTION: _eqt_sim_treatment()
// Assign treatment group indicator
//
// Arguments:
//   n       - Number of individuals
//   rcompat - R compatibility mode flag (0 = paper mode, 1 = R compat mode)
//
// Returns:
//   Column vector of length n with treatment indicators (0 or 1)
//
// Algorithm:
//   Paper mode (rcompat=0): G_i = 1 if U_i < 0.5, where U_i ~ U(0,1)
//   R compat mode (rcompat=1): G_i = 1 if mod(i, 2) == 0 (even IDs)
// ============================================================================
real colvector _eqt_sim_treatment(
    real scalar n,
    real scalar rcompat
)
{
    real colvector G, U
    real scalar i
    
    G = J(n, 1, 0)
    
    if (rcompat == 1) {
        // R compatibility mode: even IDs are treated
        for (i = 1; i <= n; i++) {
            G[i] = (mod(i, 2) == 0)
        }
    }
    else {
        // Paper mode: random assignment with Pr(G=1) = 0.5
        U = runiform(n, 1)
        G = (U :< 0.5)
    }
    
    return(G)
}


// ============================================================================
// FUNCTION: _eqt_sim_fixed_effects()
// Generate or validate fixed effects
//
// Arguments:
//   n         - Number of individuals
//   n_periods - Number of time periods (T+2)
//   eta       - Input individual fixed effects (n x 1) or empty
//   lambda    - Input time fixed effects ((T+2) x 1) or empty
//   rcompat   - R compatibility mode flag
//
// Returns:
//   Pointer to struct containing eta and lambda vectors
//   (Modified in place via pointer arguments)
//
// Algorithm:
//   If eta is empty and rcompat=0: eta_i ~ N(0,1)
//   If eta is empty and rcompat=1: eta_i = 0
//   If lambda is empty and rcompat=0: lambda_t ~ N(0,1)
//   If lambda is empty and rcompat=1: lambda_t = 0
// ============================================================================
void _eqt_sim_fixed_effects(
    real scalar n,
    real scalar n_periods,
    real colvector eta,
    real colvector lambda,
    real scalar rcompat,
    real colvector eta_out,
    real colvector lambda_out
)
{
    // Generate or copy individual fixed effects
    if (rows(eta) == 0) {
        if (rcompat == 1) {
            // R compat mode: zero fixed effects
            eta_out = J(n, 1, 0)
        }
        else {
            // Paper mode: N(0,1) random fixed effects
            eta_out = rnormal(n, 1, 0, 1)
        }
    }
    else {
        eta_out = eta
    }
    
    // Generate or copy time fixed effects
    if (rows(lambda) == 0) {
        if (rcompat == 1) {
            // R compat mode: zero fixed effects
            lambda_out = J(n_periods, 1, 0)
        }
        else {
            // Paper mode: N(0,1) random fixed effects
            lambda_out = rnormal(n_periods, 1, 0, 1)
        }
    }
    else {
        lambda_out = lambda
    }
}


// ============================================================================
// FUNCTION: _eqt_sim_dip()
// Apply Ashenfelter's dip shock to error terms
//
// Arguments:
//   u         - Error term matrix (n x n_periods)
//   G         - Treatment indicator (n x 1)
//   T         - Number of pre-treatment periods
//   dip_mean  - Mean of dip shock distribution
//
// Returns:
//   Modified error term matrix with dip shock applied
//
// Algorithm:
//   For treated individuals (G_i = 1) at base period (t = T+1):
//   u_tilde[i, T+1] = u[i, T+1] + V_i, where V_i ~ N(dip_mean, 1)
// ============================================================================
real matrix _eqt_sim_dip(
    real matrix u,
    real colvector G,
    real scalar T,
    real scalar dip_mean
)
{
    real scalar n, base_period, i
    real colvector V
    real matrix u_tilde
    
    n = rows(u)
    base_period = T + 1
    
    // Copy error matrix
    u_tilde = u
    
    // Generate dip shocks for all individuals
    V = rnormal(n, 1, dip_mean, 1)
    
    // Apply dip shock only to treated individuals at base period
    for (i = 1; i <= n; i++) {
        if (G[i] == 1) {
            u_tilde[i, base_period] = u_tilde[i, base_period] + V[i]
        }
    }
    
    return(u_tilde)
}


// ============================================================================
// FUNCTION: _eqt_sim_trend()
// Apply linear trend effect to error terms
//
// Arguments:
//   u     - Error term matrix (n x n_periods)
//   G     - Treatment indicator (n x 1)
//   psi   - Trend coefficient
//
// Returns:
//   Modified error term matrix with linear trend applied
//
// Algorithm:
//   For treated individuals (G_i = 1) at all periods:
//   u_tilde[i, t] = u[i, t] + psi * t * G_i
// ============================================================================
real matrix _eqt_sim_trend(
    real matrix u,
    real colvector G,
    real scalar psi
)
{
    real scalar n, n_periods, i, t
    real matrix u_tilde
    
    n = rows(u)
    n_periods = cols(u)
    
    // Copy error matrix
    u_tilde = u
    
    // Apply linear trend to treated individuals
    for (i = 1; i <= n; i++) {
        if (G[i] == 1) {
            for (t = 1; t <= n_periods; t++) {
                u_tilde[i, t] = u_tilde[i, t] + psi * t
            }
        }
    }
    
    return(u_tilde)
}



// ============================================================================
// FUNCTION: _eqt_sim_paneldata()
// Generate panel data following the two-way fixed effects DGP for simulations
//
// Arguments:
//   n       - Number of individuals
//   T       - Number of pre-treatment periods (generates T+2 total periods)
//   beta    - Placebo coefficients (1 x T row vector)
//   piatt   - Average treatment effect on treated
//   eta     - Individual fixed effects (n x 1) or empty for random
//   lambda  - Time fixed effects ((T+2) x 1) or empty for random
//   phi     - AR coefficients (1 x p row vector)
//   sd      - Base error standard deviation
//   het     - Heteroskedasticity parameter (treatment group multiplier)
//   dip     - Dip shock mean (. for no dip)
//   trend   - Trend coefficient (. for no trend)
//   burnins - AR process burn-in periods
//   rcompat - R compatibility mode (0 = paper, 1 = R compat)
//
// Returns:
//   Matrix of dimension (n * (T+2)) x 4 with columns: (id, period, Y, G)
//
// DGP Formula (Paper Equation 6.1):
//   Y_{i,t} = alpha_i + lambda_t + sum_{l=1}^{T} beta_l * G_i * D_l(t) 
//             + pi_ATT * G_i * D_{T+2}(t) + u_{i,t}
//
// where:
//   - D_l(t) = 1 if t == l, 0 otherwise
//   - u_{i,t} follows AR(p) process with heteroskedasticity
//
// Heteroskedasticity (Paper Section 6):
//   Error terms are drawn from AR(p) process with standard deviation (1 + het*G_i)
//   - Control group (G_i=0): sigma_i = sd
//   - Treatment group (G_i=1): sigma_i = sd * (1 + het)
// ============================================================================
real matrix _eqt_sim_paneldata(
    real scalar n,
    real scalar T,
    real rowvector beta,
    real scalar piatt,
    real colvector eta,
    real colvector lambda,
    real rowvector phi,
    real scalar sd,
    real scalar het,
    real scalar dip,
    real scalar trend,
    real scalar burnins,
    real scalar rcompat
)
{
    real scalar n_periods, N_obs, i, t, obs, sigma_i
    real colvector G, eta_gen, lambda_gen
    real matrix u, data
    real scalar Y_it, alpha_i, lambda_t, placebo_effect, att_effect
    
    // Total periods = T pre-treatment + 1 base + 1 post = T + 2
    n_periods = T + 2
    N_obs = n * n_periods
    
    // Step 1: Generate treatment assignment
    G = _eqt_sim_treatment(n, rcompat)
    
    // Step 2: Generate fixed effects
    eta_gen = J(n, 1, .)
    lambda_gen = J(n_periods, 1, .)
    _eqt_sim_fixed_effects(n, n_periods, eta, lambda, rcompat, eta_gen, lambda_gen)
    
    // Step 3: Generate error terms with AR(p) process
    // Error matrix: n x n_periods
    u = J(n, n_periods, 0)
    
    for (i = 1; i <= n; i++) {
        // Heteroskedasticity: sigma_i = sd * (1 + het * G_i)
        // Treatment group has higher variance than control when het > 0
        if (rcompat == 1) {
            // R compatibility mode: homoskedastic errors
            sigma_i = sd
        }
        else {
            // Paper mode: heteroskedastic errors per Section 6
            sigma_i = sd * (1 + het * G[i])
        }
        
        // Generate AR(p) process for this individual
        u[i, .] = _eqt_sim_ar_process(n_periods, phi, sigma_i, burnins)'
    }
    
    // Step 4: Apply PTA violation effects if specified
    if (dip != .) {
        u = _eqt_sim_dip(u, G, T, dip)
    }
    
    if (trend != .) {
        u = _eqt_sim_trend(u, G, trend)
    }
    
    // Step 5: Compute outcome variable Y using DGP formula
    // Output matrix: N_obs x 4 (id, period, Y, G)
    data = J(N_obs, 4, 0)
    
    obs = 0
    for (i = 1; i <= n; i++) {
        alpha_i = eta_gen[i]
        
        for (t = 1; t <= n_periods; t++) {
            obs = obs + 1
            
            // Individual and time fixed effects
            lambda_t = lambda_gen[t]
            
            // Placebo effect: sum_{l=1}^{T} beta_l * G_i * D_l(t)
            // D_l(t) = 1 if t == l, so only beta_t matters when t <= T
            placebo_effect = 0
            if (t <= T) {
                placebo_effect = beta[t] * G[i]
            }
            
            // ATT effect: pi_ATT * G_i * D_{T+2}(t)
            // D_{T+2}(t) = 1 if t == T+2
            att_effect = 0
            if (t == T + 2) {
                att_effect = piatt * G[i]
            }
            
            // Compute Y_{i,t}
            Y_it = alpha_i + lambda_t + placebo_effect + att_effect + u[i, t]
            
            // Store in data matrix
            data[obs, 1] = i        // id
            data[obs, 2] = t        // period
            data[obs, 3] = Y_it     // Y
            data[obs, 4] = G[i]     // G
        }
    }
    
    return(data)
}


// ============================================================================
// FUNCTION: _eqt_sim_validate_phi()
// Validate AR coefficients for stationarity using companion matrix eigenvalues
//
// Arguments:
//   phi - AR coefficients (1 x p row vector)
//
// Returns:
//   1 if valid (stationary), 0 if invalid
//
// Stationarity condition: All eigenvalues of companion matrix have modulus < 1
//
// For AR(p): X_t = phi_1*X_{t-1} + ... + phi_p*X_{t-p} + e_t
// The companion matrix is:
//   [phi_1  phi_2  ...  phi_{p-1}  phi_p]
//   [1      0      ...  0          0    ]
//   [0      1      ...  0          0    ]
//   [0      0      ...  1          0    ]
//
// This is equivalent to checking that all roots of the characteristic polynomial
// 1 - phi_1*z - phi_2*z^2 - ... - phi_p*z^p = 0 lie outside the unit circle.
// ============================================================================
real scalar _eqt_sim_validate_phi(real rowvector phi)
{
    real scalar p, max_mod
    real matrix A
    complex rowvector evals
    
    p = cols(phi)
    
    // Special case: white noise (phi = 0)
    if (p == 1 && phi[1] == 0) {
        return(1)
    }
    
    // AR(1) special case: just check |phi| < 1
    if (p == 1) {
        return(abs(phi[1]) < 1)
    }
    
    // AR(p) case: construct companion matrix and check eigenvalues
    A = J(p, p, 0)
    A[1, .] = phi
    if (p > 1) {
        A[2..p, 1..(p-1)] = I(p-1)
    }
    
    // Compute eigenvalues and check max modulus < 1
    evals = eigenvalues(A)
    max_mod = max(abs(evals))
    
    return(max_mod < 1)
}


end

