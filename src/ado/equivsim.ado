*! equivsim.ado - Panel Data Simulation for Equivalence Testing
*!
*! Generates panel data for Monte Carlo simulation studies of equivalence tests.
*!
*! Syntax:
*!   equivsim, n(integer) preperiods(integer) ///
*!       [beta(numlist) piatt(real 0)] ///
*!       [eta(numlist) lambda(numlist)] ///
*!       [dgp(string) phi(numlist) sd(real 1) het(real 1)] ///
*!       [dip(real) trend(real)] ///
*!       [burnins(integer 100) seed(integer)] ///
*!       [rcompat clear replace]
*!
*! Required Options:
*!   n()          - Number of individuals
*!   preperiods() - Number of pre-treatment periods T (generates T+2 total periods)
*!
*! Optional Parameters:
*!   beta()       - Placebo coefficients (length = preperiods), default: all zeros
*!   piatt()      - Average treatment effect on treated, default: 0
*!   eta()        - Individual fixed effects (length = n), default: N(0,1) random
*!   lambda()     - Time fixed effects (length = preperiods+2), default: N(0,1) random
*!   dgp()        - DGP type: spherical | ar1 | ar3, default: ar3
*!                  spherical: i.i.d. N(0, sigma^2) errors (no serial correlation)
*!                  ar1: AR(1) with phi=0.5 + heteroskedasticity
*!                  ar3: AR(3) with phi=(0.5, 0.3, 0.1) + heteroskedasticity
*!   phi()        - AR coefficients (overrides dgp default), default depends on dgp()
*!   sd()         - Base error standard deviation, default: 1
*!   het()        - Heteroskedasticity parameter, default: 1
*!                  Treatment group has sigma = sd * (1 + het), control has sigma = sd
*!                  Paper default: het=1 means treatment variance is 4x control variance
*!                  For spherical DGP, het is ignored (homoskedastic)
*!   dip()        - Ashenfelter's dip shock mean (omit for no dip)
*!   trend()      - Linear trend coefficient (omit for no trend)
*!   burnins()    - AR process burn-in periods, default: 100
*!   seed()       - Random seed for reproducibility
*!   rcompat      - Alternative error generation mode for cross-validation studies
*!   clear        - Clear existing data
*!   replace      - Replace existing data
*!
*! Output Variables:
*!   id     - Individual identifier (long, 1 to n)
*!   period - Time period identifier (int, 1 to T+2)
*!   Y      - Outcome variable (double)
*!   G      - Treatment group indicator (byte, 0 or 1)

program define equivsim, rclass
    version 16.0
    
    // Initialize Mata functions (compile on first use if needed)
    equitrends_init
    
    syntax, N(integer) PREperiods(integer) ///
        [Beta(numlist) PIatt(real 0)] ///
        [Eta(numlist) Lambda(numlist)] ///
        [DGP(string) PHI(numlist) SD(real 1) HET(real 1)] ///
        [DIP(real -999999) TREND(real -999999)] ///
        [BURNins(integer 100) SEED(integer -999999)] ///
        [Rcompat CLEAR REPLACE]
    
    // -------------------------------------------------------------------------
    // Step 1: Input Validation
    // -------------------------------------------------------------------------
    
    // Validate n
    if `n' <= 0 {
        display as error "n must be a positive integer"
        exit 198
    }
    
    // Validate preperiods
    if `preperiods' <= 0 {
        display as error "preperiods must be a positive integer"
        exit 198
    }
    
    // Validate sd
    if `sd' <= 0 {
        display as error "sd must be a positive number"
        exit 198
    }
    
    // Validate het
    if `het' < 0 {
        display as error "het must be a non-negative number"
        exit 198
    }
    
    // Validate burnins
    if `burnins' <= 0 {
        display as error "burnins must be a positive integer"
        exit 198
    }
    
    // Calculate total periods
    local n_periods = `preperiods' + 2
    local N_obs = `n' * `n_periods'
    
    // -------------------------------------------------------------------------
    // Step 2: Process beta parameter
    // -------------------------------------------------------------------------
    
    if "`beta'" == "" {
        // Default: all zeros
        local beta_vec ""
        forvalues i = 1/`preperiods' {
            local beta_vec "`beta_vec' 0"
        }
        local beta_vec = trim("`beta_vec'")
    }
    else {
        local beta_vec "`beta'"
        // Validate length
        local beta_len : word count `beta_vec'
        if `beta_len' != `preperiods' {
            display as error "beta must be a vector of length preperiods"
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Step 3: Process eta parameter (individual fixed effects)
    // -------------------------------------------------------------------------
    
    local eta_specified = 0
    if "`eta'" != "" {
        local eta_vec "`eta'"
        local eta_len : word count `eta_vec'
        if `eta_len' != `n' {
            display as error "eta must be a vector of length n"
            exit 198
        }
        local eta_specified = 1
    }
    
    // -------------------------------------------------------------------------
    // Step 4: Process lambda parameter (time fixed effects)
    // -------------------------------------------------------------------------
    
    local lambda_specified = 0
    if "`lambda'" != "" {
        local lambda_vec "`lambda'"
        local lambda_len : word count `lambda_vec'
        if `lambda_len' != `n_periods' {
            display as error "lambda must be a vector of length preperiods+2"
            exit 198
        }
        local lambda_specified = 1
    }
    
    // -------------------------------------------------------------------------
    // Step 5: Process dgp and phi parameters (AR coefficients)
    // -------------------------------------------------------------------------
    
    // Validate and process dgp option
    local dgp_type "ar3"  // Default DGP type
    if "`dgp'" != "" {
        local dgp_lower = lower("`dgp'")
        if !inlist("`dgp_lower'", "spherical", "ar1", "ar3") {
            display as error "dgp must be one of: spherical, ar1, ar3"
            exit 198
        }
        local dgp_type "`dgp_lower'"
    }
    
    // Set phi based on dgp type (unless user specified phi)
    if "`phi'" == "" {
        if "`dgp_type'" == "spherical" {
            // Spherical: white noise (phi = 0)
            local phi_vec "0"
        }
        else if "`dgp_type'" == "ar1" {
            // AR(1) with phi = 0.5
            local phi_vec "0.5"
        }
        else if "`dgp_type'" == "ar3" {
            // AR(3) with phi = (0.5, 0.3, 0.1)
            local phi_vec "0.5 0.3 0.1"
        }
        
        // Alternative mode defaults to AR(1) when dgp not specified
        if "`rcompat'" != "" & "`dgp'" == "" {
            local phi_vec "0.5"
            local dgp_type "ar1"
        }
    }
    else {
        local phi_vec "`phi'"
    }
    
    // For spherical DGP, force het = 0 (homoskedastic)
    local het_val = `het'
    if "`dgp_type'" == "spherical" {
        local het_val = 0
    }
    
    // Validate stationarity: check that all eigenvalues of the AR companion
    // matrix have modulus < 1
    mata: st_local("ar_stable", strofreal(_eqt_sim_validate_phi(strtoreal(tokens(st_local("phi_vec"))))))
    if `ar_stable' == 0 {
        display as error "AR coefficients do not satisfy stationarity condition"
        display as error "(companion matrix has eigenvalue with modulus >= 1)"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Step 6: Process dip and trend parameters
    // -------------------------------------------------------------------------
    
    local dip_specified = 0
    if `dip' != -999999 {
        local dip_specified = 1
        local dip_val = `dip'
    }
    else {
        local dip_val = .
    }
    
    local trend_specified = 0
    if `trend' != -999999 {
        local trend_specified = 1
        local trend_val = `trend'
    }
    else {
        local trend_val = .
    }
    
    // -------------------------------------------------------------------------
    // Step 7: Check existing data
    // -------------------------------------------------------------------------
    
    if "`clear'" == "" & "`replace'" == "" {
        quietly count
        if r(N) > 0 {
            display as error "no; data in memory would be lost"
            exit 4
        }
    }
    
    // -------------------------------------------------------------------------
    // Step 8: Set random seed if specified
    // -------------------------------------------------------------------------
    
    if `seed' != -999999 {
        if `seed' < 0 {
            display as error "seed must be non-negative"
            exit 198
        }
        set seed `seed'
    }
    
    // -------------------------------------------------------------------------
    // Step 9: Determine alternative mode
    // -------------------------------------------------------------------------
    
    local rcompat_flag = 0
    if "`rcompat'" != "" {
        local rcompat_flag = 1
    }
    
    // -------------------------------------------------------------------------
    // Step 10: Call Mata function to generate data
    // -------------------------------------------------------------------------
    
    // Clear existing data
    if "`clear'" != "" | "`replace'" != "" {
        clear
    }
    
    // Convert parameters to Mata and generate data
    mata: _equivsim_generate( ///
        `n', ///
        `preperiods', ///
        "`beta_vec'", ///
        `piatt', ///
        `eta_specified', ///
        "`eta'", ///
        `lambda_specified', ///
        "`lambda'", ///
        "`phi_vec'", ///
        `sd', ///
        `het_val', ///
        `dip_val', ///
        `trend_val', ///
        `burnins', ///
        `rcompat_flag' ///
    )
    
    // -------------------------------------------------------------------------
    // Step 11: Set variable types and labels
    // -------------------------------------------------------------------------
    
    // Set variable types
    recast long id
    recast int period
    recast double Y
    recast byte G
    
    // Set variable labels
    label variable id "Individual identifier"
    label variable period "Time period"
    label variable Y "Outcome variable"
    label variable G "Treatment group indicator"
    
    // Sort data
    sort id period
    
    // Calculate treatment group percentage
    quietly summarize G
    local mean_G = r(mean)
    
    // Display summary
    display as text ""
    display as text "Panel data generated:"
    display as text "  Individuals (n):     " as result `n'
    display as text "  Pre-periods (T):     " as result `preperiods'
    display as text "  Total periods:       " as result `n_periods'
    display as text "  Total observations:  " as result `N_obs'
    display as text "  Treatment group:     " as result %5.1f 100*`mean_G' "%"
    display as text "  DGP type:            " as result "`dgp_type'"
    if "`dgp_type'" != "spherical" {
        display as text "  Heteroskedasticity:  " as result "het = `het_val'"
    }
    if "`rcompat'" != "" {
        display as text "  Mode:                " as result "Alternative"
    }
    else {
        display as text "  Mode:                " as result "Paper (default)"
    }
    
    // -------------------------------------------------------------------------
    // Step 12: Store return values
    // -------------------------------------------------------------------------
    
    // Calculate treatment group counts
    quietly count if G == 1
    local n_treated = r(N)
    local n_control = `n' * `n_periods' - `n_treated'
    local n_treated_ind = `n_treated' / `n_periods'
    local n_control_ind = `n_control' / `n_periods'
    
    return scalar n = `n'
    return scalar preperiods = `preperiods'
    return scalar n_obs = `N_obs'
    return scalar n_treated = `n_treated_ind'
    return scalar n_control = `n_control_ind'
    if `seed' != -999999 {
        return scalar seed = `seed'
    }
    
end


// ============================================================================
// _equivsim_generate()
// Generates panel data for equivalence testing simulations
//
// This function creates a balanced panel dataset with individual and time
// fixed effects, AR(p) error structure, and optional treatment effects.
// The data generating process follows Dette and Schumann (2024).
//
// Arguments:
//   n             : scalar - Number of individuals
//   T             : scalar - Number of pre-treatment periods
//   beta_str      : string - Placebo coefficients (space-separated)
//   piatt         : scalar - Average treatment effect on treated
//   eta_specified : scalar - Flag: 1 if eta provided, 0 otherwise
//   eta_str       : string - Individual fixed effects (space-separated)
//   lambda_specified : scalar - Flag: 1 if lambda provided, 0 otherwise
//   lambda_str    : string - Time fixed effects (space-separated)
//   phi_str       : string - AR coefficients (space-separated)
//   sd            : scalar - Base error standard deviation
//   het           : scalar - Heteroskedasticity parameter
//   dip           : scalar - Ashenfelter's dip shock mean (. for none)
//   trend         : scalar - Linear trend coefficient (. for none)
//   burnins       : scalar - AR process burn-in periods
//   rcompat       : scalar - Alternative mode flag
//
// Returns:
//   Creates Stata variables: id, period, Y, G
// ============================================================================
mata:

void _equivsim_generate(
    real scalar n,
    real scalar T,
    string scalar beta_str,
    real scalar piatt,
    real scalar eta_specified,
    string scalar eta_str,
    real scalar lambda_specified,
    string scalar lambda_str,
    string scalar phi_str,
    real scalar sd,
    real scalar het,
    real scalar dip,
    real scalar trend,
    real scalar burnins,
    real scalar rcompat
)
{
    real rowvector beta, phi
    real colvector eta, lambda
    real matrix data
    real scalar n_periods
    
    // Parse beta string to row vector
    beta = strtoreal(tokens(beta_str))
    
    // Parse phi string to row vector
    phi = strtoreal(tokens(phi_str))
    
    // Parse eta if specified
    if (eta_specified) {
        eta = strtoreal(tokens(eta_str))'
    }
    else {
        eta = J(0, 1, .)
    }
    
    // Parse lambda if specified
    if (lambda_specified) {
        lambda = strtoreal(tokens(lambda_str))'
    }
    else {
        lambda = J(0, 1, .)
    }
    
    // Generate panel data
    data = _eqt_sim_paneldata(n, T, beta, piatt, eta, lambda, phi, sd, het, dip, trend, burnins, rcompat)
    
    // Create Stata variables
    n_periods = T + 2
    
    st_addobs(rows(data))
    
    (void) st_addvar("double", "id")
    (void) st_addvar("double", "period")
    (void) st_addvar("double", "Y")
    (void) st_addvar("double", "G")
    
    st_store(., "id", data[., 1])
    st_store(., "period", data[., 2])
    st_store(., "Y", data[., 3])
    st_store(., "G", data[., 4])
    
    // Store mean of G for display
    st_numscalar("r(mean_G)", mean(data[., 4]))
}

end

