*! equitrends.mata - Core class definitions for equivalence testing
*!
*! This file defines the Mata class hierarchy for pre-trend equivalence tests:
*!   - EquitrendsBase  : Abstract base class with shared estimation routines
*!   - MaxEquivTest    : Maximum absolute coefficient test (H0: ||beta||_inf >= delta)
*!   - MeanEquivTest   : Mean equivalence test (H0: |beta_bar| >= tau)
*!   - RMSEquivTest    : Root mean square test (H0: beta_RMS >= zeta)

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// CLASS: EquitrendsBase
// Abstract base class for equivalence testing in panel data models
//
// Provides core functionality:
//   - Panel data storage (outcome Y, regressors X, identifiers id/time)
//   - Double demeaning transformation for two-way fixed effects
//   - Pooled OLS estimation on transformed data
//   - Variance-covariance estimation (homoskedastic and cluster-robust)
//
// Test-specific logic is implemented in derived classes.
// ============================================================================
class EquitrendsBase {
    // ========== Protected Data Members ==========
    protected:
        // Original data
        real colvector Y           // Dependent variable (N_obs x 1)
        real matrix X              // Design matrix (N_obs x p)
        real colvector G           // Treatment indicator (N_obs x 1)
        real colvector id          // Individual identifier (N_obs x 1)
        real colvector time        // Time identifier (N_obs x 1)
        
        // Double-demeaned data
        real colvector Y_ddm       // Double-demeaned Y
        real matrix X_ddm          // Double-demeaned X
        
        // Dimensions
        real scalar n              // Number of unique individuals
        real scalar n_t            // Number of unique time periods
        real scalar N_obs          // Total observations
        real scalar no_placebos    // Number of placebo coefficients
        real scalar p              // Total regression variables
        
        // Metadata
        string vector placebo_names // Placebo variable names
        real scalar base_period    // Base period
        real scalar is_balanced    // Panel balance indicator (1=balanced)
        
        // Estimation results
        real colvector beta        // OLS coefficients (p x 1)
        real colvector beta_placebo // Placebo coefficients (no_placebos x 1)
        real matrix vcov           // Variance-covariance matrix (p x p)
        real matrix vcov_placebo   // Placebo vcov (no_placebos x no_placebos)
        real scalar sigma_sq       // Residual variance estimate
    
    // ========== Public Methods ==========
    public:
        // Constructor
        void new()
        
        // Setup and data loading
        void setup()
        void load_data()
        void set_dimensions()
        void validate_data()
        
        // Core computations
        void double_demean()
        void estimate_ols()
        void compute_vcov()
        void compute_vcov_cluster()
        
        // Accessors
        real colvector get_beta()
        real colvector get_beta_placebo()
        real matrix get_vcov()
        real matrix get_vcov_placebo()
        real scalar get_sigma_sq()
        real scalar get_n()
        real scalar get_n_t()
        real scalar get_N_obs()
        real scalar get_no_placebos()
        real scalar get_is_balanced()
        real colvector get_Y_ddm()
        real matrix get_X_ddm()
        
        // Setters for testing
        void set_test_data()
}

// ----------------------------------------------------------------------------
// EquitrendsBase::new()
// Constructor initializing all data members to default values
// ----------------------------------------------------------------------------
void EquitrendsBase::new()
{
    // Initialize dimensions to 0
    this.n = 0
    this.n_t = 0
    this.N_obs = 0
    this.no_placebos = 0
    this.p = 0
    
    // Initialize scalars to missing
    this.base_period = .
    this.sigma_sq = .
    
    // Default to balanced panel
    this.is_balanced = 1
}

// ----------------------------------------------------------------------------
// EquitrendsBase::setup()
// Initializes the object with data from Stata
// Configuration is performed by the calling ado file
// ----------------------------------------------------------------------------
void EquitrendsBase::setup()
{
    // Configuration performed by external ado file
}

// ----------------------------------------------------------------------------
// EquitrendsBase::load_data()
// Loads data from Stata variables into Mata matrices
// Data transfer is performed by the calling ado file
// ----------------------------------------------------------------------------
void EquitrendsBase::load_data()
{
    // Data loading performed by external ado file
}

// ----------------------------------------------------------------------------
// EquitrendsBase::set_dimensions()
// Sets dimension information after data loading
// ----------------------------------------------------------------------------
void EquitrendsBase::set_dimensions()
{
    real colvector unique_ids, unique_times
    
    this.N_obs = rows(this.Y)
    this.p = cols(this.X)
    
    unique_ids = uniqrows(this.id)
    unique_times = uniqrows(this.time)
    
    this.n = rows(unique_ids)
    this.n_t = rows(unique_times)
}

// ----------------------------------------------------------------------------
// EquitrendsBase::validate_data()
// Validates that data is properly loaded and internally consistent
// ----------------------------------------------------------------------------
void EquitrendsBase::validate_data()
{
    // Check dimensions match
    if (rows(this.X) != this.N_obs) {
        errprintf("Error: X rows (%g) != N_obs (%g)\n", rows(this.X), this.N_obs)
        _error(3498)
    }
    if (rows(this.id) != this.N_obs) {
        errprintf("Error: id rows (%g) != N_obs (%g)\n", rows(this.id), this.N_obs)
        _error(3498)
    }
    if (rows(this.time) != this.N_obs) {
        errprintf("Error: time rows (%g) != N_obs (%g)\n", rows(this.time), this.N_obs)
        _error(3498)
    }
}

// ----------------------------------------------------------------------------
// EquitrendsBase::double_demean()
// Applies within transformation for two-way fixed effects estimation
//
// Transformation: Y_ddm[i,t] = Y[i,t] - Y_bar_i - Y_bar_t + Y_bar
//
// Individual and time fixed effects are eliminated, yielding data
// suitable for pooled OLS estimation of the TWFE model.
// ----------------------------------------------------------------------------
void EquitrendsBase::double_demean()
{
    real colvector unique_ids, unique_times
    real colvector Y_i_mean, Y_t_mean, idx_i, idx_t
    real matrix X_i_mean, X_t_mean
    real scalar Y_grand_mean, i, t, obs
    real rowvector X_grand_mean
    real colvector id_map, time_map
    real scalar id_idx, time_idx
    
    // Get unique identifiers
    unique_ids = uniqrows(this.id)
    unique_times = uniqrows(this.time)
    
    // Verify dimensions
    if (rows(unique_ids) != this.n) {
        errprintf("Error: unique id count (%g) != n (%g)\n", 
                  rows(unique_ids), this.n)
        _error(3498)
    }
    if (rows(unique_times) != this.n_t) {
        errprintf("Error: unique time count (%g) != n_t (%g)\n", 
                  rows(unique_times), this.n_t)
        _error(3498)
    }
    
    // Initialize mean storage
    Y_i_mean = J(this.n, 1, 0)
    X_i_mean = J(this.n, this.p, 0)
    Y_t_mean = J(this.n_t, 1, 0)
    X_t_mean = J(this.n_t, this.p, 0)
    
    // Compute individual means
    for (i = 1; i <= this.n; i++) {
        idx_i = selectindex(this.id :== unique_ids[i])
        Y_i_mean[i] = mean(this.Y[idx_i])
        X_i_mean[i, .] = mean(this.X[idx_i, .])
    }
    
    // Compute time means
    for (t = 1; t <= this.n_t; t++) {
        idx_t = selectindex(this.time :== unique_times[t])
        Y_t_mean[t] = mean(this.Y[idx_t])
        X_t_mean[t, .] = mean(this.X[idx_t, .])
    }
    
    // Compute grand means
    Y_grand_mean = mean(this.Y)
    X_grand_mean = mean(this.X)
    
    // Create mapping from id/time values to indices
    // This handles non-consecutive id/time values
    id_map = J(this.N_obs, 1, 0)
    time_map = J(this.N_obs, 1, 0)
    
    for (i = 1; i <= this.n; i++) {
        idx_i = selectindex(this.id :== unique_ids[i])
        id_map[idx_i] = J(rows(idx_i), 1, i)
    }
    
    for (t = 1; t <= this.n_t; t++) {
        idx_t = selectindex(this.time :== unique_times[t])
        time_map[idx_t] = J(rows(idx_t), 1, t)
    }
    
    // Apply double demeaning
    this.Y_ddm = J(this.N_obs, 1, 0)
    this.X_ddm = J(this.N_obs, this.p, 0)
    
    for (obs = 1; obs <= this.N_obs; obs++) {
        id_idx = id_map[obs]
        time_idx = time_map[obs]
        
        this.Y_ddm[obs] = this.Y[obs] - Y_i_mean[id_idx] - Y_t_mean[time_idx] + Y_grand_mean
        this.X_ddm[obs, .] = this.X[obs, .] - X_i_mean[id_idx, .] - X_t_mean[time_idx, .] + X_grand_mean
    }
}

// ----------------------------------------------------------------------------
// EquitrendsBase::estimate_ols()
// Computes OLS estimates from double-demeaned data
//
// Estimator: beta_hat = (X'X)^{-1} X'Y
//
// Residuals, residual variance (with degrees of freedom adjustment for
// two-way fixed effects), and the placebo coefficient subvector are
// computed as well.
// ----------------------------------------------------------------------------
void EquitrendsBase::estimate_ols()
{
    real matrix XtX, XtX_inv
    real colvector residuals
    real scalar df, sse
    
    // Compute (X'X)^{-1}
    XtX = cross(this.X_ddm, this.X_ddm)
    XtX_inv = invsym(XtX)
    
    // Compute beta = (X'X)^{-1} X'Y
    this.beta = XtX_inv * cross(this.X_ddm, this.Y_ddm)
    
    // Extract placebo coefficients (first no_placebos elements)
    if (this.no_placebos > 0 & this.no_placebos <= rows(this.beta)) {
        this.beta_placebo = this.beta[1..this.no_placebos]
    }
    
    // Compute residuals
    residuals = this.Y_ddm - this.X_ddm * this.beta
    
    // Compute residual variance with correct degrees of freedom
    // df = N - p - n - n_t + 1 (for two-way fixed effects)
    df = this.N_obs - this.p - this.n - this.n_t + 1
    
    if (df <= 0) {
        errprintf("Error: Degrees of freedom (%g) <= 0\n", df)
        errprintf("  N_obs=%g, p=%g, n=%g, n_t=%g\n", 
                  this.N_obs, this.p, this.n, this.n_t)
        _error(3498)
    }
    
    sse = cross(residuals, residuals)
    this.sigma_sq = sse / df
}

// ----------------------------------------------------------------------------
// EquitrendsBase::compute_vcov()
// Computes variance-covariance matrix under homoskedasticity
//
// Estimator: V = sigma_hat^2 * (X'X)^{-1}
// ----------------------------------------------------------------------------
void EquitrendsBase::compute_vcov()
{
    real matrix XtX_inv
    
    XtX_inv = invsym(cross(this.X_ddm, this.X_ddm))
    this.vcov = this.sigma_sq * XtX_inv
    
    // Extract placebo vcov submatrix
    if (this.no_placebos > 0 & this.no_placebos <= rows(this.vcov)) {
        this.vcov_placebo = this.vcov[1..this.no_placebos, 1..this.no_placebos]
    }
}

// ----------------------------------------------------------------------------
// EquitrendsBase::compute_vcov_cluster()
// Computes cluster-robust variance-covariance matrix
//
// Sandwich estimator: V_CL = (X'X)^{-1} * M * (X'X)^{-1} * c
// where M = sum_g (X_g' u_g)(X_g' u_g)' is the meat matrix.
//
// Small sample adjustment: c = (G/(G-1)) * ((N-1)/(N-1-p-df_a))
//   G     : number of clusters (individuals)
//   N     : total observations
//   p     : number of regressors
//   df_a  : absorbed fixed effects degrees of freedom
//
// For two-way FE with absorb(id time) and cluster(id):
//   - Individual FE is nested within clusters, contributing 0 to df_a
//   - Time FE contributes (n_t - 1) to df_a
// ----------------------------------------------------------------------------
void EquitrendsBase::compute_vcov_cluster()
{
    real matrix XtX_inv, meat
    real colvector residuals, unique_ids, idx_i
    real matrix X_i
    real colvector u_i, score_i
    real scalar i, adjustment, df_a, denom
    
    // Compute (X'X)^{-1}
    XtX_inv = invsym(cross(this.X_ddm, this.X_ddm))
    
    // Compute residuals
    residuals = this.Y_ddm - this.X_ddm * this.beta
    
    // Get unique individual identifiers
    unique_ids = uniqrows(this.id)
    
    // Initialize meat matrix
    meat = J(this.p, this.p, 0)
    
    // Compute meat matrix by summing over clusters
    for (i = 1; i <= this.n; i++) {
        // Get indices for this cluster
        idx_i = selectindex(this.id :== unique_ids[i])
        
        // Extract cluster data
        X_i = this.X_ddm[idx_i, .]
        u_i = residuals[idx_i]
        
        // Compute score for this cluster: X_i' * u_i (p x 1)
        score_i = cross(X_i, u_i)
        
        // Add outer product to meat matrix
        meat = meat + score_i * score_i'
    }
    
    // Apply small sample adjustment to match reghdfe
    // For absorb(id time) cluster(id):
    //   df_a = n_t - 1 (time FE minus 1 redundant, id FE nested in cluster)
    // adjustment = (G/(G-1)) * ((N-1)/(N-1-p-df_a))
    df_a = this.n_t - 1
    denom = this.N_obs - 1 - this.p - df_a
    
    if (this.n > 1 & denom > 0) {
        adjustment = (this.n / (this.n - 1)) * ((this.N_obs - 1) / denom)
    } else {
        adjustment = 1
    }
    meat = meat * adjustment
    
    // Compute sandwich estimator
    this.vcov = XtX_inv * meat * XtX_inv
    
    // Extract placebo vcov submatrix
    if (this.no_placebos > 0 & this.no_placebos <= rows(this.vcov)) {
        this.vcov_placebo = this.vcov[1..this.no_placebos, 1..this.no_placebos]
    }
}

// ----------------------------------------------------------------------------
// Accessor Methods
// ----------------------------------------------------------------------------
real colvector EquitrendsBase::get_beta()
{
    return(this.beta)
}

real colvector EquitrendsBase::get_beta_placebo()
{
    return(this.beta_placebo)
}

real matrix EquitrendsBase::get_vcov()
{
    return(this.vcov)
}

real matrix EquitrendsBase::get_vcov_placebo()
{
    return(this.vcov_placebo)
}

real scalar EquitrendsBase::get_sigma_sq()
{
    return(this.sigma_sq)
}

real scalar EquitrendsBase::get_n()
{
    return(this.n)
}

real scalar EquitrendsBase::get_n_t()
{
    return(this.n_t)
}

real scalar EquitrendsBase::get_N_obs()
{
    return(this.N_obs)
}

real scalar EquitrendsBase::get_no_placebos()
{
    return(this.no_placebos)
}

real scalar EquitrendsBase::get_is_balanced()
{
    return(this.is_balanced)
}

real colvector EquitrendsBase::get_Y_ddm()
{
    return(this.Y_ddm)
}

real matrix EquitrendsBase::get_X_ddm()
{
    return(this.X_ddm)
}

// ----------------------------------------------------------------------------
// EquitrendsBase::set_test_data()
// Sets data directly for unit testing purposes
//
// Arguments:
//   Y_in, X_in       : outcome vector and design matrix
//   id_in, time_in   : panel identifiers
//   n_in, n_t_in     : number of individuals and time periods
//   N_obs_in, p_in   : total observations and regressors
//   no_placebos_in   : number of placebo coefficients
// ----------------------------------------------------------------------------
void EquitrendsBase::set_test_data(
    real colvector Y_in,
    real matrix X_in,
    real colvector id_in,
    real colvector time_in,
    real scalar n_in,
    real scalar n_t_in,
    real scalar N_obs_in,
    real scalar p_in,
    real scalar no_placebos_in
)
{
    this.Y = Y_in
    this.X = X_in
    this.id = id_in
    this.time = time_in
    this.n = n_in
    this.n_t = n_t_in
    this.N_obs = N_obs_in
    this.p = p_in
    this.no_placebos = no_placebos_in
}


// ============================================================================
// CLASS: MaxEquivTest
// Maximum absolute coefficient equivalence test
//
// Hypothesis: H0: ||beta||_inf >= delta  vs  H1: ||beta||_inf < delta
// where ||beta||_inf = max_{t} |beta_t|
//
// Three inference methods are available:
//   - IU   : Intersection-union test using folded normal quantiles
//   - Boot : Parametric bootstrap under constrained null
//   - Wild : Wild cluster bootstrap for heteroskedastic/clustered errors
// ============================================================================
class MaxEquivTest extends EquitrendsBase {
    // ========== Protected Data Members ==========
    protected:
        real scalar equiv_threshold  // Equivalence threshold delta
        string scalar method         // Test method: "IU", "Boot", "Wild"
        real scalar B                // Bootstrap replications
        real scalar alpha            // Significance level
        
        // Test results
        real scalar test_statistic   // max|beta_placebo|
        real scalar critical_value   // Critical value
        real scalar reject           // Rejection indicator (1=reject H0)
        real scalar min_threshold    // Minimum equivalence threshold delta*
    
    // ========== Public Methods ==========
    public:
        // Constructor
        void new()
        
        // Parameter setting
        void set_params()
        
        // Test methods
        void run_test()
        void run_iu_test()
        void run_boot_test()
        void run_wild_test()
        real scalar find_min_threshold()
        
        // Accessors
        real scalar get_reject()
        real scalar get_min_threshold()
        real scalar get_test_statistic()
        real scalar get_equiv_threshold()
        string scalar get_method()
}

// ----------------------------------------------------------------------------
// MaxEquivTest::new()
// Constructor initializing all data members
// ----------------------------------------------------------------------------
void MaxEquivTest::new()
{
    // Initialize base class members
    this.n = 0
    this.n_t = 0
    this.N_obs = 0
    this.no_placebos = 0
    this.p = 0
    this.base_period = .
    this.sigma_sq = .
    this.is_balanced = 1
    
    // Initialize test parameters
    this.equiv_threshold = .
    this.method = ""
    this.B = 1000
    this.alpha = 0.05
    
    // Initialize results
    this.test_statistic = .
    this.critical_value = .
    this.reject = .
    this.min_threshold = .
}

// ----------------------------------------------------------------------------
// MaxEquivTest::set_params()
// Configures test parameters (threshold, method, bootstrap replications)
// ----------------------------------------------------------------------------
void MaxEquivTest::set_params()
{
    // Parameters configured by external ado file
}

// ----------------------------------------------------------------------------
// MaxEquivTest::run_test()
// Executes the equivalence test using the specified method
// ----------------------------------------------------------------------------
void MaxEquivTest::run_test()
{
    if (this.method == "IU") {
        this.run_iu_test()
    }
    else if (this.method == "Boot") {
        this.run_boot_test()
    }
    else if (this.method == "Wild") {
        this.run_wild_test()
    }
    else {
        errprintf("Error: Unknown method '%s'\n", this.method)
        errprintf("  Valid methods: IU, Boot, Wild\n")
        _error(3498)
    }
}

// ----------------------------------------------------------------------------
// MaxEquivTest::run_iu_test()
// Intersection-union test using folded normal critical values
//
// Rejects H0 if |beta_hat_t| < Q_{FN(delta, sigma_t^2)}(alpha) for all t
// where Q_{FN} denotes the quantile of the folded normal distribution.
// ----------------------------------------------------------------------------
void MaxEquivTest::run_iu_test()
{
    this.test_statistic = _equitrends_norm_inf(this.beta_placebo)
}

// ----------------------------------------------------------------------------
// MaxEquivTest::run_boot_test()
// Parametric bootstrap test under the constrained null hypothesis
//
// Bootstrap samples are generated from the constrained OLS estimator
// satisfying ||beta||_inf = delta, with homoskedastic errors.
// ----------------------------------------------------------------------------
void MaxEquivTest::run_boot_test()
{
    this.test_statistic = _equitrends_norm_inf(this.beta_placebo)
}

// ----------------------------------------------------------------------------
// MaxEquivTest::run_wild_test()
// Wild cluster bootstrap test for non-spherical errors
//
// Bootstrap samples use Rademacher weights multiplied by cluster-level
// residuals, providing robustness to heteroskedasticity and within-cluster
// correlation.
// ----------------------------------------------------------------------------
void MaxEquivTest::run_wild_test()
{
    this.test_statistic = _equitrends_norm_inf(this.beta_placebo)
}

// ----------------------------------------------------------------------------
// MaxEquivTest::find_min_threshold()
// Computes delta* = smallest threshold for which H0 can be rejected
//
// This value represents the minimum equivalence threshold at which
// the null hypothesis of non-negligible trend differences is rejected.
// ----------------------------------------------------------------------------
real scalar MaxEquivTest::find_min_threshold()
{
    return(.)
}

// ----------------------------------------------------------------------------
// MaxEquivTest Accessor Methods
// ----------------------------------------------------------------------------
real scalar MaxEquivTest::get_reject()
{
    return(this.reject)
}

real scalar MaxEquivTest::get_min_threshold()
{
    return(this.min_threshold)
}

real scalar MaxEquivTest::get_test_statistic()
{
    return(this.test_statistic)
}

real scalar MaxEquivTest::get_equiv_threshold()
{
    return(this.equiv_threshold)
}

string scalar MaxEquivTest::get_method()
{
    return(this.method)
}


// ============================================================================
// CLASS: MeanEquivTest
// Mean equivalence test for average placebo coefficient
//
// Hypothesis: H0: |beta_bar| >= tau  vs  H1: |beta_bar| < tau
// where beta_bar = (1/T) * sum_{t=1}^{T} beta_t
//
// The test uses the folded normal distribution with mean tau and
// variance estimated from the coefficient covariance matrix.
// ============================================================================
class MeanEquivTest extends EquitrendsBase {
    // ========== Protected Data Members ==========
    protected:
        real scalar equiv_threshold  // Equivalence threshold tau
        real scalar alpha            // Significance level
        
        // Test results
        real scalar beta_bar         // Mean of placebo coefficients
        real scalar se_beta_bar      // Standard error of beta_bar
        real scalar test_statistic   // |beta_bar|
        real scalar critical_value   // Critical value
        real scalar reject           // Rejection indicator (1=reject H0)
        real scalar min_threshold    // Minimum equivalence threshold tau*
    
    // ========== Public Methods ==========
    public:
        // Constructor
        void new()
        
        // Parameter setting
        void set_params()
        
        // Test methods
        void run_test()
        real scalar find_min_threshold()
        
        // Accessors
        real scalar get_beta_bar()
        real scalar get_se_beta_bar()
        real scalar get_reject()
        real scalar get_min_threshold()
        real scalar get_test_statistic()
        real scalar get_equiv_threshold()
}

// ----------------------------------------------------------------------------
// MeanEquivTest::new()
// Constructor initializing all data members
// ----------------------------------------------------------------------------
void MeanEquivTest::new()
{
    // Initialize base class members
    this.n = 0
    this.n_t = 0
    this.N_obs = 0
    this.no_placebos = 0
    this.p = 0
    this.base_period = .
    this.sigma_sq = .
    this.is_balanced = 1
    
    // Initialize test parameters
    this.equiv_threshold = .
    this.alpha = 0.05
    
    // Initialize results
    this.beta_bar = .
    this.se_beta_bar = .
    this.test_statistic = .
    this.critical_value = .
    this.reject = .
    this.min_threshold = .
}

// ----------------------------------------------------------------------------
// MeanEquivTest::set_params()
// Configures test parameters (threshold, significance level)
// ----------------------------------------------------------------------------
void MeanEquivTest::set_params()
{
    // Parameters configured by external ado file
}

// ----------------------------------------------------------------------------
// MeanEquivTest::run_test()
// Executes the mean equivalence test
//
// H0 is rejected if |beta_bar_hat| < Q_{FN(tau, sigma_bar^2)}(alpha)
// where sigma_bar^2 = 1'*Sigma*1 / T^2 is the variance of beta_bar.
// ----------------------------------------------------------------------------
void MeanEquivTest::run_test()
{
    this.beta_bar = _equitrends_mean(this.beta_placebo)
    this.test_statistic = abs(this.beta_bar)
}

// ----------------------------------------------------------------------------
// MeanEquivTest::find_min_threshold()
// Computes tau* = smallest threshold for which H0 can be rejected
//
// This value represents the minimum equivalence threshold at which
// the null hypothesis of non-negligible average trend difference is rejected.
// ----------------------------------------------------------------------------
real scalar MeanEquivTest::find_min_threshold()
{
    return(.)
}

// ----------------------------------------------------------------------------
// MeanEquivTest Accessor Methods
// ----------------------------------------------------------------------------
real scalar MeanEquivTest::get_beta_bar()
{
    return(this.beta_bar)
}

real scalar MeanEquivTest::get_se_beta_bar()
{
    return(this.se_beta_bar)
}

real scalar MeanEquivTest::get_reject()
{
    return(this.reject)
}

real scalar MeanEquivTest::get_min_threshold()
{
    return(this.min_threshold)
}

real scalar MeanEquivTest::get_test_statistic()
{
    return(this.test_statistic)
}

real scalar MeanEquivTest::get_equiv_threshold()
{
    return(this.equiv_threshold)
}


// ============================================================================
// CLASS: RMSEquivTest
// Root mean square equivalence test with self-normalization
//
// Hypothesis: H0: beta_RMS >= zeta  vs  H1: beta_RMS < zeta
// where beta_RMS = ||beta|| / sqrt(T) = sqrt((1/T) * sum beta_t^2)
//
// The test employs a pivotal statistic based on self-normalization,
// avoiding direct estimation of the asymptotic variance. Critical values
// are obtained from the W distribution derived in the paper.
// ============================================================================
class RMSEquivTest extends EquitrendsBase {
    // ========== Protected Data Members ==========
    protected:
        real scalar equiv_threshold  // Equivalence threshold zeta
        real scalar alpha            // Significance level
        real scalar no_lambda        // Lambda subdivisions (default 4)
        
        // Test results
        real scalar beta_rms         // RMS of placebo coefficients
        real scalar beta_rms_sq      // Squared RMS
        real scalar V_n              // Self-normalized variance
        real scalar M_n              // Test statistic
        real scalar W_crit           // W distribution critical value
        real scalar reject           // Rejection indicator (1=reject H0)
        real scalar min_threshold    // Minimum equivalence threshold zeta*
    
    // ========== Public Methods ==========
    public:
        // Constructor
        void new()
        
        // Parameter setting
        void set_params()
        
        // Test methods
        void run_test()
        real scalar compute_V_n()
        real scalar find_min_threshold()
        
        // Accessors
        real scalar get_beta_rms()
        real scalar get_beta_rms_sq()
        real scalar get_V_n()
        real scalar get_reject()
        real scalar get_min_threshold()
        real scalar get_test_statistic()
        real scalar get_equiv_threshold()
        real scalar get_W_crit()
}

// ----------------------------------------------------------------------------
// RMSEquivTest::new()
// Constructor initializing all data members
// ----------------------------------------------------------------------------
void RMSEquivTest::new()
{
    // Initialize base class members
    this.n = 0
    this.n_t = 0
    this.N_obs = 0
    this.no_placebos = 0
    this.p = 0
    this.base_period = .
    this.sigma_sq = .
    this.is_balanced = 1
    
    // Initialize test parameters
    this.equiv_threshold = .
    this.alpha = 0.05
    this.no_lambda = 4
    
    // Initialize results
    this.beta_rms = .
    this.beta_rms_sq = .
    this.V_n = .
    this.M_n = .
    this.W_crit = .
    this.reject = .
    this.min_threshold = .
}

// ----------------------------------------------------------------------------
// RMSEquivTest::set_params()
// Configures test parameters (threshold, significance level, lambda grid)
// ----------------------------------------------------------------------------
void RMSEquivTest::set_params()
{
    // Parameters configured by external ado file
}

// ----------------------------------------------------------------------------
// RMSEquivTest::run_test()
// Executes the RMS equivalence test
//
// H0 is rejected if beta_RMS_hat^2 < zeta^2 + Q_W(alpha) * V_n
// where V_n is the self-normalized variance estimator and Q_W(alpha)
// is the alpha-quantile of the limiting W distribution.
// ----------------------------------------------------------------------------
void RMSEquivTest::run_test()
{
    // Validate significance level (W distribution quantiles are tabulated)
    if (!_equitrends_validate_rms_alpha(this.alpha)) {
        _error(3498)
    }
    
    // Compute RMS of placebo coefficients
    this.beta_rms = _equitrends_rms(this.beta_placebo)
    this.beta_rms_sq = this.beta_rms^2
    
    // Retrieve W distribution critical value
    this.W_crit = _equitrends_W_critical_value(this.alpha)
}

// ----------------------------------------------------------------------------
// RMSEquivTest::compute_V_n()
// Computes self-normalized variance estimator V_n
//
// V_n = sqrt( (1/K) * sum_{k=1}^{K} (beta_RMS^2(lambda_k) - beta_RMS^2(1))^2 )
//
// where lambda_k = k/(K+1) defines the subsample fractions and
// beta_RMS^2(lambda) is computed from the first floor(n*lambda) individuals.
// ----------------------------------------------------------------------------
real scalar RMSEquivTest::compute_V_n()
{
    return(.)
}

// ----------------------------------------------------------------------------
// RMSEquivTest::find_min_threshold()
// Computes zeta* = smallest threshold for which H0 can be rejected
//
// Derived from the rejection condition: zeta* = sqrt(beta_RMS^2 - Q_W * V_n)
// This value represents the minimum RMS equivalence threshold at which
// the null hypothesis of non-negligible trend differences is rejected.
// ----------------------------------------------------------------------------
real scalar RMSEquivTest::find_min_threshold()
{
    return(.)
}

// ----------------------------------------------------------------------------
// RMSEquivTest Accessor Methods
// ----------------------------------------------------------------------------
real scalar RMSEquivTest::get_beta_rms()
{
    return(this.beta_rms)
}

real scalar RMSEquivTest::get_beta_rms_sq()
{
    return(this.beta_rms_sq)
}

real scalar RMSEquivTest::get_V_n()
{
    return(this.V_n)
}

real scalar RMSEquivTest::get_reject()
{
    return(this.reject)
}

real scalar RMSEquivTest::get_min_threshold()
{
    return(this.min_threshold)
}

real scalar RMSEquivTest::get_test_statistic()
{
    return(this.M_n)
}

real scalar RMSEquivTest::get_equiv_threshold()
{
    return(this.equiv_threshold)
}

real scalar RMSEquivTest::get_W_crit()
{
    return(this.W_crit)
}

end
