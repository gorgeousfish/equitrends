*! equitrends_dataproc.mata - Data processing functions for EQUITRENDS
*!
*! Implements the EquitrendsDataProcessor class for panel data preparation.
*! Handles data binding, missing values, pretreatment period filtering,
*! base period specification, placebo variable construction, ID recoding,
*! and panel balance verification.

version 16.0

mata:
mata set matastrict on

// ============================================================================
// _equitrends_uniq_preserve_order()
// Returns unique values from a vector while preserving order of first appearance
//
// Uses asarray() hash set for O(N) complexity. This function is defined before
// the class that uses it due to Mata's parsing requirements.
//
// Arguments:
//   x : real colvector - Input vector
//
// Returns:
//   real colvector - Unique values in first-appearance order
// ============================================================================

real colvector _equitrends_uniq_preserve_order(real colvector x)
{
    real colvector result
    real scalar i, n, count
    transmorphic A
    
    n = rows(x)
    
    if (n == 0) {
        return(J(0, 1, .))
    }
    
    // Hash set for O(1) lookup; exact equality is used since IDs are integers
    A = asarray_create("real")
    
    // Pre-allocate result to maximum possible size, then trim
    result = J(n, 1, .)
    count = 0
    
    for (i = 1; i <= n; i++) {
        if (asarray_contains(A, x[i]) == 0) {
            // Value not yet seen; add to result and mark as seen
            count = count + 1
            result[count] = x[i]
            asarray(A, x[i], 1)
        }
    }
    
    // Trim result to actual size
    if (count > 0) {
        result = result[1::count]
    }
    else {
        result = J(0, 1, .)
    }
    
    return(result)
}

// ============================================================================
// EquitrendsDataProcessor Class Definition
// ============================================================================

class EquitrendsDataProcessor {
    // ========================================================================
    // Input data (protected)
    // ========================================================================
    protected:
        real colvector Y
        real colvector ID
        real colvector G
        real colvector period
        real matrix X
        real colvector cluster
        real colvector pretreatment_period
        real scalar base_period
        real scalar base_period_set
        string vector input_orig_x_names      // Original covariate names for output
    
    // ========================================================================
    // Processed data (protected)
    // ========================================================================
    protected:
        real colvector Y_proc
        real colvector ID_proc
        real colvector G_proc
        real colvector period_proc
        real matrix X_proc
        real colvector cluster_proc
        real matrix placebo_vars
        string vector placebo_names
        string vector orig_X_names
    
    // ========================================================================
    // Metadata (protected)
    // ========================================================================
    protected:
        real scalar n
        real scalar n_t
        real scalar N_obs
        real scalar no_placebos
        real scalar is_balanced
        real rowvector no_periods_info
        real scalar na_omitted
    
    // ========================================================================
    // Public methods
    // ========================================================================
    public:
        void new()
        void load_data()
        void process()
        void bind_and_rename()
        real scalar handle_na()
        void filter_pretreatment()
        void set_base_period()
        void create_placebos()
        void recode_id()
        void check_balance()
        void set_orig_x_names()
        real colvector get_Y()
        real colvector get_ID()
        real colvector get_G()
        real colvector get_period()
        real matrix get_X()
        real colvector get_cluster()
        real matrix get_placebo_vars()
        string vector get_placebo_names()
        string vector get_orig_X_names()
        real scalar get_base_period()
        real scalar get_n()
        real scalar get_n_t()
        real scalar get_N_obs()
        real scalar get_no_placebos()
        real scalar get_is_balanced()
        real rowvector get_no_periods_info()
        real scalar get_na_omitted()
}

// ============================================================================
// Constructor
// ============================================================================

void EquitrendsDataProcessor::new()
{
    this.n = 0
    this.n_t = 0
    this.N_obs = 0
    this.no_placebos = 0
    this.is_balanced = 0
    this.na_omitted = 0
    this.base_period = .
    this.base_period_set = 0
    
    this.Y = J(0, 1, .)
    this.ID = J(0, 1, .)
    this.G = J(0, 1, .)
    this.period = J(0, 1, .)
    this.X = J(0, 0, .)
    this.cluster = J(0, 1, .)
    this.pretreatment_period = J(0, 1, .)
    this.input_orig_x_names = J(1, 0, "")
    
    this.Y_proc = J(0, 1, .)
    this.ID_proc = J(0, 1, .)
    this.G_proc = J(0, 1, .)
    this.period_proc = J(0, 1, .)
    this.X_proc = J(0, 0, .)
    this.cluster_proc = J(0, 1, .)
    this.placebo_vars = J(0, 0, .)
    
    this.placebo_names = J(1, 0, "")
    this.orig_X_names = J(1, 0, "")
    this.no_periods_info = J(1, 0, .)
}

// ============================================================================
// load_data() - Load input data into the processor
// ============================================================================

void EquitrendsDataProcessor::load_data(
    real colvector Y,
    real colvector ID,
    real colvector G,
    real colvector period,
    | real matrix X,
    real colvector cluster,
    real colvector pretreatment_period,
    real scalar base_period
)
{
    this.Y = Y
    this.ID = ID
    this.G = G
    this.period = period
    
    if (args() >= 5 & cols(X) > 0) {
        this.X = X
    }
    
    if (args() >= 6 & rows(cluster) > 0) {
        this.cluster = cluster
    }
    
    if (args() >= 7 & rows(pretreatment_period) > 0) {
        this.pretreatment_period = pretreatment_period
    }
    
    if (args() >= 8 & !missing(base_period)) {
        this.base_period = base_period
        this.base_period_set = 1
    }
}

// ============================================================================
// process() - Main processing pipeline
// ============================================================================

void EquitrendsDataProcessor::process()
{
    this.bind_and_rename()
    this.na_omitted = this.handle_na()
    this.filter_pretreatment()
    this.set_base_period()
    this.create_placebos()
    
    // Require at least one placebo (i.e., at least 2 unique periods after filtering)
    if (this.no_placebos < 1) {
        errprintf("Error: After NA removal and pretreatment filtering, ")
        errprintf("data must contain at least 2 unique periods.\n")
        errprintf("Current data has only 1 period (the base period).\n")
        exit(198)
    }
    
    this.recode_id()
    this.check_balance()
}

// ============================================================================
// bind_and_rename() - Data binding and variable renaming
//
// Copies input data to processing variables. If original covariate names are
// provided, they are preserved; otherwise, generic names (X_1, X_2, ...) are
// generated.
// ============================================================================

void EquitrendsDataProcessor::bind_and_rename()
{
    real scalar k
    
    if (cols(this.X) > 0) {
        if (cols(this.input_orig_x_names) == cols(this.X)) {
            this.orig_X_names = this.input_orig_x_names
        }
        else {
            // Generate default names if originals not provided
            this.orig_X_names = J(1, cols(this.X), "")
            for (k = 1; k <= cols(this.X); k++) {
                this.orig_X_names[k] = "X_" + strofreal(k)
            }
        }
    }
    
    this.Y_proc = this.Y
    this.ID_proc = this.ID
    this.G_proc = this.G
    this.period_proc = this.period
    this.X_proc = this.X
    this.cluster_proc = this.cluster
}

// ============================================================================
// handle_na() - Missing value handling
// ============================================================================

real scalar EquitrendsDataProcessor::handle_na()
{
    real colvector valid_rows
    real scalar original_n, new_n, omitted
    
    original_n = rows(this.Y_proc)
    valid_rows = J(original_n, 1, 1)
    
    valid_rows = valid_rows :& !rowmissing(this.Y_proc)
    valid_rows = valid_rows :& !rowmissing(this.ID_proc)
    valid_rows = valid_rows :& !rowmissing(this.G_proc)
    valid_rows = valid_rows :& !rowmissing(this.period_proc)
    
    if (cols(this.X_proc) > 0) {
        valid_rows = valid_rows :& (rowmissing(this.X_proc) :== 0)
    }
    
    if (rows(this.cluster_proc) > 0) {
        valid_rows = valid_rows :& !rowmissing(this.cluster_proc)
    }
    
    this.Y_proc = select(this.Y_proc, valid_rows)
    this.ID_proc = select(this.ID_proc, valid_rows)
    this.G_proc = select(this.G_proc, valid_rows)
    this.period_proc = select(this.period_proc, valid_rows)
    
    if (cols(this.X_proc) > 0) {
        this.X_proc = select(this.X_proc, valid_rows)
    }
    
    if (rows(this.cluster_proc) > 0) {
        this.cluster_proc = select(this.cluster_proc, valid_rows)
    }
    
    new_n = rows(this.Y_proc)
    omitted = original_n - new_n
    
    return(omitted)
}

// ============================================================================
// filter_pretreatment() - Filter to pretreatment periods
//
// Restricts data to explicitly specified pretreatment periods for placebo
// analysis. If no pretreatment periods are specified, no filtering is applied.
// ============================================================================

void EquitrendsDataProcessor::filter_pretreatment()
{
    real colvector keep_rows
    real scalar i, j, n_pretreat
    real scalar tol
    
    if (rows(this.pretreatment_period) == 0) {
        return
    }
    
    tol = 1e-10
    n_pretreat = rows(this.pretreatment_period)
    keep_rows = J(rows(this.period_proc), 1, 0)
    
    for (i = 1; i <= rows(this.period_proc); i++) {
        for (j = 1; j <= n_pretreat; j++) {
            if (abs(this.period_proc[i] - this.pretreatment_period[j]) < tol) {
                keep_rows[i] = 1
                break
            }
        }
    }
    
    this.Y_proc = select(this.Y_proc, keep_rows)
    this.ID_proc = select(this.ID_proc, keep_rows)
    this.G_proc = select(this.G_proc, keep_rows)
    this.period_proc = select(this.period_proc, keep_rows)
    
    if (cols(this.X_proc) > 0) {
        this.X_proc = select(this.X_proc, keep_rows)
    }
    
    if (rows(this.cluster_proc) > 0) {
        this.cluster_proc = select(this.cluster_proc, keep_rows)
    }
}

// ============================================================================
// set_base_period() - Set base period (default to max if not specified)
// ============================================================================

void EquitrendsDataProcessor::set_base_period()
{
    if (this.base_period_set == 0) {
        this.base_period = max(this.period_proc)
    }
}

// ============================================================================
// create_placebos() - Construct placebo variables
//
// Creates placebo indicator variables W_{i,t,l} = G_i * D_l(t) for testing
// the parallel trends assumption. Each placebo corresponds to a pretreatment
// period relative to the base period.
// ============================================================================

void EquitrendsDataProcessor::create_placebos()
{
    real colvector unique_time, non_base_time
    real scalar l, t
    real colvector indicator
    real scalar tol
    
    tol = 1e-10
    // Preserve first-appearance order of periods for consistent output
    unique_time = _equitrends_uniq_preserve_order(this.period_proc)
    non_base_time = select(unique_time, abs(unique_time :- this.base_period) :>= tol)
    
    this.no_placebos = rows(non_base_time)
    this.placebo_vars = J(rows(this.Y_proc), this.no_placebos, 0)
    this.placebo_names = J(1, this.no_placebos, "")
    
    for (l = 1; l <= this.no_placebos; l++) {
        t = non_base_time[l]
        indicator = (abs(this.period_proc :- t) :< tol)
        this.placebo_vars[., l] = this.G_proc :* indicator
        this.placebo_names[l] = "placebo_" + strofreal(t)
    }
}

// ============================================================================
// recode_id() - Recode ID to consecutive integers (1 to N)
//
// Transforms original panel identifiers to consecutive integers for efficient
// indexing in subsequent computations. Uses hash table for O(N+n) complexity.
// ============================================================================

void EquitrendsDataProcessor::recode_id()
{
    real colvector unique_ids, new_id
    real scalar i, N_obs
    transmorphic A
    
    unique_ids = _equitrends_uniq_preserve_order(this.ID_proc)
    this.n = rows(unique_ids)
    N_obs = rows(this.ID_proc)
    new_id = J(N_obs, 1, 0)
    
    // Build hash map from original ID to consecutive index
    A = asarray_create("real")
    for (i = 1; i <= this.n; i++) {
        asarray(A, unique_ids[i], i)
    }
    
    // Apply mapping to all observations
    for (i = 1; i <= N_obs; i++) {
        new_id[i] = asarray(A, this.ID_proc[i])
    }
    
    this.ID_proc = new_id
}

// ============================================================================
// check_balance() - Check panel balance
//
// Verifies whether the panel is balanced (all individuals observed in all
// periods) and computes summary statistics about the panel structure.
// ============================================================================

void EquitrendsDataProcessor::check_balance()
{
    real colvector unique_ids, unique_periods
    real colvector first_id_periods, current_periods
    real scalar i, is_equal
    real colvector periods_count
    
    unique_ids = uniqrows(this.ID_proc)
    unique_periods = uniqrows(this.period_proc)
    
    this.n = rows(unique_ids)
    this.n_t = rows(unique_periods)
    this.N_obs = rows(this.Y_proc)
    
    first_id_periods = sort(
        select(this.period_proc, this.ID_proc :== unique_ids[1]), 1
    )
    
    this.is_balanced = 1
    periods_count = J(this.n, 1, 0)
    
    for (i = 1; i <= this.n; i++) {
        current_periods = sort(
            select(this.period_proc, this.ID_proc :== unique_ids[i]), 1
        )
        periods_count[i] = rows(current_periods)
        
        if (rows(current_periods) != rows(first_id_periods)) {
            this.is_balanced = 0
        } 
        else {
            is_equal = (sum(current_periods :== first_id_periods) == rows(first_id_periods))
            if (!is_equal) {
                this.is_balanced = 0
            }
        }
    }
    
    if (this.is_balanced) {
        this.no_periods_info = this.n_t
    } 
    else {
        this.no_periods_info = (min(periods_count), max(periods_count))
    }
}

// ============================================================================
// set_orig_x_names() - Set original covariate variable names
//
// Stores original variable names for use in output and diagnostic messages.
// Should be called before process() if custom names are desired.
// ============================================================================

void EquitrendsDataProcessor::set_orig_x_names(string vector names)
{
    this.input_orig_x_names = names
}

// ============================================================================
// Accessor methods - Processed data
// ============================================================================

real colvector EquitrendsDataProcessor::get_Y()
{
    return(this.Y_proc)
}

real colvector EquitrendsDataProcessor::get_ID()
{
    return(this.ID_proc)
}

real colvector EquitrendsDataProcessor::get_G()
{
    return(this.G_proc)
}

real colvector EquitrendsDataProcessor::get_period()
{
    return(this.period_proc)
}

real matrix EquitrendsDataProcessor::get_X()
{
    return(this.X_proc)
}

real colvector EquitrendsDataProcessor::get_cluster()
{
    return(this.cluster_proc)
}

real matrix EquitrendsDataProcessor::get_placebo_vars()
{
    return(this.placebo_vars)
}

string vector EquitrendsDataProcessor::get_placebo_names()
{
    return(this.placebo_names)
}

string vector EquitrendsDataProcessor::get_orig_X_names()
{
    return(this.orig_X_names)
}

real scalar EquitrendsDataProcessor::get_base_period()
{
    return(this.base_period)
}

// ============================================================================
// Accessor methods - Metadata
// ============================================================================

real scalar EquitrendsDataProcessor::get_n()
{
    return(this.n)
}

real scalar EquitrendsDataProcessor::get_n_t()
{
    return(this.n_t)
}

real scalar EquitrendsDataProcessor::get_N_obs()
{
    return(this.N_obs)
}

real scalar EquitrendsDataProcessor::get_no_placebos()
{
    return(this.no_placebos)
}

real scalar EquitrendsDataProcessor::get_is_balanced()
{
    return(this.is_balanced)
}

real rowvector EquitrendsDataProcessor::get_no_periods_info()
{
    return(this.no_periods_info)
}

real scalar EquitrendsDataProcessor::get_na_omitted()
{
    return(this.na_omitted)
}

end
