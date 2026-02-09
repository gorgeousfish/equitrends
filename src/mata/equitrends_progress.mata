*! equitrends_progress.mata - Progress bar display for EQUITRENDS bootstrap
*!
*! Provides a professional progress bar for long-running bootstrap computations.
*!
*! In interactive mode (default):
*!   Uses carriage return char(13) for in-place line updates, showing:
*!   - Visual progress bar with filled/empty segments
*!   - Percentage complete, iteration counter, elapsed time, ETA
*!
*! In batch mode (stata -b):
*!   Uses milestone-based display (every 10%) with newlines, since
*!   carriage returns do not overwrite in log files.
*!
*! The progress bar updates adaptively to minimize display overhead:
*!   - Updates at most once per ~0.3 seconds (interactive)
*!   - Updates at 10% milestones only (batch)
*!   - Always updates on the first and last iteration
*!
*! Timer usage: Uses Mata timer slot 99 (reserved for progress bar)
*! to avoid conflicts with user timers (typically 1-10).
*!
*! Usage:
*!   struct _eqt_progress_bar scalar pb
*!   _eqt_progress_init(pb, 1000, "Bootstrap replications", 1)
*!   for (b = 1; b <= 1000; b++) {
*!       // ... computation ...
*!       _eqt_progress_update(pb, b)
*!   }
*!   _eqt_progress_finish(pb)

version 16.0

mata:
mata set matastrict on
mata set mataoptimize on

// ============================================================================
// STRUCTURE: _eqt_progress_bar
// Container for progress bar state
// ============================================================================
struct _eqt_progress_bar {
    real scalar total              // Total number of iterations
    real scalar current            // Current iteration
    real scalar show               // Whether to display progress (1=yes, 0=no)
    real scalar bar_width          // Width of the visual bar in characters
    real scalar last_update_sec    // Elapsed seconds at last display update
    real scalar update_interval    // Minimum seconds between display updates
    real scalar last_displayed_pct // Last displayed percentage (avoid redundant updates)
    string scalar title            // Title text displayed above the bar
    real scalar timer_slot         // Mata timer slot used (99)
    real scalar batch_mode         // Whether running in batch mode (1=yes, 0=no)
    real scalar last_milestone     // Last milestone percentage displayed (batch mode)
}


// ============================================================================
// _eqt_format_time()
// Format elapsed seconds into MM:SS or HH:MM:SS string
//
// Arguments:
//   seconds : real scalar - elapsed time in seconds
//
// Returns:
//   string scalar - formatted time string
// ============================================================================
string scalar _eqt_format_time(real scalar seconds)
{
    real scalar s, m, h
    string scalar result

    if (seconds < 0) seconds = 0
    s = floor(seconds)

    if (s < 3600) {
        m = floor(s / 60)
        s = s - m * 60
        result = sprintf("%02.0f:%02.0f", m, s)
    }
    else {
        h = floor(s / 3600)
        s = s - h * 3600
        m = floor(s / 60)
        s = s - m * 60
        result = sprintf("%02.0f:%02.0f:%02.0f", h, m, s)
    }

    return(result)
}


// ============================================================================
// _eqt_format_number()
// Format an integer with comma separators (e.g. 1000 -> "1,000")
//
// Arguments:
//   n : real scalar - number to format
//
// Returns:
//   string scalar - formatted number string
// ============================================================================
string scalar _eqt_format_number(real scalar n)
{
    string scalar raw, result
    real scalar len, pos, i, digits_before_comma

    raw = strofreal(floor(n))
    len = strlen(raw)

    if (len <= 3) return(raw)

    result = ""
    // How many digits before the first comma (1-3)
    digits_before_comma = mod(len, 3)
    if (digits_before_comma == 0) digits_before_comma = 3

    pos = 1
    for (i = 1; i <= len; i = i + 1) {
        if (i > digits_before_comma && mod(i - digits_before_comma - 1, 3) == 0) {
            result = result + ","
        }
        result = result + substr(raw, i, 1)
        pos = pos + 1
    }

    return(result)
}


// ============================================================================
// _eqt_progress_init()
// Initialize progress bar state and display header
//
// Arguments:
//   pb    : struct _eqt_progress_bar scalar - progress bar state (modified)
//   total : real scalar - total number of iterations
//   title : string scalar - title text
//   show  : real scalar - whether to display (1=yes, 0=no)
// ============================================================================
void _eqt_progress_init(struct _eqt_progress_bar scalar pb,
                         real scalar total, string scalar title,
                         real scalar show)
{
    pb.total = total
    pb.current = 0
    pb.show = show
    pb.bar_width = 30
    pb.last_update_sec = 0
    pb.update_interval = 0.3
    pb.last_displayed_pct = -1
    pb.title = title
    pb.timer_slot = 99
    pb.last_milestone = 0

    // Detect batch mode: c(mode) returns "batch" when running stata -b
    if (st_global("c(mode)") == "batch") {
        pb.batch_mode = 1
    }
    else {
        pb.batch_mode = 0
    }

    if (show == 0) return

    // Start timer
    timer_clear(pb.timer_slot)
    timer_on(pb.timer_slot)

    // Display header
    if (pb.batch_mode) {
        printf("%s (%s)\n", title, _eqt_format_number(total))
    }
    else {
        printf("%s (%s)\n", title, _eqt_format_number(total))
    }
    displayflush()
}


// ============================================================================
// _eqt_progress_update()
// Update progress bar display
//
// In interactive mode: overwrites current line with carriage return
// In batch mode: prints milestone lines (every 10%) with newlines
//
// Arguments:
//   pb      : struct _eqt_progress_bar scalar - progress bar state (modified)
//   current : real scalar - current iteration number
// ============================================================================
void _eqt_progress_update(struct _eqt_progress_bar scalar pb,
                           real scalar current)
{
    real scalar elapsed, pct, filled, empty, eta
    real scalar milestone
    string scalar bar_filled, bar_empty, line, time_str, eta_str
    string scalar counter_str
    real scalar i

    pb.current = current

    if (pb.show == 0) return

    // Get elapsed time
    timer_off(pb.timer_slot)
    elapsed = timer_value(pb.timer_slot)[1]
    timer_on(pb.timer_slot)

    pct = floor(current / pb.total * 100)

    // ---- Batch mode: milestone-based display ----
    if (pb.batch_mode) {
        // Display at 10% milestones and at iteration 1
        milestone = floor(pct / 10) * 10
        if (current == 1) {
            // First iteration
            line = sprintf("  %3.0f", pct) + "%" +
                   " | " + _eqt_format_number(current) + "/" +
                   _eqt_format_number(pb.total) +
                   " | Elapsed: " + _eqt_format_time(elapsed)
            printf("%s\n", line)
            displayflush()
            pb.last_milestone = 0
        }
        else if (milestone > pb.last_milestone && milestone > 0 && milestone < 100) {
            // Milestone reached (skip 100% — handled by _eqt_progress_finish)
            time_str = _eqt_format_time(elapsed)
            if (pct > 0 && pct < 100) {
                eta = elapsed * (100 - pct) / pct
                eta_str = _eqt_format_time(eta)
            }
            else {
                eta_str = "--:--"
            }
            line = sprintf("  %3.0f", milestone) + "%" +
                   " | " + _eqt_format_number(current) + "/" +
                   _eqt_format_number(pb.total) +
                   " | Elapsed: " + time_str +
                   " | ETA: " + eta_str
            printf("%s\n", line)
            displayflush()
            pb.last_milestone = milestone
        }
        return
    }

    // ---- Interactive mode: carriage-return progress bar ----

    // Throttle updates: skip if less than update_interval since last update
    // Always display first and last iteration
    if (current > 1 && current < pb.total) {
        if (elapsed - pb.last_update_sec < pb.update_interval) {
            return
        }
    }
    pb.last_update_sec = elapsed

    // Build visual bar
    filled = round(pct / 100 * pb.bar_width)
    empty = pb.bar_width - filled

    bar_filled = ""
    for (i = 1; i <= filled; i = i + 1) {
        bar_filled = bar_filled + "#"
    }
    bar_empty = ""
    for (i = 1; i <= empty; i = i + 1) {
        bar_empty = bar_empty + "-"
    }

    // Build counter string
    counter_str = _eqt_format_number(current) + "/" + _eqt_format_number(pb.total)

    // Build time strings
    time_str = _eqt_format_time(elapsed)

    if (pct > 0 && pct < 100) {
        eta = elapsed * (100 - pct) / pct
        eta_str = " | ETA: " + _eqt_format_time(eta)
    }
    else if (pct >= 100) {
        eta_str = ""
    }
    else {
        eta_str = ""
    }

    // Assemble full line using string concatenation (safe from printf % issues)
    line = char(13) + "[" + bar_filled + bar_empty + "]" +
           sprintf(" %3.0f", pct) + "%" +
           " | " + counter_str +
           " | Elapsed: " + time_str + eta_str

    // Output using %s to prevent printf from interpreting % in the string
    printf("%s", line)
    displayflush()

    pb.last_displayed_pct = pct
}


// ============================================================================
// _eqt_progress_finish()
// Finalize progress bar display
//
// In interactive mode: shows 100% bar and moves to next line
// In batch mode: prints final 100% milestone
//
// Arguments:
//   pb : struct _eqt_progress_bar scalar - progress bar state
// ============================================================================
void _eqt_progress_finish(struct _eqt_progress_bar scalar pb)
{
    real scalar elapsed
    string scalar bar_filled, line, time_str
    real scalar i

    if (pb.show == 0) return

    // Stop timer and get final elapsed time
    timer_off(pb.timer_slot)
    elapsed = timer_value(pb.timer_slot)[1]
    time_str = _eqt_format_time(elapsed)

    if (pb.batch_mode) {
        // Final milestone line
        line = "  100" + "%" +
               " | " + _eqt_format_number(pb.total) + "/" +
               _eqt_format_number(pb.total) +
               " | Elapsed: " + time_str
        printf("%s\n", line)
        displayflush()
    }
    else {
        // Build full 100% bar
        bar_filled = ""
        for (i = 1; i <= pb.bar_width; i = i + 1) {
            bar_filled = bar_filled + "#"
        }

        line = char(13) + "[" + bar_filled + "]" +
               " 100" + "%" +
               " | " + _eqt_format_number(pb.total) + "/" +
               _eqt_format_number(pb.total) +
               " | Elapsed: " + time_str

        printf("%s", line)
        printf("\n")
        displayflush()
    }
}


// ============================================================================
// SEARCH PROGRESS FUNCTIONS
// For golden section search in _eqt_min_delta_bootstrap
// Uses timer slot 98 to avoid conflict with bootstrap progress (slot 99)
// ============================================================================

// ============================================================================
// _eqt_progress_search_init()
// Initialize search progress display
//
// Arguments:
//   max_iter : real scalar - maximum iterations
//   title    : string scalar - title text
//   show     : real scalar - whether to display (1=yes, 0=no)
// ============================================================================
void _eqt_progress_search_init(real scalar max_iter, string scalar title,
                                real scalar show)
{
    if (show == 0) return

    timer_clear(98)
    timer_on(98)

    printf("\n%s (max %s iterations)\n", title, _eqt_format_number(max_iter))
    displayflush()
}


// ============================================================================
// _eqt_progress_search_update()
// Update search progress display
//
// Arguments:
//   iter     : real scalar - current iteration
//   max_iter : real scalar - maximum iterations
//   gap      : real scalar - current search interval width
//   show     : real scalar - whether to display (1=yes, 0=no)
// ============================================================================
void _eqt_progress_search_update(real scalar iter, real scalar max_iter,
                                  real scalar gap, real scalar show)
{
    real scalar elapsed, is_batch
    string scalar line

    if (show == 0) return

    timer_off(98)
    elapsed = timer_value(98)[1]
    timer_on(98)

    // Detect batch mode
    is_batch = (st_global("c(mode)") == "batch")

    // Build line via string concatenation (safe from printf % issues)
    line = "  Iteration " + strofreal(iter) + "/" + strofreal(max_iter) +
           " | Gap: " + sprintf("%12.8f", gap) +
           " | Elapsed: " + _eqt_format_time(elapsed)

    if (is_batch) {
        // In batch mode, only print every 5 iterations or first/last
        if (iter == 1 || mod(iter, 5) == 0) {
            printf("%s\n", line)
            displayflush()
        }
    }
    else {
        printf("%s%s", char(13), line)
        displayflush()
    }
}


// ============================================================================
// _eqt_progress_search_finish()
// Finalize search progress display
//
// Arguments:
//   iter      : real scalar - final iteration count
//   converged : real scalar - whether search converged (1=yes, 0=no)
//   show      : real scalar - whether to display (1=yes, 0=no)
// ============================================================================
void _eqt_progress_search_finish(real scalar iter, real scalar converged,
                                  real scalar show)
{
    real scalar elapsed
    string scalar line, status

    if (show == 0) return

    timer_off(98)
    elapsed = timer_value(98)[1]

    if (converged) {
        status = "Converged"
    }
    else {
        status = "Max iterations reached"
    }

    line = "  " + status + " in " + strofreal(iter) +
           " iterations | Elapsed: " + _eqt_format_time(elapsed)

    printf("\n%s\n", line)
    displayflush()
}

end
