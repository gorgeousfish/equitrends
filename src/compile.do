*! compile.do - Mata compilation script for EQUITRENDS
*!
*! Author: Cai Xuanyu, Xu Wenli
*! Version: 0.1.0

version 16.0
clear all
set more off

* ============================================
* Step 0: Change to src directory
* ============================================
* Try to detect the src directory from the do-file path
* In batch mode, c(filename) may be empty, so we use a fallback approach

local src_dir ""

* Method 1: Try c(filename) first
local this_file "`c(filename)'"
if "`this_file'" != "" {
    * Extract directory from full path
    local src_dir = substr("`this_file'", 1, strrpos("`this_file'", "/") - 1)
    if "`src_dir'" == "" {
        local src_dir = substr("`this_file'", 1, strrpos("`this_file'", "\") - 1)
    }
}

* Method 2: If c(filename) is empty, check if we're already in src directory
if "`src_dir'" == "" {
    capture confirm file "mata/equitrends_utils.mata"
    if _rc == 0 {
        local src_dir "`c(pwd)'"
    }
}

* Method 3: Try common project paths (removed - use relative paths instead)
if "`src_dir'" == "" {
    display as error "ERROR: Could not find src directory"
    display as error "Please run this script from the src directory or specify the path"
    exit 601
}

* Change to src directory if found
quietly cd "`src_dir'"
display as text "Working directory: `src_dir'"

display as text ""
display as text "============================================"
display as text "EQUITRENDS Mata Compilation"
display as text "============================================"
display as text ""

* ============================================
* Step 1: Set compilation options
* ============================================
display as text "Step 1: Setting compilation options..."

set matastrict on
set mataoptimize on
set matalnum off

display as text "  matastrict: on"
display as text "  mataoptimize: on"
display as text "  matalnum: off"

* ============================================
* Step 2: Create build directory
* ============================================
display as text ""
display as text "Step 2: Creating build directory..."

capture mkdir "build"
display as text "  build/ directory ready."

* ============================================
* Step 3: Clear any existing Mata library
* ============================================
display as text ""
display as text "Step 3: Clearing existing library..."

capture mata: mata mlib index
capture erase "build/lequitrends.mlib"
display as text "  Previous library cleared."

* ============================================
* Step 4: Compile Mata source files in order
* ============================================
display as text ""
display as text "Step 4: Compiling Mata source files..."

local compile_error = 0

* Define compilation order (dependency order)
* Note: equitrends_validate must come before equitrends (used by validation)
* Note: equitrends_dataproc must come after equitrends_validate (uses validation functions)
* Note: equitrends_demean must come before equitrends_rms (used by RMS functions)
* Note: equitrends_foldnorm provides folded normal distribution functions
* Note: equitrends_multicollin provides multicollinearity detection
* Note: equitrends_progress provides progress bar display (used by bootstrap and min_delta)
* Note: equitrends_min_delta provides minimum threshold calculation
* Note: equitrends_mean provides mean test functions
* Note: equitrends_maxtest_core provides max test core functions (IU and Bootstrap)
* Note: equitrends_meantest_core provides mean test core functions
* Note: equitrends_rms must come after equitrends_dataproc (uses data processing functions)
* Note: equitrends_sim has no dependencies on other equitrends modules
local mata_files "equitrends_utils equitrends_validate equitrends_demean equitrends_dataproc equitrends_vcov equitrends_foldnorm equitrends_multicollin equitrends_progress equitrends_bootstrap equitrends_min_delta equitrends_mean equitrends_maxtest_core equitrends_meantest_core equitrends_rms equitrends_sim equitrends"

foreach mfile of local mata_files {
    local mata_path "mata/`mfile'.mata"
    
    capture confirm file "`mata_path'"
    if _rc == 0 {
        display as text "  Compiling `mfile'.mata..."
        capture noisily do "`mata_path'"
        if _rc != 0 {
            display as error "    ERROR: Failed to compile `mfile'.mata"
            local compile_error = 1
        }
        else {
            display as text "    OK"
        }
    }
    else {
        display as text "  Skipping `mfile'.mata (not found - will be created in later stories)"
    }
}

* ============================================
* Step 5: Create Mata library
* ============================================
display as text ""
display as text "Step 5: Creating Mata library..."

if `compile_error' == 0 {
    * Check if any Mata functions were compiled
    capture mata: mata describe Equitrends*()
    local has_functions = (_rc == 0)
    
    capture mata: mata describe _equitrends*()
    if _rc == 0 {
        local has_functions = 1
    }
    
    if `has_functions' {
        * Create the library with exported functions
        * Export patterns: Equitrends*() and _equitrends*()
        capture noisily mata: mata mlib create lequitrends, dir("build") replace
        
        if _rc == 0 {
            * Add functions matching the export patterns
            capture noisily mata: mata mlib add lequitrends _equitrends*(), dir("build")
            
            * Add _eqt_* functions (utility functions)
            capture noisily mata: mata mlib add lequitrends _eqt_*(), dir("build")
            
            * Add _maxequivtest_* functions (max test core functions)
            capture noisily mata: mata mlib add lequitrends _maxequivtest_*(), dir("build")
            
            * Add _meanequivtest_* functions (mean test core functions)
            capture noisily mata: mata mlib add lequitrends _meanequivtest_*(), dir("build")
            
            * Add classes - must be added explicitly
            capture noisily mata: mata mlib add lequitrends EquitrendsBase(), dir("build")
            capture noisily mata: mata mlib add lequitrends EquitrendsDataProcessor(), dir("build")
            capture noisily mata: mata mlib add lequitrends MaxEquivTest(), dir("build")
            capture noisily mata: mata mlib add lequitrends MeanEquivTest(), dir("build")
            capture noisily mata: mata mlib add lequitrends RMSEquivTest(), dir("build")
            
            display as text "  Library created: build/lequitrends.mlib"
        }
        else {
            display as error "  ERROR: Failed to create library"
            local compile_error = 1
        }
    }
    else {
        display as text "  No Mata functions to export (skeleton build)."
        display as text "  Library will be created when Mata code is added."
    }
}
else {
    display as error "  Skipping library creation due to compilation errors."
}

* ============================================
* Step 6: Install library to PLUS directory
* ============================================
display as text ""
display as text "Step 6: Installing library to PLUS directory..."

* Get the PLUS directory path
local plus_dir "`c(sysdir_plus)'"
display as text "  PLUS directory: `plus_dir'"

* Copy the library to PLUS directory (so it's found after clear all)
capture copy "build/lequitrends.mlib" "`plus_dir'l/lequitrends.mlib", replace
if _rc == 0 {
    display as text "  Library installed to: `plus_dir'l/lequitrends.mlib"
}
else {
    * Try creating the l/ subdirectory first
    capture mkdir "`plus_dir'l"
    capture copy "build/lequitrends.mlib" "`plus_dir'l/lequitrends.mlib", replace
    if _rc == 0 {
        display as text "  Library installed to: `plus_dir'l/lequitrends.mlib"
    }
    else {
        display as text "  Warning: Could not install to PLUS directory (rc=`=_rc')"
        display as text "  Library remains in: build/lequitrends.mlib"
    }
}

* ============================================
* Step 7: Rebuild library index
* ============================================
display as text ""
display as text "Step 7: Rebuilding library index..."

capture noisily mata: mata mlib index
if _rc == 0 {
    display as text "  Library index updated."
    
    * Verify lequitrends is in the index
    capture mata: mata mlib query lequitrends
    if _rc == 0 {
        display as text "  [OK] lequitrends library is now in the search path"
    }
    else {
        display as text "  Warning: lequitrends not found in library index"
    }
}
else {
    display as text "  Warning: Could not update library index."
}

* ============================================
* Compilation Complete
* ============================================
display as text ""
display as text "============================================"
if `compile_error' == 0 {
    display as result "Compilation completed successfully!"
    display as text "  Library location: `plus_dir'l/lequitrends.mlib"
    display as text "  Functions will persist after 'clear all'"
}
else {
    display as error "Compilation completed with errors."
}
display as text "============================================"
display as text ""

* List compiled functions (if any)
display as text "Exported functions:"
capture noisily mata: mata describe Equitrends*()
capture noisily mata: mata describe _equitrends*()

if _rc != 0 {
    display as text "  (No functions exported yet - skeleton build)"
}

display as text ""
