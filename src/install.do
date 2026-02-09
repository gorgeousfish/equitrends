*! install.do - Local development installation script for EQUITRENDS
*! 
*! This script installs the equitrends package from the local development
*! directory. It should be run from the src/ directory.
*!
*! Usage:
*!   cd "/path/to/equitrends-stata/src"
*!   do install.do
*!
*! Requirements:
*!   - Stata 16.0 or higher
*!   - Python 3.6+ with numpy and scipy
*!
*! Author: EQUITRENDS Development Team
*! Version: 1.0.0
*! Date: 2026-01-08

version 16.0
clear all
set more off

display as text ""
display as text "============================================"
display as text "EQUITRENDS Local Development Installation"
display as text "============================================"
display as text ""

* ============================================
* Step 1: Uninstall any existing installation
* ============================================
display as text "Step 1: Removing existing installation..."

capture ado uninstall equitrends
if _rc == 0 {
    display as text "  Previous installation removed."
}
else {
    display as text "  No previous installation found."
}

* ============================================
* Step 2: Set Mata compilation options
* ============================================
display as text ""
display as text "Step 2: Setting Mata compilation options..."

set matastrict on
set mataoptimize on
display as text "  matastrict: on"
display as text "  mataoptimize: on"

* ============================================
* Step 3: Compile Mata source files
* ============================================
display as text ""
display as text "Step 3: Compiling Mata source files..."

* Check if Mata source files exist
local mata_files_exist = 0
capture confirm file "mata/equitrends_utils.mata"
if _rc == 0 {
    local mata_files_exist = 1
}

if `mata_files_exist' {
    * Run the compile script
    capture noisily do compile.do
    if _rc != 0 {
        display as error "  Mata compilation failed!"
        display as text "  Continuing with installation..."
    }
    else {
        display as text "  Mata compilation successful."
    }
}
else {
    display as text "  No Mata source files found (skeleton installation)."
    display as text "  Mata library will be created in later stories."
    
    * Create empty mlib placeholder
    capture mkdir "build"
    * Note: We can't create an empty mlib, so we skip this for now
}

* ============================================
* Step 4: Install package from local directory
* ============================================
display as text ""
display as text "Step 4: Installing package..."

* Get the parent directory (package root)
local current_dir = c(pwd)
local parent_dir = subinstr("`current_dir'", "/src", "", 1)
local parent_dir = subinstr("`parent_dir'", "\src", "", 1)

* Install from parent directory
capture noisily net install equitrends, from("`parent_dir'") replace force
if _rc != 0 {
    display as error "  Installation failed!"
    display as text ""
    display as text "  Trying alternative installation method..."
    
    * Alternative: manually copy files to PLUS directory
    local plus_dir : sysdir PLUS
    
    * Copy ado files
    local ado_files : dir "ado" files "*.ado"
    foreach f of local ado_files {
        capture copy "ado/`f'" "`plus_dir'`f'", replace
    }
    
    * Copy Python module
    capture copy "python/equitrends_python.py" "`plus_dir'equitrends_python.py", replace
    
    display as text "  Files copied to PLUS directory."
}
else {
    display as text "  Package installed successfully."
}

* ============================================
* Step 5: Verify installation
* ============================================
display as text ""
display as text "Step 5: Verifying installation..."

* Check if main programs are available
local verification_ok = 1

capture which _equitrends_check_deps
if _rc != 0 {
    display as error "  _equitrends_check_deps not found!"
    local verification_ok = 0
}
else {
    display as text "  _equitrends_check_deps: OK"
}

capture which _equitrends_load_python
if _rc != 0 {
    display as error "  _equitrends_load_python not found!"
    local verification_ok = 0
}
else {
    display as text "  _equitrends_load_python: OK"
}

* ============================================
* Step 6: Check dependencies
* ============================================
display as text ""
display as text "Step 6: Checking dependencies..."

capture noisily _equitrends_check_deps
if _rc != 0 {
    display as text ""
    display as error "Warning: Some dependencies are missing."
    display as text "Please install missing dependencies before using the package."
}

* ============================================
* Installation Complete
* ============================================
display as text ""
display as text "============================================"
if `verification_ok' {
    display as result "Installation completed successfully!"
}
else {
    display as error "Installation completed with warnings."
}
display as text "============================================"
display as text ""
display as text "To verify the installation, run:"
display as input "  . _equitrends_check_deps"
display as text ""
