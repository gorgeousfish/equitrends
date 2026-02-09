*! _equitrends_load_python.ado - Load EQUITRENDS Python module
*!
*! Syntax: _equitrends_load_python [, Force Quiet]
*!
*! Options:
*!   force - Force reload of equitrends_python module
*!   quiet - Suppress status messages
*!
*! Returns (r-class):
*!   r(python_available) - 1 if Python is available in Stata, 0 otherwise
*!   r(module_loaded)    - 1 if equitrends_python was imported, 0 otherwise
*!   r(python_path)      - Path added to sys.path for module import
*!
*! Exit codes:
*!   0   - Success
*!   198 - Python not available or module load failed

capture program drop _equitrends_load_python
program define _equitrends_load_python, rclass
    version 16.0
    syntax [, Force Quiet]

    local is_quiet = ("`quiet'" != "")
    local force_reload = ("`force'" != "")

    // -------------------------------------------------------------------------
    // Check Python availability
    // -------------------------------------------------------------------------
    capture python query
    if _rc != 0 {
        if !`is_quiet' {
            display as error "Error: Python is not available in Stata."
            display as text "Please configure Python: . python set exec <path>"
        }
        return scalar python_available = 0
        return scalar module_loaded = 0
        exit 198
    }

    // -------------------------------------------------------------------------
    // Locate Python module path
    // Strategy:
    //   1. Search S_ADO paths for development installations
    //   2. Check PLUS directory for formal installations
    //   3. Check current directory for development use
    // -------------------------------------------------------------------------
    local python_path ""
    
    // Method 1: Search S_ADO paths (for development: .../src/ado/ structure)
    foreach path of global S_ADO {
        capture confirm file "`path'/_equitrends_load_python.ado"
        if _rc == 0 {
            // Found ADO in S_ADO, check for python subdir
            local try_python = subinstr("`path'", "/ado", "/python", 1)
            capture confirm file "`try_python'/equitrends_python.py"
            if _rc == 0 {
                local python_path "`try_python'"
                continue, break
            }
        }
    }
    
    // Method 2: Check PLUS directory (for formal installation)
    if "`python_path'" == "" {
        local plus_dir : sysdir PLUS
        // Python module is installed to <PLUS>p/ directory
        capture confirm file "`plus_dir'p/equitrends_python.py"
        if _rc == 0 {
            local python_path "`plus_dir'p"
        }
    }
    
    // Method 3: Check current directory (for development from project root)
    if "`python_path'" == "" {
        local cwd "`c(pwd)'"
        capture confirm file "`cwd'/src/python/equitrends_python.py"
        if _rc == 0 {
            local python_path "`cwd'/src/python"
        }
    }
    
    // Method 4: Check if we're in src directory
    if "`python_path'" == "" {
        local cwd "`c(pwd)'"
        capture confirm file "`cwd'/python/equitrends_python.py"
        if _rc == 0 {
            local python_path "`cwd'/python"
        }
    }

    if "`python_path'" == "" {
        if !`is_quiet' {
            display as error "Error: Could not locate equitrends_python.py."
            display as error ""
            display as error "The module was searched in:"
            display as error "  1. S_ADO paths (development structure)"
            display as error "  2. PLUS directory: `c(sysdir_plus)'p/"
            display as error "  3. Current directory: `c(pwd)'/src/python/"
            display as error "  4. Current directory: `c(pwd)'/python/"
            display as error ""
            display as error "Solutions:"
            display as error "  - Run install.do from the package src/ directory"
            display as error "  - Or cd to the package root before running commands"
        }
        return scalar python_available = 1
        return scalar module_loaded = 0
        exit 198
    }

    // Ensure absolute path for Python
    if substr("`python_path'", 1, 1) != "/" {
        local python_path "`c(pwd)'/`python_path'"
    }

    // -------------------------------------------------------------------------
    // Load module using Python block syntax
    // -------------------------------------------------------------------------
    
    capture python: import sys
    if _rc != 0 {
        display as error "Error: Python sys module unavailable."
        exit 198
    }

    capture python: py_path = "`python_path'"
    if _rc != 0 {
        display as error "Error: Failed to set Python path string."
        exit 198
    }

    capture python: sys.path.insert(0, py_path) if (py_path and py_path not in sys.path) else None
    if _rc != 0 {
        display as error "Error: Failed to set Python path."
        exit 198
    }

    capture python: force_reload = (`force_reload' == 1)
    if _rc != 0 {
        display as error "Error: Failed to read reload flag."
        exit 198
    }

    capture python: sys.modules.pop('equitrends_python', None) if force_reload else None
    if _rc != 0 {
        display as error "Error: Failed to reset equitrends_python module."
        exit 198
    }

    capture python: import equitrends_python

    if _rc != 0 {
        if !`is_quiet' {
            display as error "Error: Failed to import equitrends_python."
            display as error "Python path attempted: `python_path'"
        }
        return scalar python_available = 1
        return scalar module_loaded = 0
        exit 198
    }

    if !`is_quiet' {
        display as text "Python module loaded: " as result "`python_path'"
    }

    return scalar python_available = 1
    return scalar module_loaded = 1
    return local python_path "`python_path'"
end
