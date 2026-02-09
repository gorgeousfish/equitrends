*! equitrends_init.ado - Initialize EQUITRENDS Mata functions
*! Version 2.0 - On-demand compilation for net install compatibility
*!
*! This program ensures EQUITRENDS Mata functions are loaded by compiling
*! source files on-demand. This strategy:
*!   - Eliminates dependency on mlib indexing issues
*!   - Works reliably with net install
*!   - Compiles quickly (1-2 seconds) on first use
*!   - Functions persist in memory until clear mata/clear all

program define equitrends_init
    version 16.0
    
    // =========================================================================
    // Check if Mata functions already available (from mlib or memory)
    // =========================================================================
    
    capture mata: mata describe _equitrends_is_binary()
    if _rc == 0 {
        // Functions already available - no action needed
        exit 0
    }
    
    // =========================================================================
    // Try to use pre-compiled mlib (silent, fast)
    // =========================================================================
    
    // The mlib file should be in PLUS/l/ directory after net install
    local plus_dir "`c(sysdir_plus)'"
    capture confirm file "`plus_dir'l/lequitrends.mlib"
    if _rc == 0 {
        // mlib exists, rebuild index to make it available
        quietly mata: mata mlib index
        
        // Check if functions now available
        capture mata: mata describe _equitrends_is_binary()
        if _rc == 0 {
            // mlib loaded successfully, no compilation needed
            exit 0
        }
    }
    
    // =========================================================================
    // Locate Mata source files in installed package
    // =========================================================================
    
    // Search for .mata files
    // net install puts .mata files in e/ directory (first letter of filename)
    local plus_dir "`c(sysdir_plus)'"
    local mata_base ""
    
    // Strategy 1: Check PLUS/e/ (where net install puts equitrends*.mata)
    capture confirm file "`plus_dir'e/equitrends_utils.mata"
    if _rc == 0 {
        local mata_base "`plus_dir'e"
    }
    
    // Strategy 2: Check PLUS/src/mata/ (explicit directory in pkg)
    if "`mata_base'" == "" {
        capture confirm file "`plus_dir'src/mata/equitrends_utils.mata"
        if _rc == 0 {
            local mata_base "`plus_dir'src/mata"
        }
    }
    
    // Strategy 3: Search in all S_ADO paths
    if "`mata_base'" == "" {
        foreach path of global S_ADO {
            capture confirm file "`path'/equitrends_utils.mata"
            if _rc == 0 {
                local mata_base "`path'"
                continue, break
            }
        }
    }
    
    // Verify mata source found
    if "`mata_base'" == "" {
        display as error "EQUITRENDS: Mata source files not found"
        display as error ""
        display as error "Searched in:"
        display as error "  - `plus_dir'e/"
        display as error "  - `plus_dir'src/mata/"
        display as error "  - All S_ADO paths"
        display as error ""
        display as error "The package may not be correctly installed."
        display as error "Please reinstall using: net install equitrends"
        exit 601
    }
    
    // =========================================================================
    // Compile Mata source files silently (NO output unless error)
    // =========================================================================
    
    local cwd "`c(pwd)'"
    local compile_error = 0
    local compiled_count = 0
    
    quietly {
        capture cd "`mata_base'"
        if _rc != 0 {
            noisily display as error "EQUITRENDS: Could not find Mata source files"
            exit 601
        }
        
        // Compile in dependency order - COMPLETELY SILENT
        local mata_files "equitrends_utils equitrends_validate equitrends_demean"
        local mata_files "`mata_files' equitrends_dataproc equitrends_vcov"
        local mata_files "`mata_files' equitrends_foldnorm equitrends_multicollin"
        local mata_files "`mata_files' equitrends_progress equitrends_bootstrap equitrends_min_delta"
        local mata_files "`mata_files' equitrends_mean equitrends_maxtest_core"
        local mata_files "`mata_files' equitrends_meantest_core equitrends_rms"
        local mata_files "`mata_files' equitrends_sim equitrends"
        
        foreach mfile of local mata_files {
            capture confirm file "`mfile'.mata"
            if _rc == 0 {
                capture do "`mfile'.mata"
                if _rc != 0 {
                    local compile_error = 1
                    continue, break
                }
                local ++compiled_count
            }
        }
        
        cd "`cwd'"
    }
    
    if `compile_error' {
        display as error "EQUITRENDS: Mata compilation failed"
        exit 3000
    }
    
    // =========================================================================
    // Verify compilation success (silent if OK)
    // =========================================================================
    
    capture mata: mata describe _equitrends_is_binary()
    if _rc != 0 {
        display as error "EQUITRENDS: Initialization failed"
        exit 3499
    }
    
    // SUCCESS - No output, functions are ready to use
    
end
