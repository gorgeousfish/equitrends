*! equitrends_data.ado - Load bundled example datasets for equitrends package
*! Version 2.0.0  2026-02-08

program define equitrends_data
    version 16.0
    
    syntax [anything(name=dataset)] [, clear]
    
    // Default dataset
    if `"`dataset'"' == "" {
        local dataset "MonthlyPanel"
    }
    
    // Normalize dataset name (remove .dta extension if present)
    local dataset = subinstr(`"`dataset'"', ".dta", "", .)
    
    // Locate the dataset file from installed package directories
    _equitrends_find_data `"`dataset'"'
    local filepath `"`s(filepath)'"'
    local searched `"`s(searched)'"'
    
    if `"`filepath'"' == "" {
        display as error `"dataset `dataset'.dta not found"'
        display as error ""
        display as error "Searched in:"
        // Parse semicolon-delimited searched paths
        local remaining `"`searched'"'
        while `"`remaining'"' != "" {
            gettoken path remaining : remaining, parse(";")
            if `"`path'"' != ";" & `"`path'"' != "" {
                display as error `"  - `path'"'
            }
        }
        display as error ""
        display as error "Available datasets: MonthlyPanel"
        display as error ""
        display as error "The package may not be correctly installed."
        display as error "Please reinstall using: net install equitrends"
        exit 601
    }
    
    // Load the dataset
    use `"`filepath'"', `clear'
    
    display as text "(`dataset'.dta loaded from equitrends package)"
end


// =========================================================================
// Internal helper: locate dataset file across installation contexts
// =========================================================================

program define _equitrends_find_data, sclass
    version 16.0
    
    args dataset
    
    sreturn clear
    local filepath ""
    local searched ""
    
    local plus_dir `"`c(sysdir_plus)'"'
    
    // Strategy 1: PLUS/d/ (net install places data/*.dta files under d/)
    mata: st_local("try_path", pathjoin(pathjoin(st_local("plus_dir"), "d"), st_local("dataset") + ".dta"))
    local searched `"`try_path'"'
    capture confirm file `"`try_path'"'
    if _rc == 0 {
        local filepath `"`try_path'"'
    }
    
    // Strategy 2: PLUS/data/ (explicit data subdirectory)
    if `"`filepath'"' == "" {
        mata: st_local("try_path", pathjoin(pathjoin(st_local("plus_dir"), "data"), st_local("dataset") + ".dta"))
        local searched `"`searched';`try_path'"'
        capture confirm file `"`try_path'"'
        if _rc == 0 {
            local filepath `"`try_path'"'
        }
    }
    
    // Strategy 3 & 4: Search all adopath entries (flat and data/ subdirectory)
    // S_ADO is semicolon-delimited and may contain codewords (BASE, SITE,
    // PLUS, PERSONAL, OLDPLACE) that must be resolved to actual paths.
    // All path construction uses Mata pathjoin() for cross-platform safety.
    if `"`filepath'"' == "" {
        local ado_raw `"$S_ADO"'
        while `"`ado_raw'"' != "" {
            gettoken token ado_raw : ado_raw, parse(";")
            if `"`token'"' == ";" {
                continue
            }
            
            // Resolve codewords to actual directory paths
            local path `"`token'"'
            if `"`path'"' == "BASE"     local path `"`c(sysdir_base)'"'
            if `"`path'"' == "SITE"     local path `"`c(sysdir_site)'"'
            if `"`path'"' == "PLUS"     local path `"`c(sysdir_plus)'"'
            if `"`path'"' == "PERSONAL" local path `"`c(sysdir_personal)'"'
            if `"`path'"' == "OLDPLACE" local path `"`c(sysdir_oldplace)'"'
            
            // Flat search: path/dataset.dta
            mata: st_local("try_path", pathjoin(st_local("path"), st_local("dataset") + ".dta"))
            local searched `"`searched';`try_path'"'
            capture confirm file `"`try_path'"'
            if _rc == 0 {
                local filepath `"`try_path'"'
                continue, break
            }
            
            // Subdirectory search: path/data/dataset.dta
            mata: st_local("try_path", pathjoin(pathjoin(st_local("path"), "data"), st_local("dataset") + ".dta"))
            local searched `"`searched';`try_path'"'
            capture confirm file `"`try_path'"'
            if _rc == 0 {
                local filepath `"`try_path'"'
                continue, break
            }
        }
    }
    
    // Strategy 5: Locate via the ado file's own directory
    // Use -findfile- to find where equitrends_data.ado lives, then search
    // sibling data/ and ancestor ../data/, ../../data/ directories.
    // This handles development (adopath ++) and non-standard installs.
    // All path manipulation uses Mata pathjoin/pathgetparent for
    // cross-platform compatibility (Windows backslash, macOS/Linux slash).
    if `"`filepath'"' == "" {
        capture findfile equitrends_data.ado
        if _rc == 0 {
            local ado_location `"`r(fn)'"'
            
            // ado_dir = directory containing equitrends_data.ado
            mata: st_local("ado_dir", pathgetparent(st_local("ado_location")))
            
            // Search: ado_dir/data/dataset.dta
            mata: st_local("try_path", pathjoin(pathjoin(st_local("ado_dir"), "data"), st_local("dataset") + ".dta"))
            local searched `"`searched';`try_path'"'
            capture confirm file `"`try_path'"'
            if _rc == 0 {
                local filepath `"`try_path'"'
            }
            
            // Search: ado_dir/../data/dataset.dta
            if `"`filepath'"' == "" {
                mata: st_local("try_path", pathjoin(pathjoin(pathgetparent(st_local("ado_dir")), "data"), st_local("dataset") + ".dta"))
                local searched `"`searched';`try_path'"'
                capture confirm file `"`try_path'"'
                if _rc == 0 {
                    local filepath `"`try_path'"'
                }
            }
            
            // Search: ado_dir/../../data/dataset.dta
            if `"`filepath'"' == "" {
                mata: st_local("try_path", pathjoin(pathjoin(pathgetparent(pathgetparent(st_local("ado_dir"))), "data"), st_local("dataset") + ".dta"))
                local searched `"`searched';`try_path'"'
                capture confirm file `"`try_path'"'
                if _rc == 0 {
                    local filepath `"`try_path'"'
                }
            }
        }
    }
    
    sreturn local filepath `"`filepath'"'
    sreturn local searched `"`searched'"'
end
