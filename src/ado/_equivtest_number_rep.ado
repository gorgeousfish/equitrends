*! _equivtest_number_rep.ado - Number representation for display output
*!
*! Formats numeric values for consistent display output:
*!   - Values with |x| < 2e-16 are displayed as "<2e-16"
*!   - Values with |x| < 1e-4 use scientific notation (4 decimal places)
*!   - Other values use 4 significant digits

program define _equivtest_number_rep, rclass
    version 16.0
    args x
    
    local abs_x = abs(`x')
    
    if `abs_x' < 2e-16 {
        return local result "<2e-16"
    }
    else if `abs_x' < 1e-4 {
        // Scientific notation with 4 decimal places
        local result = string(`x', "%10.4e")
        local result = strtrim("`result'")
        return local result "`result'"
    }
    else {
        // 4 significant digits
        local result = string(`x', "%10.4g")
        local result = strtrim("`result'")
        return local result "`result'"
    }
end
