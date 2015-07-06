# -*- mode: cmake -*-

# 
# Jali Print variable
#

function(PRINT_VARIABLE VAR_NAME)
    message("==> " "${VAR_NAME}=${${VAR_NAME}}")
endfunction()    

