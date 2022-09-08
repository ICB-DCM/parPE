function(split_version versionstring libname major minor patch)
    string(REGEX REPLACE "([A-Za-z0-9_]*)-[vV].*" "\\1" locallibname ${versionstring} )
    set(libname ${locallibname} PARENT_SCOPE)

    string(REGEX REPLACE "^([A-Za-z0-9_]*-[vV])([0-9]*)([.][0-9]*[.][0-9]*-?.*)$" "\\2" numbers ${versionstring} )
    set(major ${numbers} PARENT_SCOPE)

    string(REGEX REPLACE "^([A-Za-z0-9_]*-[vV][0-9]*[.])([0-9]*)([.][0-9]*-?.*)$" "\\2" numbers ${versionstring} )
    set(minor ${numbers} PARENT_SCOPE)

    string(REGEX REPLACE "^([A-Za-z0-9_]*-[vV][0-9]*[.][0-9]*[.])([0-9]*)(-?.*)$" "\\2" numbers ${versionstring} )
    set(patch ${numbers} PARENT_SCOPE)
endfunction()


