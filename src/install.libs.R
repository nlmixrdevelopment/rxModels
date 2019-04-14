## Similar to default behavior according to R extension manual
if (FALSE){
    ## No binaries are created in this package
    files <- Sys.glob(paste0("*", SHLIB_EXT))
    dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    file.copy(files, dest, overwrite = TRUE)
    if(file.exists("symbols.rds"))
        file.copy("symbols.rds", dest, overwrite = TRUE)
}
## This pre-compiled RxODE binaries to the correct location.
RxODE::rxUse()
