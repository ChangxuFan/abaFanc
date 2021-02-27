liftover.core <- function(in.dir, in.pattern, chain.file, out.dir, zip = T, run = T, 
                          liftOver = "/opt/apps/kentUCSC/334/liftOver") {
  dir.create(out.dir, recursive = T, showWarnings = F)
  in.vec <- list.files(in.dir) %>% .[grepl(in.pattern, .)]
  lapply(in.vec, function(x) {
    in.file <- file.path(in.dir, x)
    suffix <- chain.file %>% basename() %>% sub(".over.+", "", .) %>% sub(".to.", "_to_", .)
    ext <- tools::file_ext(x)
    root <- tools::file_path_sans_ext(x)
    out.file <- paste0(out.dir, "/", root, "_", suffix, ".", ext)
    cmd <- paste0(liftOver, " ", in.file, " ", chain.file, " ", out.file, " /dev/null")
    print(cmd)
    if (run == T) {
      system(cmd)
      if (zip == T) {
        cmd <- paste0("~/scripts/bed_browser_v2.sh ", out.file)
        print(cmd)
        system(cmd)
      }
    }
      
    
    return()
  }) 
  return()
}
