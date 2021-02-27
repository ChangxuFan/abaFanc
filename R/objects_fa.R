add.fa <- function(df, outdir) {
  system(paste0("mkdir -p ", outdir))
  df$fa <- paste0(outdir, "/", df$acc, ".fasta")
  for (i in 1:nrow(df)) {
    
    fa <- entrez_fetch("nuccore", id = df[i, "acc"], rettype = "fasta", api_key ="4d2374d16cee57dded1296ac1a48ba9c3b09")
    write(fa, df[i, "fa"])
  }
  
  trash <- lapply(df$fa, function(x) {
    if(!file.exists(x)) {
      stop(paste0(x, " was not successfully created"))
    }
  })
  return(df)
}