# a generalized version of blastn dispatching and parsing.
blast.table.to.bed <- function(blast.table, out.root.name = NULL, expand = F, squish = F, shift= 0,
                               genome = NULL, filter.func = NULL, filter.out.name = NULL, 
                               igv.gen = T, push.2.igv = F, simple = F, ...) {
  if (is.null(out.root.name))
    out.root.name <- sub(".txt", "", blast.table)
  if (is.character(blast.table))
    blast.table <- read.table(blast.table, header = F, sep = "\t", quote = "", as.is = T)
  if (!is.null(genome))
    blast.table <- blast.table %>%  filter(! is.na (chr2refseqID(V2, genome = genome)))
  
  if (!is.null(filter.func))
    blast.table <- filter.func(blast.table)
  
  if (!is.null(filter.out.name))
    write.table(blast.table, filter.out.name, sep = "\t", row.names = F, col.names = F, quote = F)
  
  sub.f.print <- function(b.tbl, out.bed, squish) {
    bed <- data.frame(chr = b.tbl$V2,
                      left = pmin(b.tbl$V10, b.tbl$V9) -1 + shift,
                      right = pmax(b.tbl$V10, b.tbl$V9) + shift,
                      text = b.tbl %>% 
                        mutate(text = if (simple == F) paste(V3, V4, V5, V6, V7, V8, V12, sep = "_") else sub(".+$","",V3)) %>% 
                        pull(text),
                      fifth = ".", 
                      strand = sapply(b.tbl$V10 - b.tbl$V9, function(x) {
                        return(if_else(x > 0, "+", "-"))
                      }))
    if (squish == T) {
      bed$text <- paste0(b.tbl$V1, ":", bed$text)
    }
    if (!is.null(genome))
      bed$chr <- bed$chr %>% chr2refseqID(genome = genome)
    
    write.table(bed, out.bed, sep = "\t", quote = F, row.names = F, col.names = F)
    system(paste0("/bar/cfan/scripts/bed_browser_v2.sh ", out.bed))
    return(bed)
  }
  
  out.bed.list <- c()
  
  if (squish == T) {
    out.bed <- paste0(out.root.name, ".bed")
    # print(out.bed)
    bed <- sub.f.print(b.tbl = blast.table, out.bed = out.bed, squish = T)
    out.bed.list <- c(out.bed.list, out.bed)
    # return(bed)
  }
  
  if (expand ==  T) {
    out.beds <- blast.table %>% split(., f=factor(.$V1, levels = unique(.$V1))) %>% 
      sapply(function(bt.sub) {
        query <- bt.sub$V1[1]
        out.bed <- paste0(out.root.name, "_", query, ".bed")
        sub.f.print(b.tbl = bt.sub, out.bed = out.bed, squish = F)
        return(out.bed)
      })
    out.bed.list <- c(out.bed.list, out.beds)
  }
  
  if (igv.gen == T) {
    igv.hubgen(files.vec = out.bed.list, no.writing = !push.2.igv, ...)
  }
  
  return(out.bed.list)
}

blastn.fanc <- function(query.fa.vec, subject.fa.vec, outtable, 
                        blastn = "/opt/apps/blast/2.10.1/bin/blastn", other.params = "") {
  system(paste0("mkdir -p ", dirname(outtable)))
  if (length(query.fa.vec) > 1) {
    query <- paste0(outtable, ".prep.query")
    cmd <- paste0("cat ", paste0(query.fa.vec, collapse = " "), " > ", query)
    print(cmd)
    system(cmd)
  } else {
    query <- query.fa.vec
  }
  
  if (length(subject.fa.vec) > 1) {
    subject <- paste0(outtable, ".prep.subject")
    cmd <- paste0("cat ", paste0(subject.fa.vec, collapse = " "), " > ", subject)
    print(cmd)
    system(cmd)
  } else {
    subject <- subject.fa.vec
  }
  
  cmd <- paste0(blastn, " -query ",query, " -subject ", subject, " -outfmt 7 ", " -out ", outtable, " ", other.params )
  print(cmd)
  system(cmd)
  return(outtable)
}



blast.filter.longest <- function(x) {
  res <- x %>% group_by(V1, V2) %>% slice_max(order_by = V4, n = 1) %>% ungroup()
  return(res)
}

blast.filter.best.match <- function(x, force.n = NULL) {
  res <- x %>% group_by(V1, V2) %>% slice_max(order_by = V3, n = 1) %>% ungroup()
  return(res)
}

blast.filter.length <- function(x, threshold=100) {
  res <- x %>% filter(V4 > threshold)
  return(res)
}