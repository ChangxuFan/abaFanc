# #
SHORT.ALIGN.PARAMS <- " -word_size 7 -gapopen 5 -gapextend 2 -max_target_seqs 5000 -evalue 10"

fasta.wrap.fanc <- function(in.fa, out.fa=NULL) {
  if (is.null(out.fa))
    out.fa <- in.fa
  fa <- seqinr::read.fasta(in.fa, forceDNAtolower = F)
  headers <- lapply(fa, attr, "Annot") %>% unlist() %>% sub("^>", "", .)
  seqinr::write.fasta(fa, headers, out.fa)
  return(out.fa)
}

# fasta.mask.fanc <- function(in.fa, out.fa=NULL) {
#   if (is.null(out.fa))
#     out.fa <- in.fa
#   fa <- seqinr::read.fasta(in.fa, forceDNAtolower = F)
#   headers <- lapply(fa, attr, "Annot") %>% unlist() %>% sub("^>", "", .)
#   seqs <- lapply(fa, function(x) paste0(x[1:length(x)], collapse = "") ) %>%
#     unlist() %>% `names<-`(NULL) %>%
#     gsub("[a-z]", "", .)
#   seqinr::write.fasta(strsplit(seqs, ""), headers, out.fa)
#   return(out.fa)
# }


# bash2ftp <- function(filename) {
#   ftp <- sub("^~", "https://wangftp.wustl.edu/~cfan", filename) %>%
#     sub("/bar/cfan", "https://wangftp.wustl.edu/~cfan", .)
#   return(ftp)
# }

# bed.align <- function(bed, root.name=NULL, genome, align=F, mafft.options = "--auto --reorder",
#  add.coordinate=T, mask =F) {

#   if (is.null(root.name)) {
#     if (is.character(bed))
#       root.name <- sub(".bed", "", bed)
#     else
#       stop("root.name not specified")
#   }

#   if (is.character(bed))
#     bed <- read.table(bed, as.is = T, sep = "\t")

#   bed.out <- bed %>% mutate(id=1:nrow(.)) %>% split(., f=.$id) %>% lapply(function(x) {
#     x$id <- NULL
#     name.col <- ""
#     if (add.coordinate == T)
#       name.col <- paste0(x[1,1], ":", x[1,2], "-", x[1,3])
#     if (ncol(x) >= 4) {
#       if (name.col != "")
#         name.col <- paste0(name.col, "|")
#       name.col <- name.col %>% paste0(paste0(x[1,4:ncol(x)], collapse = "|"))
#     }
#     x[1,4] <- name.col
#     return(x[,1:4])
#   }) %>% Reduce(rbind,.)
  
#   write.table(bed.out, paste0(root.name, ".pre.bed"), sep = "\t", quote = F, col.names = F, row.names = F)
#   cmd <- paste0("/bar/cfan/anaconda2/envs/jupyter/bin/bedtools getfasta -fi /bar/cfan/genomes/",
#                 genome, "/", genome, ".fa -fo ", root.name, ".pre.fa -bed ", root.name, ".pre.bed -name")
#   print(cmd)
#   system(cmd)
  
#   fasta.wrap.fanc(in.fa = paste0(root.name, ".pre.fa"))

#   to.mafft <- paste0(root.name, ".pre.fa")

#   if (mask == T) {
#     masked <- fasta.mask.fanc(in.fa = to.mafft, out.fa = paste0(to.mafft, ".masked"))
#     to.mafft <- paste0(to.mafft, ".masked")
#   } 

#   if (align == T) {
#     after.mafft <- paste0(root.name, ".aligned.fa")
#     cmd <- paste0("/opt/apps/mafft/7.427/mafft ", mafft.options, " ", to.mafft, " > ", after.mafft)
#     print(cmd)
#     system(cmd)
#   } else {
#     after.mafft <- to.mafft
#   }

#   return(after.mafft)

# }

get.fasta.bed <- function(bed, root.name=NULL, genome=NULL, source.fasta=NULL, add.coordinate=T, add.additional.columns=T, return.fasta = F) {

  if (is.null(root.name)) {
    if (is.character(bed))
      root.name <- sub(".bed", "", bed)
    else
      stop("root.name not specified")
  }

  if (is.character(bed))
    bed <- read.table(bed, as.is = T, sep = "\t")

  bed.out <- bed %>% mutate(id=1:nrow(.)) %>% split(., f=.$id) %>% lapply(function(x) {
    x$id <- NULL
    name.col <- ""
    if (add.coordinate == T)
      name.col <- paste0(x[1,1], ":", x[1,2], "-", x[1,3])
    if (ncol(x) >= 4 && add.additional.columns==T) {
      if (name.col != "")
        name.col <- paste0(name.col, "|")
      name.col <- name.col %>% paste0(paste0(x[1,4:ncol(x)], collapse = "|"))
    }
    x[1,4] <- name.col
    return(x[,1:4])
  }) %>% Reduce(rbind,.)
  
  write.table(bed.out, paste0(root.name, ".pre.bed"), sep = "\t", quote = F, col.names = F, row.names = F)
  out.fa <- paste0(root.name, ".fa")
  if (!is.null(genome)) {
    fi <- paste0("/bar/cfan/genomes/", genome, "/", genome, ".fa")
  }
  if (!is.null(source.fasta)) {
    fi <- source.fasta
  }
  cmd <- paste0("/bar/cfan/anaconda2/envs/jupyter/bin/bedtools getfasta -fi ",fi,
    " -fo ", out.fa, " -bed ", root.name, ".pre.bed -name")
  print(cmd)
  system(cmd)
  
  fasta.wrap.fanc(in.fa = out.fa)
  if (return.fasta == T) {
    return(seqinr::read.fasta(out.fa))
  }
  return(out.fa)

}

# mpw <- function(from.fa=NULL, from.bed=NULL, ref.pattern, chr.name, shift=0,
#                 genome=NULL, tempt.dir, SNV.dir, aligner, dry=T, out.json=F) {
#   # note: shift can be set at equal to the left of the reference region in the bed file.
#   if (is.null(from.fa) && is.null(from.bed))
#     stop("at least one of fa or bed should be supplied")
#   if (!is.null(from.fa) && !is.null(from.bed))
#     stop("only one of fa and bed can be supplied, you supplied both somehow")
#   system(paste0("rm -rf ", tempt.dir, " ", SNV.dir))
#   system(paste0("mkdir -p ", tempt.dir, " ", SNV.dir))
#   if (!is.null(from.bed)) {
#     if (is.null(genome))
#       stop("if bed is supplied, genome must also be specified")
#     from.fa <- bed.align(bed = from.bed,  root.name = paste0(tempt.dir, "/mpw"), genome = genome, align = F)
#   }

#   # now extract reference from the fasta file:
#   fa <- seqinr::read.fasta(from.fa)
#   fa.ref <- grep(ref.pattern,sapply(fa, attr, "Annot"))
#   if (length(fa.ref) > 1)
#     stop("ref.pattern matched to more than one sequence, they are: " %>% paste0(paste0(fa.ref, collapse = ", ")))
#   if (length(fa.ref) == 0)
#     stop("ref.pattern didn't match any thing")

#   fa.ref <- fa[[fa.ref]]
#   seqinr::write.fasta(fa.ref, sub("^>", "", attr(fa.ref, "Annot")),
#                       paste0(tempt.dir, "/mpw.ref.fa"))
#   formatted.header <- sub("^>", "", sapply(fa, attr, "Annot"))  %>% gsub("[^0-9a-zA-Z]", "_",.)
#   seqinr::write.fasta(fa, formatted.header,
#                       paste0(tempt.dir, "/mpw.strain.fa"))

#   align.cmd <- paste0("/bar/cfan/anaconda2/envs/jupyter/bin/python3 ",
#                       "/bar/cfan/viralBrowser/release_github_4_22/publicAlignment.py",
#                       " --script_dir ", "/bar/cfan/viralBrowser/release_github_4_22/",
#                       " --ref_fa ", tempt.dir, "/mpw.ref.fa",
#                       " --ref_name ", chr.name,
#                       " --shift ", shift,
#                       " --strain_fa ", tempt.dir, "/mpw.strain.fa",
#                       " --tempt_dir ", tempt.dir,
#                       " --SNV_dir ", SNV.dir, 
#                       # " --email changxu.fan@gmail.com",
#                       " --aligner ", aligner)
#   print(align.cmd)
#   if (dry == F)
#     system(align.cmd)
#   if (out.json==T) {
#     # snvs <- system(paste0("ls ", SNV.dir, "/*bed.gz"), intern = T) %>% basename()
#     jsongen <- lapply(formatted.header, function(x) {
#       df <- data.frame(name = x, url = paste0(x, ".bed.gz"), type = "pairwise")
#     }) %>% Reduce(rbind,.)
#     jsongen %>% jsonlite::toJSON() %>% jsonlite::prettify() %>% write(paste0(SNV.dir, "/snv.json"))
#   }
#   names(formatted.header) <- NULL
#   return(formatted.header)
# }

# readGTF.fanc <- function(gtf.file, get.fields=NULL) {
#   gtf <- read.table(gtf, sep = "\t", quote = "", header = F, as.is = T)
#   colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
#   if (!is.null(get.fields))
#     gtf <- gtf.get.fields(gtf, fields = get.fields)
#   return(gtf)
# }

# gtf.get.fields <- function(gtf, fields) {
#   for (field in fields) {
#     gtf[, field] <- str_extract(string = gtf$attribute, pattern = paste0(field,"[^;]+")) %>% gsub("\"", "",.) %>% sub(paste0(field, " "), "",. )
#   }
#   return(gtf)
# }

# # t <- gtf.get.fields(gtf.df, "gene_name")

# stitch.fa <- function(bed, genome, out.file) {
#   # bed has to be written in a way that columns starting from the 4th one represents "group"
#   # sequences are grouped based on all columns (starting from 4th) cat together ("|" separated)
#   tempt.fa <- bed.align(bed = bed, aligned.fa = out.file, genome = genome, add.coordinate = F)

#   fas <- seqinr::read.fasta(tempt.fa)
#   fa.df <- data.frame(header = lapply(fas, attr, "name") %>% unlist() %>% `names<-`(NULL),
#                       seq = lapply(fas, function(x) paste0(x[1:length(x)], collapse = "") ) %>% unlist())
#   fa.df <- fa.df %>% group_by(header) %>% summarise(seq = paste0(seq, collapse = ""))
#   seqinr::write.fasta(sequences = fa.df$seq, names = fa.df$header, file.out = out.file)
#   return(list(fa.df = fa.df, fa.file = fa.file))

# }

# # stitch.fa("~/test/random/test.bed", genome = "mm10", out.file = "~/test/random/test_miao.fa")

# bowtie2.genome.gen <- function(map.bed, genome, bowtie2.build = "/opt/apps/bowtie2/2.3.4.1/bowtie2-build",
#                                out.dir, root.name, thread = 8, run = T) {
#   # map.bed: bed file. only first 3 columns used.currently only supports mono-region.
#   # get fasta file:
#   system(paste0("mkdir -p ", out.dir))
#   map.bed <- map.bed[1, 1:3]
#   map.bed$V4 <- map.bed[, 1]
#   fa <- bed.align(bed = map.bed, aligned.fa = paste0(out.dir, "/", root.name), genome = genome, mask = F, add.coordinate = F)
#   # now make a mini genome:
#   cmd <- paste0("cd ",out.dir," && ", bowtie2.build, " --thread ", thread,
#                 " ", basename(fa), " ", root.name)
#   print(cmd)
#   if (run == T)
#     system(cmd)

#   return(list(bowtie.genome = paste0(out.dir, "/", root.name), map.bed = map.bed, genome = genome))
# }

# mini.align <- function(fastq.vec, R1.pattern="1.fastq.gz", R2.pattern = "2.fastq.gz", k = NULL, a = F,
#                        bowtie.genome, map.bed, out.bam, bowtie.log.file=NULL, thread = 16, run = T,
#                        bowtie2 = "/opt/apps/bowtie2/2.3.4.1/bowtie2", no.unal = F) {
#   system(paste0("mkdir -p ", dirname(out.bam)))

#   R1 <- fastq.vec[grepl(R1.pattern, fastq.vec)][1]
#   R2 <- fastq.vec[grepl(R2.pattern, fastq.vec)][1]
#   # cmd <- paste0("/bar/cfan/scripts/dna-seq/bowtie2_4_22_20.sh ", " -p ", thread,
#   #               " -x ", bowtie.genome, " -i ", R1, " -I ", R2, " -o ", out.bam, " -k ", 50,
#   #               " -s ", map.bed[1,2])
#   cmd <- paste0(bowtie2, " -X2000 --reorder --very-sensitive --xeq --seed 42 ",
#                 " -p ", thread, " -x ", bowtie.genome, " -1 ", R1, " -2 ", R2)

#   # if (!is.null(bowtie.log.file))
#   #   cmd <- paste0(cmd, " --met-file ", bowtie.log.file)

#   if (!is.null(k))
#     cmd <- paste0(cmd, " -k ", k)
#   if (a == T)
#     cmd <- paste0(cmd, " -a ")
#   if (no.unal == T)
#     cmd <- paste0(cmd, " --no-unal")

#   shift <- map.bed[1,2]
#   cmd <- paste0(cmd , " | ", "awk -F \"\\t\" 'BEGIN {OFS = \"\\t\"} {$4 = $4+",shift,"; print $0}' | ")
#   cmd <- paste0(cmd, " /bar/cfan/anaconda2/envs/jupyter/bin/samtools sort -O bam -m 4G -@ ",thread," -  > ", out.bam)

#   utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run)
#   if (!file.exists(out.bam))
#     stop(paste0(out.bam, " was not successfully generated"))

#   cmd <- paste0("samtools index ", out.bam)
#   utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run, stdout.file = bowtie.log.file)
#   if (!file.exists(out.bam %>% paste0(".bai")))
#     stop(paste0(out.bam, ".bai was not successfully generated"))

#   return(out.bam)
# }

# trimmed.fastq.gen <- function(fastq.vec, thread, run=T) {
#   cmd <- paste0("python3 ", "~/software/atac/encode_pipeline/src/encode_task_trim_adapter.py ", paste0(fastq.vec, collapse = " "),
#                 " --paired-end ", " --auto-detect-adapter ", " --cutadapt-param ' -e 0.1 -m 5' ", " --nth ", thread)
#   utilsFanc::cmd.exec.fanc(cmd = cmd, run = run, intern = F)
# }
