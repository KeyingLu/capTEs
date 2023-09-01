args <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
srcdir <- dirname(script.name)
softwaredir <- dirname(srcdir)
library('rtracklayer')
library('biomartr')
library('dplyr')
library('data.table')
library('Rsamtools')
library('reshape2')
library('scales')
library('plyr')

ID = args[6]
gtfFile = args[7]


# RM + uniquely mapping res
readAnno_func <- function(ID, RM_file, bamFile, gtfFile, div=18, w=50, sw_score=225){
  # w : the region length of head or tail
  # loading gtf file
  gtf <- rtracklayer::import(gtfFile)
  gtf <- as.data.frame(gtf)
  transcript_gtf <- gtf[gtf$type == "transcript",]
  coding_gene =  gtf[gtf$type=="gene" & gtf$gene_type=="protein_coding", ]
  
  # bam file
  bam <- scanBam(bamFile)
  bam.tb <- data.frame(bam[[1]])
  bam.tb <- bam.tb[bam.tb$flag %in% c(0, 16), c("qname", "rname")]
  colnames(bam.tb) <- c("qry_id", "transcript_id")
  bam.anno <- merge(bam.tb, transcript_gtf[, c("transcript_id", "gene_id", "gene_name", "gene_type")], by="transcript_id")
  rownames(bam.anno) <- bam.anno$qry_id
  
  # RepeatMasker results
  RM.tb <- read_rm(RM_file)
  RM.tb <- data.frame(RM.tb)
  RM.tb <- RM.tb[as.numeric(RM.tb$sw_score)>=sw_score & as.numeric(RM.tb$perc_div)<=div, ]
  RM.tb <- RM.tb[grepl("L1HS", RM.tb$repeat_id) | 
                   grepl("L1P", RM.tb$repeat_id) | 
                   RM.tb$matching_class %in% c("SINE/Alu"), ]

  #
  TEs <- c("SINE/Alu", "LINE/L1")
  readTE.list <- lapply(setNames(TEs, TEs), function(TE){
    
    sub <- RM.tb[RM.tb$matching_class %in% TE, ]
    label <- apply(sub, 1, function(x){
      start <- as.integer(x[6])
      end <- as.integer(x[7])
      left <- x[8]
      left <- gsub("\\)", "", gsub("\\(", "", left))
      left <- as.integer(left)
      if(start <= w & left <= w){
        label0 = "cross"
      }else if(start <= w & left > w){
        label0 = "head"
      }else if(start > w & left <= w){
        label0 = "tail"
      }else{
        label0 = "middle"
      }
      return(label0)
    })
    sub$label <- label
    sub.list  <- split(sub, factor(sub$qry_id))
    
    a <- lapply(sub.list, function(df){
      read_name <- df[1,5]
      p <- paste(df$qry_start, df$qry_end, df$qry_left,sep=",")
      rep  <- paste(df$matching_class, df$repeat_id, sep=",")
      label <- df$label
      
      p <- paste(p, collapse = ":")
      rep <- paste(rep, collapse = ":")
      label <- paste(label, collapse = ":")
      
      
      #
      read_info <- data.frame(read_name, p, rep, label, dim(df)[1])
      colnames(read_info) <- c("read_name", "position", "repeat", "label", "number")
      
      return(read_info)
    })
    a <- rbindlist(a)
    a <- data.frame(a)
    rownames(a) <- a$read_name
    # a <- merge(a, talon.tb, by="row.names", all.x=TRUE)
    a <- merge(a, bam.anno, by="row.names", all.x=TRUE)
    a$gene_novelty <- a$gene_id
    a$gene_novelty[!is.na(a$gene_id)] <- "Known"
    a$protein_coding[a$gene_type =="protein_coding"] <- "YES"
    
    return(a)
  })
  
  
  return(readTE.list)
}





Reads_head_PCT_func <- function(ID, Test_anno){
  # total reads
  pa <- system(sprintf("grep '>' %s.trimming.Q7.L300.fa | wc -l", ID), intern=T)
  pa <- as.integer(tstrsplit(pa, " ")[[1]])
  
  ## Alu
  tb <- Test_anno[["SINE/Alu"]]
  Alu_target <- tb[grepl("head",tb$label), "read_name"]
  
  ## L1
  tb <- Test_anno[["LINE/L1"]]
  L1_target <- tb[grepl("head",tb$label), "read_name"]
  
  # 
  Alu_tartget_reads = length(setdiff(Alu_target, L1_target))
  L1_tartget_reads = length(setdiff(L1_target, Alu_target))
  sgRNA_num <- c("Alu-target Sequence"=Alu_tartget_reads, "L1-target Sequence"=L1_tartget_reads, 
                 "Alu&L1-target Sequence"=length(intersect(Alu_target, L1_target)), 
                 "Not Found"=pa-length(union(Alu_target, L1_target)))
  
  return(sgRNA_num)
  
}


Reads_PCT_func <- function(ID, Test_anno){
    # total reads
    pa <- system(sprintf("grep '>' %s.trimming.Q7.L300.fa | wc -l", ID), intern=T)
    pa <- as.integer(tstrsplit(pa, " ")[[1]])
  
    ## Alu
    tb <- Test_anno[["SINE/Alu"]]
    Alu_target <- tb$read_name

    ## L1
    tb <- Test_anno[["LINE/L1"]]
    L1_target <- tb$read_name
  
    # pie
    Alu_tartget_reads = length(setdiff(Alu_target, L1_target))
    L1_tartget_reads = length(setdiff(L1_target, Alu_target))
    sgRNA_num <- c("Alu-target Sequence"=Alu_tartget_reads, "L1-target Sequence"=L1_tartget_reads, 
                 "Alu&L1-target Sequence"=length(intersect(Alu_target, L1_target)), 
                 "Not Found"=pa-length(union(Alu_target, L1_target)))
      
    return(sgRNA_num)
  
}


####
anno <- readAnno_func(ID, 
                    RM_file=sprintf("%s.trimming.Q7.L300.fa.out", ID), 
                    bamFile <- sprintf("%s.transcriptom.sort.unique.bam", ID), 
                    gtfFile=gtfFile,
                    div=18, w=50, sw_score=225)

a <- Reads_head_PCT_func(ID, anno)
aa <- data.frame(Facet="TE in head (<50bp)",
                 Count=c(sum(a)-a["Not Found"], a["Not Found"]),
                 Class=c("TE in head (<50bp)", "Not Found")
                 )
b <- Reads_PCT_func(ID, anno)
bb <- data.frame(Facet="with TE",
                 Count=c(sum(b)-b["Not Found"], b["Not Found"]),
                 Class=c("with TE", "Not Found")
                 )
    
c <- data.frame(rbind(aa, bb), ID=ID) 
write.table(c, 'target_efficiency.txt', sep='\t', quote=F, row.names=F, col.names=T)

