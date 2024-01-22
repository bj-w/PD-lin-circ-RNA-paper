library(data.table) # fast reading and writing of files

# set path of BSJ tools outputs
output_dir <- "path/to/outputdir"

# functions ---------------------------------------------------------------
successfulSamples <- function(filesPresent) {
  # remove file suffix
  filesPresent$file <- sub("\\.bed", "", filesPresent$file)
  # put files into a table
  filesPresent <- data.frame(table(filesPresent))
  # sum number of occurences by tool
  filesPresent <- aggregate(Freq ~ file, filesPresent, sum)
  # remove samples not succesfully run through all tools
  success <- subset(filesPresent, Freq == max(filesPresent$Freq))
  # pull out sampled IDs
  success <- as.character(success$file)
  return(success)
}

importAndFormatToolOutput <- function(tool) {
  samples <- paste0(output_dir, "/", tool, "/", usable_samples, ".bed")
  print(paste0("Importing and formatting counts from ", tool))
  if (tool == "ciri") {
    print("Reading in files...")
    counts <- lapply(samples, fread, header = TRUE)
    names(counts) <- usable_samples
    print("Formatting")
    counts <- lapply(counts, function(x) {
      x$circRNA_start <- x$circRNA_start - 1
      reformat <- data.frame(coord_id = paste0(x$chr, ":", x$circRNA_start, "-", x$circRNA_end, ":", x$strand), count = x[, 5])
      return(reformat)
    })
  } else if (tool == "ce"){
    print("Reading in files...")
    counts <- lapply(samples, fread, header = FALSE)
    names(counts) <- usable_samples
    print("Formatting")
    counts <- lapply(counts, function(x) {
      reformat <- data.frame(coord_id = paste0(x$V1, ":", x$V2, "-", x$V3, ":", x$V6), count = x$V13)
      return(reformat)
    })
  } else if (tool == "pf") {
    print("Reading in files...")
    counts <- lapply(samples, fread, header = FALSE)
    names(counts) <- usable_samples
    print("Formatting")
    counts <- lapply(counts, function(x) {
      reformat <- data.frame(coord_id = paste0(x$V1, ":", x$V2, "-", x$V3, ":", x$V6), count = x$V5)
      return(reformat)
    })}
    else {
      print("Tool should be one of: ciri, ce or pf")
    }
  # add name of sample as a new column
  counts <- Map(cbind, counts, id = names(counts))
  # merge into one df
  print("Merging into one df")
  counts <- do.call(rbind, counts)
  # rbind adds rownames so remove these
  rownames(counts) <- NULL
  # pivot wider
  print("Pivot into wider format")
  countsWide <- reshape(counts, direction = "wide", idvar = "coord_id", timevar = "id")
  # reshape adds annoying prefix so remove it
  names(countsWide) <- gsub("count.", "", names(countsWide))
  # also remove X.junction_reads. from ciri table
  names(countsWide) <- gsub("X.junction_reads.", "", names(countsWide))
  # convert NAs to 0
  countsWide[is.na(countsWide)] <- 0
  # keep BSJs with at least two reads in two or more samples
  print("Keeping BSJs present in at least 2 samples..")
  keep <- rowSums(countsWide[, -1] >= 2) >= 2
  countsWide <- countsWide[keep, ]
  return(countsWide)
}

# which samples have been successfully run through all tools? --------------
files <- read.table(text = list.files(path = output_dir, recursive = TRUE, full.names = FALSE), sep = "/", col.names = c("tool", "file"))
usable_samples <- successfulSamples(files)

# import and reformat counts ----------------------------------------------
list_of_tools <- list.dirs(path = output_dir, full.names = FALSE, recursive = FALSE)
combined_bsj <- lapply(list_of_tools, importAndFormatToolOutput)
names(combined_bsj) <- list_of_tools

# find high confidence BSJs -----------------------------------------------
# add on tool name as column
high_conf_bsj <- Map(cbind, combined_bsj, tool = names(combined_bsj))
# how many times was each BSJ detected by each tool?
high_conf_bsj <- as.data.frame.matrix(table(do.call(rbind, lapply(high_conf_bsj, function(x) {
  x[, c("coord_id", "tool")]
}))))
# keep BSJs detected by at least 2 tools
print("Present in at least 2 tools?")
print(table(rowSums(high_conf_bsj >= 1) >= 2))
high_conf_bsj <- high_conf_bsj[rowSums(high_conf_bsj >= 1) >= 2, ]
high_conf_bsj <- rownames(high_conf_bsj)
# make BED file of these high confidence BSJs
high_conf_bsj_bed <- read.table(text = high_conf_bsj, sep = ":", col.names = c("chr", "coords", "strand"))
high_conf_bsj_bed <- cbind(high_conf_bsj_bed, read.table(text = high_conf_bsj_bed$coords, sep = "-", col.names = c("start", "end")))
high_conf_bsj_bed <- data.frame(high_conf_bsj_bed, name = high_conf_bsj)
high_conf_bsj_bed$score <- "ignore"
high_conf_bsj_bed <- unique(high_conf_bsj_bed[, c("chr", "start", "end", "name", "score", "strand")])
# format for ciriquant
ciriquant_bed <- data.frame(chr = high_conf_bsj_bed$chr, 
           start = as.numeric(high_conf_bsj_bed$start) + 1,
           end = high_conf_bsj_bed$end, 
           coord = paste0(high_conf_bsj_bed$chr, ":", high_conf_bsj_bed$start, "|", high_conf_bsj_bed$end),
           score = ".",
           strand = high_conf_bsj_bed$strand)
write.table(high_conf_bsj_bed, "highConfidenceBSJ.bed", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(ciriquant_bed, "ciriquant.bed", quote = F, sep = "\t", row.names = F, col.names = F)

# filter BSJ count matrices for high confidence junctions -----------------
combined_bsj <- lapply(combined_bsj, function(x) {
  subset(x, coord_id %in% high_conf_bsj)
})

# export count BSJ count tables -------------------------------------------
sapply(names(combined_bsj), function(x) {
  fwrite(combined_bsj[[x]], file = paste0(x, "_BSJ_aggregate.csv"))
})