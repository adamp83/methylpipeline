#!/usr/bin/Rscript

# All-in-one script to analyse methlation status of BiSEQ reads.
# 
# Required input files:
# - A reference genome in FASTA format in the file reference_genome/reference.fa
# - Input FASTA files for n samples. For each sample, put the read files (usually 2 .fastq files) in the directory input/<sample_name>. The script will recurse the input directory looking for .fastq files.
# 
# What this script does:
# 1. Calls runbwameth.sh script, which uses BWA open source aligner to align reads and generates SAM files using samtools.
# 2. Reads the CpG positions from the reference genome
# 3. Counts how much overlap there is between reads and defines how many bases should be trimmed from the reads based on this.
# 4. Creates a large matrix of base pairs with paired reads on the y-axis and reference positions on the x-axis
# 5. Analyses methylation of the CpG islands by comparing C to T ratios. Outputs to an Excel file.
# 6. Creates heatmaps to visualise this


##################################################
# SETTINGS
# 
# Change the following to tweak the below script.
#
##################################################

# How many bases to trim from both forward and reverse reads
# If set to NA, the script will guess based on average read overlaps
trim <- NA

##################################################
# END of settings
##################################################


# Check if the user has set a trim value
user_trim <- !is.na(trim)

# First, run bwameth script:
system("./runbwameth.sh")


# Install required packages if required
list.of.packages <- c("stringr", "WriteXLS", 'pheatmap', 'gplots', 'cluster')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(new.packages)
}

library(stringr)
library(WriteXLS)
library(pheatmap)
library("gplots")
library("cluster")

# First find positions of interest from a specified reference file

ref_file <- "./reference_genome/reference.fa"
ref <- scan(ref_file, what=character())[2]

# Pick out the CpG islands
cpgs <- str_locate_all(ref, "CG")[[1]][,1]

if(length(cpgs) == 0){
  print("ERROR: NO CpG ISLANDS IN REFERENCE GENOME")
}

outputs <- list()

# Read in SAM files
sam_dir <- "./output"
sam_files <- list.files(sam_dir, pattern=".sam", full.names=T)

for(sam_file in sam_files){
  print(paste("Analysing", sam_file))
  
  sam <- read.table(sam_file, fill=T, sep="\t")
  
  if(!user_trim){
    trim <- NA
  }
  
  # As per SAM specification:
  #   $V4 -> starting position
  #   $V5 -> mapping quality
  #   $V10 -> sequence
  
  # Make a direction flag which we will fill in
  # 1 = forward, -1 = reverse
  sam$direction = NA
  
  # Create a large matrix with read pair IDs as rows and reference positions as columns
  # Then fill in from the SAM data
  
  # Pair IDs are in column 1. Those with 3 letters are headers.
  pairs <- unique(as.character(sam[,1]))
  pairs <- pairs[which(nchar(pairs) > 7)]
  
  data <- matrix(nrow=length(pairs), ncol=nchar(ref))
  rownames(data) <- pairs
  colnames(data) <- 1:nchar(ref)
  
  # Check the average overlap distance - will help inform trimming
  # Defined as start of forward read (S1) + length of forward read (L1) - start of reverse read (S2)
  overlaps = data.frame(pair=pairs, s1=NA, l1=NA, s2=NA, overlap=NA)
  for(pair in pairs){
    #print(paste("Analysing pair", pair))
    
    reads <- sam[which(sam[,1] == pair),]
    
    
    
    reads[,4] <- as.numeric(as.character(reads[,4]))
    mapq = as.numeric(as.character(reads[,5]))
    if(!all(mapq >= 60)){
      next
    }
    read1 <- reads[which(reads[,4] == min(reads[,4])),]
    read2 <- reads[which(reads[,4] == max(reads[,4])),]
    
    # If there are duplicates, ignore this pair
    if(dim(read1)[1] > 1 | dim(read2)[1] > 1){
      print("SKIPPING READ DUE TO DUPLICATE")
      next
    }
    
    sam$direction[which(rownames(sam) == rownames(read1))] = 1
    sam$direction[which(rownames(sam) == rownames(read2))] = -1
    
    o_sel <- which(overlaps$pair == pair)
    s1 <- read1[,4]
    s2 <- read2[,4]
    l1 <- nchar(as.character(read1[,10]))
    
    overlaps[o_sel,]$s1 <- s1
    overlaps[o_sel,]$s2 <- s2
    overlaps[o_sel,]$l1 <- l1
    overlaps[o_sel,]$overlap <- s1 + l1 - s2
  }
  if(dim(overlaps)[1] == 0){
    # Don't trim at all if there aren't any overlapping bases
    if(is.na(trim)){
      trim <- 0
    }
  }else{
    overlaps = overlaps[which(!is.na(overlaps$overlap)),]
    print(paste("Read overlaps: mean", mean(overlaps$overlap), ", SD", sd(overlaps$overlap)))
    
    # Figure out a trim value, if not set
    # Base this on the mean and SD of the overlaps
    if(is.na(trim)){
      trim <- round( (mean(overlaps$overlap) - 2* sd(overlaps$overlap)) / 2 )
      trim <- max(trim, 0) # Don't allow negative trimming!
      print(paste("Trimming", trim, "bases from each read based on above values"))
    }
  }
  
  
  
  for(pair in pairs){
    reads <- sam[which(sam[,1] == pair),]
    sel = which(rownames(data) == pair)
    
    for(i in 1:dim(reads)[1]){
      start = as.numeric(as.character(reads[i,4]))
      seq = as.character(reads[i,10])
      mapq = as.numeric(as.character(reads[i,5]))
      
      # Skip if MAPQ < 60
      if(!is.na(mapq) & mapq < 60){
        next
      }
      # Make sure we know the direction
      if(is.na(reads$direction[i])){
        next
      }
      
      if(reads$direction[i] == 1){
        range <- 1:(nchar(seq) - trim)
      }else{
        if(reads$direction[i] == -1){
          range <- (1 + trim):nchar(seq)
        }
      }
      
      # Count along seq. If a value exists, remove it if it clashes, else leave it.
      for(j in range){
        pos <- start + j -1
        val = substr(seq, j, j)
        
        if(pos <= dim(data)[2] & pos > 0){
          if(is.na(data[sel, pos])){
            data[sel, pos] <- val
          }else{
            if(data[sel, pos] != val){
              data[sel, pos] <- NA
            }
          }
        }
      }
    }
  }
  
  #For each reference CpG, count the number of Cs, Ts and other
  counts <- matrix(nrow=length(cpgs), ncol=5, data=0)
  rownames(counts) <- cpgs
  colnames(counts) <- c("C", "T", "other", "total", "percentC")
  
  for(i in 1:length(cpgs)){
    cpg <- cpgs[i]
    counts[i,1] <- length(which(data[,cpg] == "C"))
    counts[i,2] <- length(which(data[,cpg] == "T"))
    counts[i,3] <- length(which(data[,cpg] == "A" | data[,cpg] == "G"))
    counts[i,4] <- length(which(!is.na(data[,cpg])))
    counts[i,5] <- 100*counts[i,1] / (counts[i,1] + counts[i,2])
  }
  counts <- data.frame(counts)
  
  # Write counts data to output file
  opfile <- gsub("sam", "methyldata.xls", sam_file)
  WriteXLS("counts", ExcelFileName=opfile, row.names=T, col.names=T)
  
  
  # Make some plots:
  # On the x-axis, the cpg sites
  # On the y-axis, the methylation status of each read
  plotdata <- matrix(nrow=dim(data)[1], ncol=length(cpgs))
  for(i in 1:dim(data)[1]){
    for(j in 1:length(cpgs)){
      cpg <- cpgs[j]
      #print(paste(i, cpg, data[i,cpg], sep="; "))
      if(is.na(data[i,cpg])){
        next
      }
      if(data[i,cpg] == "C"){
        plotdata[i,j] <- 1
      }else{
        if(data[i,cpg] == "T"){
          plotdata[i,j] <- 2
        }
      }
    }
  }
  colnames(plotdata) <- cpgs
  
  
  # Get rid of rows which have no data:
  sel <- which(apply(plotdata, 1, function(x){
    !all(is.na(x))
  }))
  if(length(sel) == 0){
    print("No coverage. Moving to next sample.")
    next
  }
  
  # Write heatmaps to a PDF file
  pdf_file = gsub(".sam", "_heatmap.pdf", sam_file)
  pdf(pdf_file)
  
  plotdata <- plotdata[sel,]
  
  # Make an unclustered heatmap first
  heatmap.2(as.matrix(plotdata),dendrogram="none",trace="none", margin=c(8,9), Rowv=F, Colv=F, key=F, labRow=F, main=paste(sam_file, "Unclustered:", dim(plotdata)[1], "reads"));
  
  
  # For clustering, need at least half to not be NA...
  sel <- which(apply(plotdata, 1, function(x){
    length(which(!is.na(x))) >= round(length(cpgs) / 2)
  }))
  plotdata <- plotdata[sel,]
  
  
  
  # Attempt to use pheatmap method - but does not respond well to NAs
  #pdf(pdf_file)
  #pheatmap(plotdata, treeheight_row=0, treeheight_col=0, cluster_cols=F, cluster_rows=F, annotation_col=collab, main=sam_file, scale="none", col=c("Red", "Green"))
  #dev.off()
  
  
  # Clustering using daisy method
  hclustfunc <- function(x) hclust(x, method="complete")
    
  # Try using daisy GOWER function 
  # which suppose to work with NA value
  distfunc <- function(x) daisy(x,metric="gower")
  
  
  plotdata <- apply(plotdata, 2, as.numeric)
  
  # For clustering we need to avoid too many of the same number in each column in large data sets
  # Therefore replace with random variation
  # See: http://stackoverflow.com/questions/16559250/error-in-heatmap-2-gplots
  plotdata2 <- plotdata
  vary <- function(n){
    if(is.na(n)){
      return(NA)
    }
    runif(1, n-0.00001, n+0.00001)
  }
  plotdata2 <- apply(plotdata2, c(1,2), vary)
  
  
  # Plot the heatmap
  d <- distfunc(plotdata2)
  fit <- hclustfunc(d)
  
  # Perform clustering heatmap
  heatmap.2(as.matrix(plotdata2),dendrogram="none",trace="none", margin=c(8,9), hclust=hclustfunc,distfun=distfunc, Colv=F, key=F, labRow=F, main=paste(sam_file, "clustered: ", dim(plotdata2)[1], "reads"));
  dev.off()
  
  
}

# Now we can load the IGV viewer
# Wait a few seconds to complete the file writing
Sys.sleep(2)
bam_files <- list.files("./output", pattern=".bam$", full.names=T)
igv <- paste("./dependencies/igv/igv.sh --genome", ref_file, paste(bam_files, collapse=","), sep=" " )
system(igv)