# Analyse BAM files

library(stringr)

# First find positions of interest from a specified reference file

ref_file <- "./bwameth-method/reference_genome/reference.fa"
ref <- scan(ref_file, what=character())[2]

# Pick out the CpG islands
cpgs <- str_locate_all(ref, "CG")[[1]][,1]

outputs <- list()

# Read in SAM files
sam_dir <- "./bwameth-method/output"
sam_files <- list.files(sam_dir, pattern=".sam", full.names=T)

for(sam_file in sam_files){
  
  sam <- read.table(sam_file, fill=T, sep="\t")
  
  # As per SAM specification:
  #   $V4 -> starting position
  #   $V5 -> mapping quality
  #   $V10 -> sequence
  
  # For each reference CpG, count the number of Cs, Ts and other
  counts <- matrix(nrow=length(cpgs), ncol=4, data=0)
  rownames(counts) <- cpgs
  colnames(counts) <- c("C", "T", "other", "total")
  
  for(i in 1:dim(sam)[1]){
    
    #print(paste("Checking row", i))
    start = as.numeric(as.character(sam[i,4]))
    seq = as.character(sam[i,10])
    mapq = as.numeric(as.character(sam[i,5]))
    
    # Skip header rows
    if(is.na(seq) | seq == ""){
      next
    }
    # Skip if MAPQ < 60
    if(mapq < 60){
      next
    }
    
    # Find out if this read is paired
    paired <- F
    if(duplicated(sam[,1])[i]){
      paired <- T
      # Find the paired line
      # This will always return 2 rows, the second being the current row. Select the first.
      pair <- sam[which(sam[,1] == sam[i,1]),][1,]
    }
    
    for(n in 1:length(cpgs)){
      pos <- cpgs[n]
      #print(paste("Checking position ", pos))
      index <- pos - start + 1
      skip = F
      if(index < 1){
        #print("Index < 0")
        skip = T
      }
      if(index > nchar(seq)){
        #print("Index > seq length")
        skip = T
      }
      
      if(!skip){
        
        read <- substr(seq, index, index)
        counts[n, 4] = counts[n, 4] + 1 # Total count
        
        if(read == "C"){
          counts[n, 1] = counts[n, 1] + 1
        } else{
          if(read == "T"){
            counts[n, 2] = counts[n, 2] + 1
          } else{
            counts[n, 3] = counts[n, 3] + 1
          }
        }
      }
    }
    
  }
  
  # Post-processing - identify percentage of Cs
  df <- data.frame(counts)
  df$meth <- 100*df$C / (df$C + df$T)
 
  outputs[length(outputs) + 1] <- df
  
}