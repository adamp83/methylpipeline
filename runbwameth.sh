#!/bin/sh

# Run BWAmeth software to extract methylation data from reads.
# Dependencies are listed below.

# Ensure correct java runtime environment
export JAVA_HOME=`/usr/libexec/java_home -v '1.6*'`

############################################
# Input file definitions
############################################

# Define positions of the reference genome
# This can be any fasta file, indexing is performed below
REF=./reference_genome/reference.fa

# List samples. Each should be a folder within ./input containing 2 .fastq files
# UPDATED TO COVER ALL FILES IN input/ DIRECTORY
cd input
SAMPLES=""
for d in */ ; do
    SAMPLES="$SAMPLES ${d:0:${#d}-1}"
done
cd ..
echo $SAMPLES
# If you want to use specific samples, uncomment this:
#SAMPLES="81-24966001 91-24966011 92-24966012"

# Define bases to trim from 5' and 3' ends respectively
# Doesn't appear to affect results significantly...
TRIM=0,0

# Define required mapping quality.
# Does not produce .bed files if set to >60.
MAPQ=60

# The BisSNP java package is a SNP caller designed for bisulfite methylation
# See http://doi.org/10.1186/gb-2012-13-7-r61
# This should point to the .jar file
BISSNP=./dependencies/BisSNP-0.82.2.jar

############################################
# End of input
############################################


echo "samples"
echo ${SAMPLES[*]}

# Loop over the given samples
for SAMPLE in $SAMPLES

do
  # Print which sample is being processes
	echo "Processing $SAMPLE"

  # Define the input directory and select all .fastq files in that directory
  IPDIR=./input/$SAMPLE
  echo $IPDIR
  FILES=$IPDIR/*.fastq
  echo $FILES

  # Call the bwameth and samtools indexing functions on the reference genome
  bwameth.py index $REF
  samtools faidx $REF

  # Align files -> produces .bam file
  bwameth.py --reference $REF -p output/$SAMPLE $FILES

  # Tabulate methylation -> produces .bed file
  bwameth.py tabulate --trim $TRIM --map-q $MAPQ --reference $REF --bissnp $BISSNP --prefix output/$SAMPLE -t 12 output/$SAMPLE.bam

done
wait

# Open any generated .bed files
open -e ./output/*.meth.bed

# Launch the IGV viewer to inspect BAM files
# (comment this if you don't want to inspect the alignments)
# Input should be a comma-delimited list of BAM files - convert to output/$SAMPLE.bam
# (all of the next 5 lines are just converting the $SAMPLES string into a list of input bam files)
ARR=($SAMPLES)
ARR=( "${ARR[@]/%/.bam}" )
ARR=( "${ARR[@]/#/output/}" )
SAMPLES_STR=$(printf ",%s" "${ARR[@]}")
SAMPLES_STR=${SAMPLES_STR:1}

echo "Loading BAM viewer"
./dependencies/igv/igv.sh --genome $REF $SAMPLES_STR



############################################
# Dependencies
############################################
#
# Python (installed by default on Mac systems)
#
# toolshed:
#   wget https://pypi.python.org/packages/source/t/toolshed/toolshed-0.4.0.tar.gz
#    tar xzvf toolshed-0.4.0.tar.gz
#    cd toolshed-0.4.0
#    sudo python setup.py install
#
# samtools:
#   Download stable binary from http://samtools.sourceforge.net
#   You must be able to run > samtools at the command line without error
#   tested with v1.2
#
# bwa-mem:
#   http://bio-bwa.sourceforge.net
#   tested with version 0.7.12
#
# bwa-meth:
#   Download from https://github.com/brentp/bwa-meth/blob/master/README.md
#   sudo python setup.py install
#   Should be able to run bwameth.py from the command line
#   Tested with master branch revision d5c6fd89eb from github
