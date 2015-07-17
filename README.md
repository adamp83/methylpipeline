# methylpipeline
Simple wrapper for BWA-meth and IGV viewer for multi-sample methylation analysis

# Quick Start

1. Install dependencies below
2. Create a folder 'input' in the same folder as this script. This should contain n subfolders, with each subfolder containing FASTA files for analysis (typically 2 paired samples). These must have a .fastq extension.
3. Create a folder 'reference_genome'. Put a FASTA file here to be used as a reference - this should be called reference.fa by default.
4. Run the runbwameth.sh script.

Output will be:
* Alignments will be loaded in IGV viewer
* Methlation summaries (as .bed files) will be opened in a text editor.

############################################
# Dependencies
############################################

Python (installed by default on Mac systems)

toolshed:

* wget https://pypi.python.org/packages/source/t/toolshed/toolshed-0.4.0.tar.gz
* tar xzvf toolshed-0.4.0.tar.gz
* cd toolshed-0.4.0
* sudo python setup.py install

samtools:
* Download stable binary from http://samtools.sourceforge.net
* You must be able to run > samtools at the command line without error
* tested with v1.2

bwa-mem:
* http://bio-bwa.sourceforge.net
* tested with version 0.7.12

bwa-meth:
* Download from https://github.com/brentp/bwa-meth/blob/master/README.md
* sudo python setup.py install
* Should be able to run bwameth.py from the command line
* Tested with master branch revision d5c6fd89eb from github

Copies of BisSNP and IGV binaries are included in this repository. Further information can be found here:
* BisSNP: http://people.csail.mit.edu/dnaase/bissnp2011/
* IGV: https://www.broadinstitute.org/igv/
