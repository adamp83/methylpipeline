# methylpipeline
Simple wrapper for BWA-meth and IGV viewer for multi-sample methylation analysis

############################################
# Dependencies
############################################

Python (installed by default on Mac systems)

toolshed:
  wget https://pypi.python.org/packages/source/t/toolshed/toolshed-0.4.0.tar.gz
   tar xzvf toolshed-0.4.0.tar.gz
   cd toolshed-0.4.0
   sudo python setup.py install

samtools:
  Download stable binary from http://samtools.sourceforge.net
  You must be able to run > samtools at the command line without error
  tested with v1.2

bwa-mem:
  http://bio-bwa.sourceforge.net
  tested with version 0.7.12

bwa-meth:
  Download from https://github.com/brentp/bwa-meth/blob/master/README.md
  sudo python setup.py install
  Should be able to run bwameth.py from the command line
  Tested with master branch revision d5c6fd89eb from github
