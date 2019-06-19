# contfilter
Utility for removing likely contaminant alignments from a BAM file using BAM files mapping the same set of reads to suspected contaminant genomes.

    usage: contfilter [options] cont1.bam cont2.bam
      -edit-penalty float
        	multiple for how to penalize edit distance (default 2)
      -ercc
        	exclude ERCC mappings from sample before filtering
      -limit int
        	limit the number of sample reads considered (0 = no limit)
      -log string
        	write parameters and stats to a log file
      -margin float
        	how much better sample needs to be matched (default 1)
      -max-edit-dist int
        	max edit distance for a sample match (default 5)
      -min-len int
        	min length for an alignment (default 60)
      -output string
        	output bam file (required)
      -sample string
        	BAM file of the sample you want to filter (sorted by name, required)
      -verbose
        	keep a record of what happens to each read in the log (must give -log name)
