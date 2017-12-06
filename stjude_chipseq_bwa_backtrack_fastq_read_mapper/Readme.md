# BWA-backtrack FASTQ Read Mapper

## What does this app do?

This app maps FASTQ reads (paired or unpaired) to a reference genome with the BWA-backtrack algorithm.

## What are typical use cases for this app?

This app can be used to map reads to a reference genome, which is a typical step in most bioinformatics
analyses. It is suitable as a first step if you are planning on performing downstream variation calling with GATK,
or any other analysis which requires mappings.

This app uses the BWA (Burrows-Wheeler Alignment Tool) software package. BWA includes three algorithms (BWA-backtrack,
BWA-SW and BWA-MEM), but this app is specifically using the BWA-backtrack algorithm (commonly referred to as "bwa aln").

The BWA-backtrack algorithm works well for for Illumina sequence reads up to 100bp. If you have longer reads, consider
using the BWA-MEM app instead.

## What data are required for this app to run?

This app requires reads files in gzipped FASTQ format (`*.fastq.gz` or `*.fq.gz`), such as those typically produced by Illumina
instruments. A single file is needed for unpaired reads, and two files (left and right read mates) are required for paired reads.

The app also requires a reference genome sequence index. This must be a gzipped tar archive file (`*.bwa-index.tar.gz`) containing
all the sequence index files as previously output by the BWA indexer. (Indexing is an one-time operation that needs to be performed to a
reference genome sequence in order for it to be usable by BWA. If you have created your own BWA index outside of DNAnexus,
place all the index files in a gzipped tar archive and provide that as the input; if you have a reference genome sequence in FASTA
format, you can index it with the BWA FASTA Indexer app). Some pre-indexed genomes are also available as suggested inputs. See
also ['which human reference sequence should I use?'](https://answers.dnanexus.com/p/183/) on DNAnexus Answers.

If you will be using the mappings to perform variation calling, we encourage you to provide the correct read group information,
and in particular to enter a read group sample. (The default behavior is to add a read group named after the input file, assign
the value `ILLUMINA` to the read group platform, and the value `1` to the read group sample).

## What does this app output?

This app outputs the mappings, as a coordinate-sorted BAM file (`*.bam`). An associated BAM index file (`*.bai`) is also generated.

## How does this app work?

This app uses the BWA-backtrack algorithm from the BWA software package. For general information, consult the BWA manual at:

http://bio-bwa.sourceforge.net/bwa.shtml

This app performs the following steps:

- Mapping with BWA-backtrack (for unpaired reads: `bwa aln` and `bwa samse`; for paired reads: `bwa aln` and `bwa sampe`).
- Conversion to BAM with `samtools view`.
- Sorting with `samtools sort`.
- Indexing with `samtools index`.

