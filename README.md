# ChIP-Seq (St Jude)

**Background**

ChIP-seq technology are widely used to study DNA binding proteins such as transcription factors and histone modifications. The workflow we provided here calls peaks from ChIP-seq data following a similar protocol to ENCODE. First we perform quality checks over raw sequencing reads, then we map reads to a reference genome, then we preprocesses alignments by removing multiple mapped reads and PCR duplicates and perform enrichment QC. Finally, we call peaks and filter out those that overlap with ENCODE reported black lists. Predicted peak coordinates and visualization files are automatically generated for users to view using our embedded visualizers.


**Repository**

Here you'll find the source code for the applets implementing individual pipeline stages of the ChIP-seq workflow. 

  The stjude_chipseq_parameter_wrapper applet provides an interface to allow users to input all the parameters required by the downstream applets.
  
  The stjude_chipseq_fastqc applet is used to run FastQC to check the sequencing quality.

  The stjude_chipseq_bwa_backtrack_fastq_read_mapper applet is used to run BWA to map short reads onto the reference reference genome.
  
  The stjude_macs2 applet is used to run MACS2 to call narrow peaks based on reads mapped onto the genome. A control input file is required.
  
  The stjude_sicer applet is used to run SICER to call broad peaks based on reads mapped onto the genome. A control input file is required.
  
  The stjude_macs2_noControl applet is used to run MACS2 to call narrow peaks based on reads mapped onto the genome. Control file is not required.
  
  The stjude_sicer_noControl applet is used to run SICER to call broad peaks based on reads mapped onto the genome. Control file is not required.
  
  The stjude_chipseq_report applet is used to summarize output from all upstream applets and generate a summary report.
  
  The stjude_chipseq_create_viewers applet is used to generate short cut for genome browser visulization files.
  
  
**Prerequisite**  

The workflow incoporated the following softwares:

FastQC 0.11.3 [1]

BWA 0.7.15 [2]  

samtools 0.1.18 [3]

MACS2 2.1.1 [4]

SICER 1.1 [5]

SPP 1.11 [6]

These softwares are already included in the resource folder of each applet which implements them. Users don't have to download them again.

**References**  

[1]	Andrews, S., FastQC: a quality control tool for high throughput sequence data. 2010.

[2] Li, H. and R. Durbin, Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics, 2009. 25(14): p. 1754-1760.

[3]	Li, H., B. Handsaker, A. Wysoker, T. Fennell, J. Ruan, N. Homer, G. Marth, G. Abecasis, and R. Durbin, The sequence alignment/map format and SAMtools. Bioinformatics, 2009. 25(16): p. 2078-2079.

[4]	Feng, J., T. Liu, B. Qin, Y. Zhang, and X.S. Liu, Identifying ChIP-seq enrichment using MACS. Nature protocols, 2012. 7(9): p. 1728-1740.

[5]	Zang, C., D.E. Schones, C. Zeng, K. Cui, K. Zhao, and W. Peng, A clustering approach for identification of enriched domains from histone modification ChIP-Seq data. Bioinformatics, 2009. 25(15): p. 1952-1958.

[6]	Kharchenko, P.V., M.Y. Tolstorukov, and P.J. Park, Design and analysis of ChIP-seq experiments for DNA-binding proteins. Nature biotechnology, 2008. 26(12): p. 1351-1359.

