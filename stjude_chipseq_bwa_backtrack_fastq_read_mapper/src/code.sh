#!/bin/bash
#

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

# read parameters from file
if [ "$parameter_file" != "" ]; then 
echo "  [*] parsing parameters from parameter file ..."
dx download "$parameter_file" -o parameters.txt
source parameters.txt
echo "$genome"
echo "$out_prefix"
#if [ "$genome" == "mm9" ]; then genomeindex_targz=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-BQbb7g80f9Vp7JJVX0vQ9BKb ; fi
#if [ "$genome" == "mm10" ]; then genomeindex_targz=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-BQbb8z00FKKB7JJVX0vQ9Bb5 ; fi
#if [ "$genome" == "hg19" ]; then genomeindex_targz=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-B6qq53v2J35Qyg04XxG0000V ; fi
#if [ "$genome" == "GRCh38" ]; then genomeindex_targz=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-BFBy4G805pXZKqV1ZVGQ0FG8 ; fi
#if [ "$genome" == "dm3" ]; then genomeindex_targz=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-BQbYZv009Byk9g9V1Vf004jJ ; fi
#if [ "$genome" == "ce10" ]; then genomeindex_targz=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-BQbYJpQ09j3x9Fj30kf003JG ; fi
if [ "$genome" == "hg19" ]; then dx download -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/Homo_sapiens/GRCh37-lite/BWA ; fi
if [ "$genome" == "GRCh38" ]; then dx download -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/Homo_sapiens/GRCh38_no_alt/BWA ; fi
if [ "$genome" == "mm9" ]; then dx download -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/Mus_musculus/MGSCv37/BWA/ ; fi
if [ "$genome" == "mm10" ]; then dx download -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/Mus_musculus/GRCm38/BWA ; fi
if [ "$genome" == "dm3" ]; then dx download -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/Drosophila_melanogaster/BDGPr5/BWA ; fi
fi

#
# Fetch and uncompress genome
#
genome_file=`ls BWA/*.bwt`     # Locate a file called <ref>.bwt
genome_file="${genome_file%.bwt}" # Remove the bwt suffix to keep the <ref>

#
# Fetch reads
#
dx-download-all-inputs --except genomeindex_targz --parallel

#
# Set up options for "bwa aln"
#
opts=""
if [ "$phred64" == "true" ]; then
  opts="$opts -I"
fi
if [ "$skip_casava_filtered" == "true" ]; then
  opts="$opts -Y"
fi
if [ "$advanced_aln_options" != "" ]; then
  opts="$advanced_aln_options"
fi

#
# Run bwa aln
#
bwa aln -t `nproc` -f out.sai "$genome_file" ./in/reads_fastqgz/* $opts
echo "bwa aln -t `nproc` -f out.sai "$genome_file" ./in/reads_fastqgz/* $opts" >> logfile
if [ "$reads2_fastqgz" != "" ]; then
  bwa aln -t `nproc` -f out2.sai "$genome_file" ./in/reads2_fastqgz/* $opts
  echo "bwa aln -t `nproc` -f out2.sai "$genome_file" ./in/reads2_fastqgz/* $opts" >> logfile
fi

#
# Set up options for "bwa samse" or "bwa sampe"
#
opts=""
if [ "$add_read_group" == "true" ]; then
  opts="$opts -r @RG\\tID:${read_group_id}\\tPL:${read_group_platform}\\tPU:${read_group_platform_unit}\\tLB:${read_group_library}\\tSM:${read_group_sample}"
fi
if [ "$preload_index" == "true" -a "$reads2_fastqgz" != "" ]; then
  opts="$opts -P"
fi
if [ "$advanced_samse_options" != "" -a "$reads2_fastqgz" == "" ]; then
  opts="$advanced_samse_options"
fi
if [ "$advanced_sampe_options" != "" -a "$reads2_fastqgz" != "" ]; then
  opts="$advanced_sampe_options"
fi

#
# Run "bwa samse" or "bwa sampe", and samtools
#
if [ "$reads2_fastqgz" == "" ]; then
  bwa samse $genome_file out.sai ./in/reads_fastqgz/* $opts | samtools view -u -S - | samtools sort -@ `nproc` - output
  echo "bwa samse $genome_file out.sai ./in/reads_fastqgz/* $opts | samtools view -u -S - | samtools sort -@ `nproc` - output" >> logfile 
else
  bwa sampe $genome_file out.sai out2.sai ./in/reads_fastqgz/* ./in/reads2_fastqgz/* $opts | samtools view -u -S - | samtools sort -m 256M -@ `nproc` - output
  echo "bwa sampe $genome_file out.sai out2.sai ./in/reads_fastqgz/* ./in/reads2_fastqgz/* $opts | samtools view -u -S - | samtools sort -m 256M -@ `nproc` - output" >> logfile 
fi
#############################################
# add "chr" to chromosme name
#############################################
samtools view -H output.bam | sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' | sed -e 's/SN:MT/SN:chrM/'|sed -e 's/SN:M/SN:chrM/' | sed -e 's/SN:U/SN:chrU/' | samtools reheader - output.bam > output.reheader.bam
mv output.reheader.bam output.bam

samtools index output.bam

# output software version to file
echo "bwa Version: 0.7.12-r1039" >> logfile 
echo "samtools Version: 1.1 (using htslib 1.1)" >> logfile

#
# Upload result
#
name="$reads_fastqgz_prefix"
# Remove any _R1 / _1 suffixes
name="${name%_1}"
name="${name%_R1}"
mkdir -p ~/out/sorted_bam/ ~/out/sorted_bai/ ~/out/bwa_log/
mv output.bam ~/out/sorted_bam/"$name".bam
mv output.bam.bai ~/out/sorted_bai/"$name".bam.bai
mv logfile ~/out/bwa_log/"$name".bwa.log
echo "before upload"
dx-upload-all-outputs --parallel
echo "uploaded"
dx mkdir -p $out_folder/Results/${out_prefix}/BWA/
dx mv "$name".bam $out_folder/Results/${out_prefix}/BWA/
dx mv "$name".bam.bai $out_folder/Results/${out_prefix}/BWA/
dx mv "$name".bwa.log $out_folder/Results/${out_prefix}/BWA/
