#!/bin/bash
#

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

apt-get update
apt-get -y install \
   bwa \
   samtools

# read parameters from file
if [ "$parameter_file" != "" ]; then 
  echo "  [*] parsing parameters from parameter file ..."
  dx download "$parameter_file" -o parameters.txt
  source parameters.txt
  echo "$genome"
  echo "$out_prefix"
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
  echo "bwa samse $genome_file out.sai ./in/reads_fastqgz/* $opts | samtools view -u -S - | samtools sort -@ `nproc` - > output.bam" >> logfile
  bwa samse $genome_file out.sai ./in/reads_fastqgz/* $opts | samtools view -u -S - | samtools sort -@ `nproc` - > output.bam
else
  echo "bwa sampe $genome_file out.sai out2.sai ./in/reads_fastqgz/* ./in/reads2_fastqgz/* $opts | samtools view -u -S - | samtools sort -m 256M -@ `nproc` - > output.bam" >> logfile
  bwa sampe $genome_file out.sai out2.sai ./in/reads_fastqgz/* ./in/reads2_fastqgz/* $opts | samtools view -u -S - | samtools sort -m 256M -@ `nproc` - > output.bam
fi
#############################################
# add "chr" to chromosme name
#############################################
samtools view -H output.bam | sed -re 's/SN:([0-9,X,Y,M,T,U,E]+)/SN:chr\1/' | sed -e 's/SN:chrMT/SN:chrM/' | samtools reheader - output.bam > output.reheader.bam
mv output.reheader.bam output.bam

samtools index output.bam

# output software version to file
set +e
bwa 2>> logfile
set -e
samtools --version >> logfile

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
