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
  echo "$out_prefix"
fi

#
# Stream input data
#
echo "streaming input data ..."
mkfifo ./"$reads_name"
dx cat "$reads" > ./"$reads_name" &
cat_pid=$!

#
# Set up some options
#
echo "running FastQC ..."
opts=""
if [ "$contaminants_txt" != "" ]; then
  opts="$opts -c ./in/contaminants_txt/*"
fi
if [ "$adapters_txt" != "" ]; then
  opts="$opts -a ./in/adapters_txt/*"
fi
if [ "$limits_txt" != "" ]; then
  opts="$opts -l ./in/limits_txt/*"
fi
if [ "$format" != "auto" ]; then
  opts="$opts -f $format"
fi
if [ "$nogroup" == "true" ]; then
  opts="$opts --nogroup"
fi

#
# Run FastQC
#
mkdir results
/FastQC/fastqc -t `nproc` --extract -k $kmer_size -o results $opts $extra_options ./"$reads_name"
wait "$cat_pid"

#
# Upload results
#
echo "uploading results ..."
mkdir -p ~/out/report_html/ ~/out/stats_txt/
mv results/*/fastqc_report.html ~/out/report_html/"$reads_prefix".stats-fastqc.html
mv results/*/fastqc_data.txt ~/out/stats_txt/"$reads_prefix".stats-fastqc.txt
dx-upload-all-outputs
dx mkdir -p $out_folder/Results/${out_prefix}/FASTQC/
dx mv "$reads_prefix".stats-fastqc.html $out_folder/Results/${out_prefix}/FASTQC/
dx mv "$reads_prefix".stats-fastqc.txt $out_folder/Results/${out_prefix}/FASTQC/
