#!/bin/bash
# stjude_macs2 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.

set -e -x -o pipefail

main() {

    #############################################
    # write log
    #############################################

    echo "Value of ChIP_bam: '$ChIP_bam'" 
    echo "Value of genome: '$genome'"
    echo "Value of out_prefix: '$out_prefix'"

    
    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    echo ""
    echo "=== Setup ==="
    sudo apt-get update
    sudo apt-get install libsigsegv2
    export PATH=/user/bin:$PATH
    gawk --version
    sudo apt-get --yes install openjdk-8-jdk
    java -version
    sudo apt-get install tabix

    echo "install spp ..."
    R --version
    R CMD INSTALL /home/bitops_1.0-6.tar.gz
    R CMD INSTALL /home/caTools_1.17.1.tar.gz
    R CMD INSTALL /home/snow_0.4-1.tar.gz
    R CMD INSTALL /home/phantompeakqualtools/spp_1.10.1.tar.gz
    echo "R CMD INSTALL /home/phantompeakqualtools/spp_1.10.1.tar.gz" >> logfile

    echo "install bedtools ..."
    sudo apt-get install bedtools
    echo "install macs2 ..."
    sudo apt-get install python-numpy 

    pip install MACS2
    echo "pip install MACS2" >> logfile

    echo "  [*] Downloading input files..." 
    dx download "$ChIP_bam" -o ${ChIP_bam_prefix}.bam
        
    # Path setup 
    export PATH=$PATH:/home/dnanexus/samtools-0.1.18
    echo "export PATH=$PATH:/home/dnanexus/samtools-0.1.18" >> logfile

    #############################################
    # parse parameters from file
    #############################################
    if [ "$parameter_file" != "" ]; then 
    echo "  [*] parsing parameters from parameter file ..."
    dx download "$parameter_file" -o parameters.txt
    source parameters.txt
    echo "$genome"
    echo "$out_prefix"
    fi
    if [ "$genome" == "mm9" ] || [ "$genome" == "mm10" ] ; then export specie=mm; fi
    if [ "$genome" == "hg19" ] || [ "$genome" == "GRCh38" ] ; then export specie=hs; fi
    if [ "$genome" == "dm3" ] ; then export specie=dm; fi


    #############################################
    # remove multiple mapped reads and duplicates
    #############################################
    echo "removing multiple mapped reads and duplicates ..."
    echo "removing multiple mapped reads and duplicates ..." >> logfile
    java -Xmx4g -jar picard.jar SortSam VALIDATION_STRINGENCY=LENIENT I=${ChIP_bam_prefix}.bam O=${ChIP_bam_prefix}.sort.bam SORT_ORDER=coordinate &
    echo "java -Xmx4g -jar picard.jar SortSam VALIDATION_STRINGENCY=LENIENT I=${ChIP_bam_prefix}.bam O=${ChIP_bam_prefix}.sort.bam SORT_ORDER=coordinate &" >> logfile
    wait
    java -Xmx4g -jar picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=${ChIP_bam_prefix}.sort.bam O=${ChIP_bam_prefix}_processed.1.bam M=marked_dup_metrics.txt REMOVE_DUPLICATES=TRUE &
    echo "java -Xmx4g -jar picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=${ChIP_bam_prefix}.sort.bam O=${ChIP_bam_prefix}_processed.1.bam M=marked_dup_metrics.txt REMOVE_DUPLICATES=TRUE &" >> logfile
    wait
    rm ${ChIP_bam_prefix}.sort.bam 

    samtools view -hb -q 1 ${ChIP_bam_prefix}_processed.1.bam  > ${ChIP_bam_prefix}_processed.bam &
    echo "samtools view -hb -q 1 ${ChIP_bam_prefix}_processed.1.bam  > ${ChIP_bam_prefix}_processed.bam &" >> logfile
    wait
    ls
    


    #############################################
    # run spp cross correlation
    #############################################
    echo "run spp cross correlation ..."
    echo "run spp cross correlation ..." >> logfile
    Rscript /home/phantompeakqualtools/run_spp_nodups.R -c=${ChIP_bam_prefix}_processed.bam -savp=${ChIP_bam_prefix}_phantomPeak.pdf -out=${ChIP_bam_prefix}_phantomPeak.out &
    echo "Rscript /home/phantompeakqualtools/run_spp_nodups.R -c=${ChIP_bam_prefix}_processed.bam -savp=${ChIP_bam_prefix}_phantomPeak.pdf -out=${ChIP_bam_prefix}_phantomPeak.out &" >> logfile
    wait

    # set up output folder
    dx mkdir -p $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/
    echo "upload cross correlation plots..."
    ChIP_cc=$(dx upload --tag sjcp-result-file --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ ${ChIP_bam_prefix}_phantomPeak.pdf --brief  )
    dx-jobutil-add-output ChIP_cc "$ChIP_cc" --class=file

    echo "parse fragment size ..."
    if [ "$fragment_size" == "NA" ];then
    fragment_size=`cat ${ChIP_bam_prefix}_phantomPeak.out | cut -f 3 | sed 's/,/\t/g' | cut -f 1`
    if [ "$fragment_size" -lt 50 ]; then dx-jobutil-report-error "The estimated fragment size for ${ChIP_bam_prefix} is smaller than 50bp based on cross correlation. please manually setup a fix fragment length and rerun the analysis." AppError; fi
    fi


    #############################################
    # run MACS2
    #############################################
    echo "Calling peaks with MACS2 ..."
    echo "Calling peaks with MACS2 ..." >> logfile
    macs2 callpeak -t ${ChIP_bam_prefix}_processed.bam -g $specie -f BAM --outdir macs2  -n ${ChIP_bam_prefix} --nomodel --extsize $fragment_size
    echo "macs2 callpeak -t ${ChIP_bam_prefix}_processed.bam -g $specie -f BAM --outdir macs2  -n ${ChIP_bam_prefix} --nomodel --extsize $fragment_size" >> logfile
    nPeaks=`wc -l macs2/${ChIP_bam_prefix}_peaks.narrowPeak | awk '{print $1}' `
    if [ "$nPeaks" == "0" ]; then dx-jobutil-report-error "No peaks are called." AppError; fi 
   

    #############################################
    # remove black list
    #############################################
    echo "rm_blackList=$rm_blackList"
    if [ "$rm_blackList" == "true" ]; then
    echo "Removing peaks overlapped with black list ..."
    black_list=${genome}-Blacklist.bed
    subtractBed -a <(gawk '{if($1 !~ /^chr/ && $1 ~ /^[1-9XYM]/ ) {print "chr"$_} else {print $_} }' macs2/${ChIP_bam_prefix}_peaks.narrowPeak) -b $black_list | cut -f 1-4 > macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.bed
    subtractBed -a <(gawk '{if($1 !~ /^chr/ && $1 ~ /^[1-9XYM]/ ) {print "chr"$_} else {print $_} }' macs2/${ChIP_bam_prefix}_summits.bed) -b $black_list | cut -f 1-4 > macs2/${ChIP_bam_prefix}_summits.clean.bed
    head macs2/${ChIP_bam_prefix}_summits.clean.bed
    tail macs2/${ChIP_bam_prefix}_summits.clean.bed
    echo "Removing peaks overlapped with black list ..." >> logfile
    echo "subtractBed -a <(gawk '{if($1 !~ /^chr/ && $1 ~ /^[1-9XYM]/ ) {print "chr"$_} else {print $_} }' macs2/${ChIP_bam_prefix}_peaks.narrowPeak) -b $black_list | cut -f 1-4 > macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.bed" >> logfile
    echo "subtractBed -a <(gawk '{if($1 !~ /^chr/ && $1 ~ /^[1-9XYM]/ ) {print "chr"$_} else {print $_} }' macs2/${ChIP_bam_prefix}_summits.bed) -b $black_list | cut -f 1-4 > macs2/${ChIP_bam_prefix}_summits.clean.bed" >> logfile
    fi

    #############################################
    # metrics
    #############################################
    echo "Summarizing metrics"
    echo "TotalReads MappedReads MapRate FinalRead DupRate" > ${ChIP_bam_prefix}.metrics.txt
    for i in "${ChIP_bam_prefix}" 
    do
        echo $i
        sambamba flagstat ${i}.bam > $i.flagstat
        TotalReads=`head -n 1 $i.flagstat | awk '{print $1}'`
        MappedReads=`head -n 5 $i.flagstat | tail -n 1 | awk '{print $1}'`
        AfterDeDup=`sambamba flagstat ${i}_processed.1.bam | head -n 5 | tail -n 1 | awk '{print $1}'`
        FinalRead=`sambamba flagstat ${i}_processed.bam | head -n 5 | tail -n 1 | awk '{print $1}'`
        MapRate=`echo "$MappedReads/$TotalReads" | bc -l `
        DupRate=`echo "1-$AfterDeDup/$MappedReads" | bc -l`
        echo "$TotalReads $MappedReads $MapRate $FinalRead $DupRate"
    done | paste - - >> ${ChIP_bam_prefix}.metrics.txt
    echo >> ${ChIP_bam_prefix}.metrics.txt
    echo "`wc -l macs2/${ChIP_bam_prefix}_peaks.narrowPeak` peaks" >> ${ChIP_bam_prefix}.metrics.txt
    if [ "$rm_blackList" == "true" ]; then
    echo "`wc -l macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.bed` peaks after subtracting peaks from blacklist" >> ${ChIP_bam_prefix}.metrics.txt
    fi
    echo "fragment size used: $fragment_size" >> ${ChIP_bam_prefix}.metrics.txt


    #############################################
    # output bigwig file
    #############################################
    echo "output bigwig file ..."
    if [ "$bw_out" == "true" ]; then
function bam2bw {
i=$1
nReads=`head -n 3 ${ChIP_bam_prefix}.metrics.txt| grep -w $i | awk '{print $5}'`
scale=`echo "15000000/$nReads" | bc -l `
bamToBed -i ${i}_processed.bam > ${i}_processed.bed
extendTag.pl ${genome}.chrom.sizes ${i}_processed.bed $fragment_size 
echo "genomeCoverageBed -bg -i ${i}_processed.extended.bed -g ${genome}.chrom.sizes -scale "$scale" | sort -k1,1 -k2,2n"
genomeCoverageBed -bg -i ${i}_processed.extended.bed -g ${genome}.chrom.sizes -scale "$scale" | sort -k1,1 -k2,2n > ${i}_processed.bedg
bedGraphToBigWig ${i}_processed.bedg ${genome}.chrom.sizes ${i}.bw
echo "genomeCoverageBed -bg -i ${i}_processed.extended.bed -g ${genome}.chrom.sizes -scale "$scale" | sort -k1,1 -k2,2n" >> logfile
echo "genomeCoverageBed -bg -i ${i}_processed.extended.bed -g ${genome}.chrom.sizes -scale "$scale" | sort -k1,1 -k2,2n > ${i}_processed.bedg" >> logfile
echo "bedGraphToBigWig ${i}_processed.bedg ${genome}.chrom.sizes ${i}.bw" >> logfile
rm ${i}_processed.bedg
}
    for prefix in "${ChIP_bam_prefix}"
    do
	bam2bw $prefix & 
    done
    wait
    fi

    # Fill in your application code here.
    #
    # To report any recognized errors in the correct format in
    # $HOME/job_error.json and exit this script, you can use the
    # dx-jobutil-report-error utility as follows:
    #
    # dx-jobutil-report-error "My error message"
    #
    # Note however that this entire bash script is executed with -e
    # when running in the cloud, so any line which returns a nonzero
    # exit code will prematurely exit the script; if no error was
    # reported in the job_error.json file, then the failure reason
    # will be AppInternalError with a generic error message.

    # The following line(s) use the dx command-line tool to upload your file
    # outputs after you have created them on the local file system.  It assumes
    # that you have used the output field name for the filename for each output,
    # but you can change that behavior to suit your needs.  Run "dx upload -h"
    # to see more options to set metadata.

    echo "uploading files ..."


    if [ "$rm_blackList" == "true" ]; then
    peak_file=$(dx upload --tag sjcp-result-file --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.bed --brief )
    LC_COLLATE=C sort -k1,1 -k2,2n macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.bed > macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.sorted.bed
    head macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.sorted.bed
    tail macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.sorted.bed
    bedToBigBed macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.sorted.bed ${genome}.chrom.sizes macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.bb
    peak_bb=$(dx upload --tag sjcp-result-file --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ macs2/${ChIP_bam_prefix}_peaks.narrowPeak.clean.bb --brief  )
    else 
    cut -f 1-4 macs2/${ChIP_bam_prefix}_peaks.narrowPeak > macs2/${ChIP_bam_prefix}_peaks.narrowPeak.bed
    echo "upload peak file ..."
    peak_file=$(dx upload --tag sjcp-result-file --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ macs2/${ChIP_bam_prefix}_peaks.narrowPeak.bed --brief  )
    LC_COLLATE=C sort -k1,1 -k2,2n macs2/${ChIP_bam_prefix}_peaks.narrowPeak.bed > macs2/${ChIP_bam_prefix}_peaks.narrowPeak.sorted.bed
    bedToBigBed macs2/${ChIP_bam_prefix}_peaks.narrowPeak.sorted.bed ${genome}.chrom.sizes macs2/${ChIP_bam_prefix}_peaks.narrowPeak.bb
    peak_bb=$(dx upload --tag sjcp-result-file --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ macs2/${ChIP_bam_prefix}_peaks.narrowPeak.bb --brief  )
    fi

    dx-jobutil-add-output peak_file "$peak_file" --class=file
    dx-jobutil-add-output peak_bb "$peak_bb" --class=file

    if [ "$bw_out" == "true" ]; then
    ChIP_bw=$(dx upload --tag sjcp-result-file --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ ${ChIP_bam_prefix}.bw --brief  )
    dx-jobutil-add-output ChIP_bw "$ChIP_bw" --class=file
    fi

    mv logfile ${ChIP_bam_prefix}.macs2.log
    macs2_log=$(dx upload --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ ${ChIP_bam_prefix}.macs2.log --brief  )
    dx-jobutil-add-output macs2_log "$macs2_log" --class=file
    metrics_file=$(dx upload --tag sjcp-result-file --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ ${ChIP_bam_prefix}.metrics.txt --brief  )
    dx-jobutil-add-output metrics_file "$metrics_file" --class=file
    peak_file_raw=$(dx upload --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ macs2/${ChIP_bam_prefix}_peaks.xls --brief  )
    dx-jobutil-add-output peak_file_raw "$peak_file_raw" --class=file
    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    # output json bed files
    source bed2jsonbed.sh
    if [ "$rm_blackList" == "true" ]; then
        bed2convert=${ChIP_bam_prefix}_peaks.narrowPeak.clean.bed
    else 
        bed2convert=${ChIP_bam_prefix}_peaks.narrowPeak.bed
    fi
    run_bed2jsonbed macs2/$bed2convert
    out_bedj=$(dx upload -r --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ macs2/${bed2convert%.bed}.bedj.gz --brief  )
    out_bedj_tbi=$(dx upload -r --path $DX_PROJECT_CONTEXT_ID:$out_folder/Results/${out_prefix}/MACS2/ macs2/${bed2convert%.bed}.bedj.gz.tbi --brief  )
    dx-jobutil-add-output bedj --class=file $out_bedj
    dx-jobutil-add-output bedj_tbi --class=file $out_bedj_tbi

}
