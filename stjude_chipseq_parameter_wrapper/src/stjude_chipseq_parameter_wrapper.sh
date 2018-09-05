#!/bin/bash

set -e -x -o pipefail

main() {
    if [ "$genome"  == "dm3" ]; then rm_blackList="false"; fi
    out_folder=$(jq --raw-output ".folder" <<< $(dx describe --json $DX_JOB_ID))

    echo "Value of out_prefix: '$out_prefix'"
    out_prefix=${out_prefix// /_}
    echo "Sanitized value of out_prefix: '$out_prefix'" 
    echo "Value of genome: '$genome'"
    echo "Value of bw_out: '$bw_out'"
    echo "Value of rm_blackList: '$rm_blackList'"
    echo "Value of fragment_size: '$fragment_size'"
    echo "Value of out_folder: '$out_folder'"
    out_folder=${out_folder// /_}
    echo "Sanitized value of out_folder: '$out_folder'"

    # check input data
    if [ "$ChIP_fastq" == "$Control_fastq" ]; then
	dx-jobutil-report-error "ChIP fastq file is the same with Control fastq file" AppError 
    fi
    echo "Input fastq tested."    

    unset DX_WORKSPACE_ID
    dx cd $DX_PROJECT_CONTEXT_ID:
    new_file_sig=0
    for i in `echo "$out_folder/Results/$out_prefix/" | perl -e '$line=<STDIN>; print join("\n", split(/\//, "$line")), "\n";' | grep -vP "^$"`
    do
	dx pwd
	dx ls
	has_file=0;
	for j in `dx ls`
	do
	    j=$(echo $j|sed 's/\/$//')
	    if [ "$j" == "$i" ]; then has_file=1; fi
	done
	if [ "$has_file" == "1" ]; then 
	    dx cd $i; 
	else 
	    echo "$i is new"; 
	    new_file_sig=1; 
	    break; 
	fi 
    done
    echo "Output directory tested: $new_file_sig"
    
    dx cd $DX_PROJECT_CONTEXT_ID:/
    if [ "$new_file_sig" == "0" ]; then
        dCont=$(dx ls "$DX_PROJECT_CONTEXT_ID":"$out_folder/Results/$out_prefix/") 
        echo "Number of files in this folder: ${#dCont}"
        if [ "${#dCont}" -gt 0 ];then
            # Make output path unique by adding the Job ID to the out_prefix
            # The pipeline needs an empty directory downstream
            out_prefix="${out_prefix}_${DX_JOB_ID:0:10}"
        fi
    fi

    dx mkdir -p "$DX_PROJECT_CONTEXT_ID":"$out_folder/Results/$out_prefix/"

    # Writing the parameters file
    echo "export out_folder=$out_folder" > $out_prefix.parameters.txt
    echo "export out_prefix=$out_prefix" >> $out_prefix.parameters.txt
    echo "export genome=$genome" >> $out_prefix.parameters.txt
    echo "export bw_out=$bw_out" >> $out_prefix.parameters.txt
    echo "export rm_blackList=$rm_blackList" >> $out_prefix.parameters.txt
    echo "export fragment_size=$fragment_size" >> $out_prefix.parameters.txt
    echo "export ChIP_fastq_name=$ChIP_fastq_prefix" >> $out_prefix.parameters.txt
    echo "export Control_fastq_name=$Control_fastq_prefix" >> $out_prefix.parameters.txt

    parameter_file=$(dx upload --path $out_folder/Results/$out_prefix/ $out_prefix.parameters.txt --brief)
#    mkdir -p $out_folder/Results/$out_prefix/
#    mv $out_prefix.parameters.txt $out_folder/Results/$out_prefix/
#    parameter_file=$(dx upload -r $out_folder/Results/$out_prefix/$out_prefix.parameters.txt --brief)

    dx-jobutil-add-output parameter_file "$parameter_file" --class=file
}
