#!/bin/bash
set -x -e

main() {

    echo "Value of view_files: '${view_files[@]}'"
    echo "Value of parameter_file: '$parameter_file'"
    echo "Value of genome: '$genome'"

    sudo apt-get update
    sudo apt-get install -y python3

    if [[ -n "$parameter_file" ]]; then
        dx download "$parameter_file" -o parameter_file
        source parameter_file
    fi

    input_files=$(echo $view_files | jq .["$dnanexus_link"] -r | tr '\n' ',');
    input_files=${input_files::-1}
    echo "Input files: $input_files"

    if [[ "$genome" == "GRCh38" ]]; then
        genome=hg38;
    fi

    VIEWER_PROJECT="project-F5444K89PZxXjBqVJ3Pp79B4"
    viewer=`dx ls $VIEWER_PROJECT:/viewers/ProteinPaint\ \($genome\)\ \(VCF,\ bigWig\) --brief`

    shortcut=$(./make_dnanexus_shortcut -v $viewer \
                                        -f $input_files \
                                        -o $out_prefix.coverage \
                                        -p $VIEWER_PROJECT \
                                        --verbose)

    dx describe $out_prefix.coverage
    dx-jobutil-add-output viewer_shortcut --class=record $shortcut
}
