
function run_bed2jsonbed {
input_bed=$1

awk '{print $1"\t"$2"\t"$3"\t{\"name\":\""$4"\"}"}' $input_bed |sort -k1,1 -k2,2n > ${input_bed%.bed}.bedj
bgzip ${input_bed%.bed}.bedj
tabix -p bed ${input_bed%.bed}.bedj.gz
}
