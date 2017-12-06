function run_rescue {
input_bam=$1
output_prefix=$2
bed_file=$3

samtools view -hb -q 1 $input_bam  > ${input_bam%.bam}_rescue.1.bam
samtools view -H $input_bam > ${input_bam%.bam}_rescue.2.sam
intersectBed -abam $input_bam -b $bed_file -wa | samtools view - | awk -F "\t" '{if($5==0){print}}' >> ${input_bam%.bam}_rescue.2.sam
samtools view -Sb ${input_bam%.bam}_rescue.2.sam > ${input_bam%.bam}_rescue.2.bam
nRescued=`samtools view ${input_bam%.bam}_rescue.2.bam | wc -l `
if [ "$nRescued" -gt 0 ];then
samtools merge ${input_bam%.bam}_rescue.bam ${input_bam%.bam}_rescue.1.bam ${input_bam%.bam}_rescue.2.bam
samtools sort ${input_bam%.bam}_rescue.bam $output_prefix
samtools index $output_prefix.bam
fi
}
