cpus=$1
read_group_tag=$2
ref_genome=$3
read1=$4
read2=$5
bedfile=$6
outbam=$7


bwa mem -R $read_group_tag -t $cpus $ref_genome $read1 $read2 |samtools view -b -h -L $bedfile - |samtools sort -@ $cpus - > $outbam
