// This script converts a position-sorted BAM file into fastq.gz files
params.input_bam = null
params.output_dir = null 
params.mem = "16GB"
params.cpus = 1


in_bam_files = Channel.fromPath(params.input_bam)


process bam2fastq{
    publishDir params.output_dir 
    container = 'broadinstitute/picard:latest'
    cpus = params.cpus
    memory = params.mem 
    queue = 'nf-500gb'


    input:
        path inBam from in_bam_files
    
    output:
        file "*.fq.gz"

    script:
        """
        infile=$inBam 
        ls . > listfiles.txt

        java -jar /usr/picard/picard.jar SamToFastq INPUT=$inBam \\
                FASTQ=\${infile/.bam/.R1.fq} \\
                SECOND_END_FASTQ=\${infile/.bam/.R2.fq}    
  
        gzip *.fq 
        """
}