params.input_bam = "s3://honeybee-s3/data/run27-29/dragen_enrichment_germline/*HYB*/*.bam"
params.ds_factor = null 
params.outdir = null 
//params.scratch_dir = "/data/"

// -w s3://nfwork 
in_bam_files = Channel.fromPath(params.input_bam)


process downsample_bam{
    publishDir = params.outdir 
    container = 'broadinstitute/picard:latest'
    memory = "16GB"
    cpus = 4
    queue = 'nf-500gb'
    errorStrategy { sleep( 200 as long); return 'retry' } // This can deal with network issues that leadto failure, sleep for 3 mins everytime it fails
    maxRetries 30 // max 5 total, not per sample --> ideally 5*number_of_samples

    
    //scratch = params.scratch_dir
    
    input:
        path in_bam from in_bam_files
    output:
        file "*.ds.bam"
    script:
        """
        infile=$in_bam
        java -jar /usr/picard/picard.jar  DownsampleSam  I=\$infile  O=\${infile/.bam/.ds.bam} P=${params.ds_factor}
       
        """
}