// This is the nextflow script to run the simulator

// TODOs:
// MAYBEs:
// - add a mode to create a paired normal sample BAM file

params.input = "examples/main_input.txt"
params.out_dir = "test_output/"
params.roi_bed = "/home/vngo/Work/liquid_biopsy_pipeline/mnv_calling_module/simulator/examples/roi_2.bed"
params.ref_fa = "/home/vngo/genomes/hg19_ambry.fa"
params.ART_PATH = "/home/vngo/tools/art_bin_MountRainier/art_illumina"
params.read1_prof = "art_models/ctDNAFusion_profileR1.txt"
params.read2_prof = "art_models/ctDNAFusion_profileR2.txt"
params.noise = 0.0 // this is for the input_splitter


params.BWA_INDEX_BASE = "/home/vngo/genomes/hg19_ambry.fa" // location of the ref index base for BWA

SENTIEON='/mnt/clinical/bin/sentieon/sentieon-genomics-202010/bin/sentieon'
GENOME='/mnt/swu/binary/LynchPanelPipeline_142genes/ref/whole_genome.fa'



process input_split{
    publishDir  "$params.out_dir/hap_gen_inputs/"

    input:
        path main_input from params.input
    
    output:
        file "HAP*.txt" optional true into hap_gen_input_ch
    
    script:
        """
        input_splitter.py -i $main_input -n ${params.noise} -o "./"
        """

}

process haplotype_generate{
    // generate seqs for each hap
    // convert the coordinates to the new format (rename chromosomes)
    echo true
    publishDir "$params.out_dir/haplotypes/"
    clusterOptions = "-S /bin/bash -V"
    conda "/home/vngo/anaconda3/envs/mnv_simulator"

    input:
        path hap_input from hap_gen_input_ch.flatten() // process each file separately
        path roi_bed from params.roi_bed
        path ref_fa from params.ref_fa
        
    output:
        set tag, file(hap_fasta), file(converted_hap_input) into hap_ch
    
    script:
        tag = hap_input.toString().replace("_input.txt", "")
        hap_fasta = tag+'.sim.fa'
        converted_hap_input = tag + '_conv_input.txt'
        //println(tag)
        //println(out_file)
        """
        haplotype_generator.py -i $hap_input -r $ref_fa -b $roi_bed -t $tag -o $hap_fasta -c $converted_hap_input
        """
}



process fastq_generate{
    // this is the main script to generate reads
    echo true
    publishDir "$params.out_dir/fastq/"
    clusterOptions = "-S /bin/bash -V"
    conda "/home/vngo/anaconda3/envs/mnv_simulator"
    
    input:
        set tag, file(hap_fasta), file(converted_hap_input) from hap_ch
        path read1_prof from params.read1_prof
        path read2_prof from params.read2_prof
    output:
        file("*_R?.fq") into hap_fastq_ch
    
    script:
        """
        generate_fastq_reads.py -a ${params.ART_PATH} -1 $read1_prof -2 $read2_prof -i $converted_hap_input -r $hap_fasta --out-prefix $tag
        """
}


process bwa_mem{
    // align the reads to the reference
    echo true
    publishDir "$params.out_dir/bam/"
    clusterOptions = "-S /bin/bash -V"
    conda "/home/vngo/anaconda3/envs/mnv_simulator"
    penv = 'smp'
    cpus = 4
    memory = "16G"
    
        
    input:
        path fastq_file from hap_fastq_ch.collect() // load all fastq files above to this process
        path ref_fa from params.BWA_INDEX_BASE
        path ref_fa_pac from params.BWA_INDEX_BASE + ".pac"
        path ref_fa_ann from params.BWA_INDEX_BASE + ".ann"
        path ref_fa_amb from params.BWA_INDEX_BASE + ".amb"
        path ref_fa_0123 from params.BWA_INDEX_BASE + ".0123"
        path ref_fa_bwt from params.BWA_INDEX_BASE + ".bwt.2bit.64"
        
    output:
        set file("sim.bam"), file("sim.bam.bai") into bam_ch
        
    script:
        """
        # combine fastq files into only two 
        cat *_R1.fq > combine_R1.fq
        cat *_R2.fq > combine_R2.fq
        
        # mapping
        bwa-mem2 mem -R '@RG\\tID:foo\\tSM:bar' -t 4 $ref_fa combine_R1.fq combine_R2.fq |samtools view -S -h -b |samtools sort -@ 4 -m 3G - > sim.bam
        
        # make index
        samtools index sim.bam
        """
}

process sentieon{
    echo true
    publishDir "$params.out_dir/sentieon_reports/", pattern:"*.{csv,txt,pdf}"
    publishDir "$params.out_dir/bam/", pattern: "*.bam"
    publishDir "$params.out_dir/vcf/", pattern: "*.vcf"
    
    clusterOptions = "-S /bin/bash -V"
    conda "/home/vngo/anaconda3/envs/mnv_simulator"
    penv = 'smp'
    cpus = 4
    memory = "16G"
    
    input:
        set file(sim_bam), file(sim_bai) from bam_ch
    
    output:
        file "*.txt" 
        file "*.bam"
        file "*.pdf"
        file "*.vcf"
        file "*.csv"
        
    script:
        """
        $SENTIEON driver -t 4 -r $GENOME -i $sim_bam \
          --algo GCBias --summary GC_SUMMARY.txt GC_METRIC.txt \
          --algo MeanQualityByCycle MQ_METRIC.txt \
          --algo QualDistribution QD_METRIC.txt \
          --algo InsertSizeMetricAlgo IS_METRIC.txt  \
          --algo AlignmentStat ALN_METRIC.txt
          
        # Mark Duplicates
        $SENTIEON driver -t 4 -i $sim_bam \
            --algo LocusCollector --fun score_info SCORE.gz
            
        $SENTIEON driver -t 4 -i $sim_bam \
            --algo Dedup --score_info SCORE.gz  \
            --metrics DEDUP_METRIC.txt sim.deduped.bam
        
        # In-del realignment 
        $SENTIEON driver -t 4 -r $GENOME \
            -i sim.deduped.bam --algo Realigner sim.realigned.bam
        
        # make recal_data.table
        $SENTIEON driver -t 4 -r $GENOME \
            -i sim.realigned.bam --algo QualCal RECAL_DATA.TABLE
        
        # Base Quality score recalibration
        $SENTIEON driver -t 4 -r $GENOME -i sim.realigned.bam \
            -q RECAL_DATA.TABLE --algo QualCal \
            RECAL_DATA.TABLE.POST
            
        $SENTIEON driver -t 4 --algo QualCal --plot \
            --before RECAL_DATA.TABLE --after RECAL_DATA.TABLE.POST RECAL_RESULT.csv
            
        $SENTIEON plot QualCal -o BQSR.pdf RECAL_RESULT.csv
        
        # calling variants
        $SENTIEON driver -t 4 -r $GENOME -i sim.realigned.bam \
            -q RECAL_DATA.TABLE --algo Haplotyper sim.vcf
        
        """   
}

/*
process filter{
    // this process is to filter for only reads in region of interest
    // this might or might not be useful, leave it unused for now


}
*/

