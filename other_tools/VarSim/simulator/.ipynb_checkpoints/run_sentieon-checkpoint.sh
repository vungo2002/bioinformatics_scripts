#!/bin/bash


/mnt/clinical/bin/sentieon/sentieon-genomics-202010/bin/sentieon driver -t 4 -r /mnt/swu/binary/LynchPanelPipeline_142genes/ref/whole_genome.fa -i sim.bam           --algo GCBias --summary GC_SUMMARY.txt GC_METRIC.txt           --algo MeanQualityByCycle MQ_METRIC.txt           --algo QualDistribution QD_METRIC.txt           --algo InsertSizeMetricAlgo IS_METRIC.txt            --algo AlignmentStat ALN_METRIC.txt

# Mark Duplicates
/mnt/clinical/bin/sentieon/sentieon-genomics-202010/bin/sentieon driver -t 4 -i sim.bam             --algo LocusCollector --fun score_info SCORE.gz

/mnt/clinical/bin/sentieon/sentieon-genomics-202010/bin/sentieon driver -t 4 -i sim.bam             --algo Dedup --score_info SCORE.gz              --metrics DEDUP_METRIC_TXT sim.deduped.bam

# In-del realignment 
/mnt/clinical/bin/sentieon/sentieon-genomics-202010/bin/sentieon driver -t 4 -r /mnt/swu/binary/LynchPanelPipeline_142genes/ref/whole_genome.fa             -i sim.deduped.bam --algo Realigner sim.realigned.bam

# make recal_data.table
/mnt/clinical/bin/sentieon/sentieon-genomics-202010/bin/sentieon driver -t 4 -r /mnt/swu/binary/LynchPanelPipeline_142genes/ref/whole_genome.fa             -i sim.realigned.bam --algo QualCal RECAL_DATA.TABLE

# Base Quality score recalibration
/mnt/clinical/bin/sentieon/sentieon-genomics-202010/bin/sentieon driver -t 4 -r /mnt/swu/binary/LynchPanelPipeline_142genes/ref/whole_genome.fa -i sim.realigned.bam             -q RECAL_DATA.TABLE --algo QualCal             RECAL_DATA.TABLE.POST

/mnt/clinical/bin/sentieon/sentieon-genomics-202010/bin/sentieon driver -t 4 --algo QualCal --plot             --before RECAL_DATA.TABLE --after RECAL_DATA.TABLE.POST RECAL_RESULT.csv

/mnt/clinical/bin/sentieon/sentieon-genomics-202010/bin/sentieon plot QualCal -o BQSR.pdf RECAL_RESULT.csv

# calling variants
/mnt/clinical/bin/sentieon/sentieon-genomics-202010/bin/sentieon driver -t 4 -r /mnt/swu/binary/LynchPanelPipeline_142genes/ref/whole_genome.fa -i sim.realigned.bam             -q RECAL_DATA.TABLE --algo Haplotyper sim.vcf