#!/usr/bin/env python3
## ---------------------------
##
## Script name: simFusion.py
##
## Purpose of script: Use RSVSim + InsilicoSeq to insert simulated translocations into hg19 reference genome and generate Illumina NovaSeq Paired-End Reads
##
## Author: Ankit Jambusaria
##
## Date Created: 10-07-21
##
## Email: ajambusaria@ambrygen.com
##
## ---------------------------
##
## conda_env: simFusion (/mnt/ajambusaria/anaconda3/envs/simFusion)
##
## ---------------------------
##
## Inputs:
## 1. geneAchr 
## 2. geneAstart 
## 3. geneAend
## 4. geneBchr  
## 5. geneBstart  
## 6. geneBend  
## 7. fusionName 
## 8. threads 
## 9. bedFile 
## 10. nReads 
## 11. outputPath
##
## Outputs:
## 1. reads_R1.fastq.gz
## 2. reads_R1.fastq.gz
## 3. reads_abundance.csv
## 4. genome_rearranged.fasta
## 5. genome_rearraned.fasta.fai
## 6. subset_genome_rearranged.fasta
## 7. translocations.csv
##
## ---------------------------
## Test 

'''
python3 simFusion.py \
--geneAchr chr2 \
--geneAstart 42396477 \
--geneAend 42559688 \
--geneBchr chr2 \
--geneBstart 29416090 \
--geneBend 30144477 \
--fusionName EML4-ALK \
--threads 32 \
--bedFile /mnt/ajambusaria/projects/liquid_biopsy/simulateFusions/simulations/EML4-ALK-2/ref/eml4_alk.bed \
--nReads 100k \
--outputPath /mnt/ajambusaria/projects/liquid_biopsy/simFusion/output/EML4-ALK/
'''

## ---------------------------

import argparse
import os
import logging
import time

# Paths 
simTranslocationPath = "/mnt/ajambusaria/projects/liquid_biopsy/github/useful_bioinfo_tools/other_tools/simFusion////"
issPath = "/mnt/ssamadian/anaconda3/envs/stats/bin/"
bedtoolsPath = "/mnt/ajambusaria/anaconda3/envs/simFusion/bin/"

# Functions

def parse_args():
    '''
    Parse user inputs
    '''
    
    parser = argparse.ArgumentParser(description='Generate fusion positive sample paired-end reads.')
    # geneA
    parser.add_argument('--geneAchr', type=str,
                        help='Chromosome of geneA (e.g. chr2)', required=True)
    parser.add_argument('--geneAstart', type=str,
                        help='gene A hg19 start coordinate', required=True)
    parser.add_argument('--geneAend', type=int,
                        help='gene A hg19 end coordinate', required=True)
    #geneB
    parser.add_argument('--geneBchr', type=str,
                        help='Chromosome of geneB (e.g. chr2)', required=True)
    parser.add_argument('--geneBstart', type=str,
                        help='gene B hg19 start coordinate', required=True)
    parser.add_argument('--geneBend', type=int,
                        help='gene B hg19 end coordinate', required=True)    
    # script params
    parser.add_argument('--fusionName', type=str,
                        help='Fusion Name (e.g EML4-ALK)', required=True)
    parser.add_argument('--nReads', type=str,
                        help='Number of Reads (i.e. 100k, 1m, 10m, 25m, 50m', required=True)
    parser.add_argument('--bedFile', type=str,
                        help='Path to bedfile to subset rearranged genome', required=True)
    parser.add_argument('--threads', type=int,
                        help='Number of threads', required=True)
    parser.add_argument('--outputPath', type=str,
                        help='Directory of output', required=True)
    return parser.parse_args()

def job_submitter(threads, jobpath, jobname, cmd, outdir):
    with open(jobpath+'.job', 'w') as out:
        out.write("#! /bin/bash\n"\
        "#$ -j y\n"\
        "#$ -cwd\n"\
        "#$ -S /bin/bash\n"\
        "#$ -q test.q\n"\
        "#$ -M ajambusaria@ambrygen.com\n")
        out.write("#$ -pe mpi {}\n\n".format(threads))
        out.write("#$ -N {}\n\n".format(jobname))
        out.write('source ~/.bashrc\n')
        out.write('cd {}\n'.format(outdir))
        out.write(cmd)
    qsubCmd = 'qsub {jobpath}.job'.format(jobpath=jobpath)
    print(qsubCmd)
    os.system(qsubCmd)

def run_RSVSim(geneAchr, geneAstart, geneAend, geneBchr, geneBstart, geneBend, fusionName, outputPath, simTranslocationPath, threads):
    '''
    run simTranslocations.R which takes in user inputs and generates fusion positive reference fasta
    Inputs: 
        geneA: chr, start, end
        geneB: chr, start, end 
        fusionName
    Outputs:
        genome_rearraged.fasta
        translocations.csv
    '''
    print(os.path.isfile(outputPath))
    if os.path.isfile(outputPath)==True:
        os.mkdir(outputPath)
    
    rsvsimCmd = "Rscript {simTranslocationPath}simTranslocation.R \
    -a {geneAchr} \
    -b {geneAstart} \
    -c {geneAend} \
    -d {geneBchr} \
    -e {geneBstart} \
    -f {geneBend} \
    -n {fusionName} \
    -o {outputPath}".format(
        geneAchr=geneAchr, geneAstart=geneAstart, geneAend=geneAend, 
        geneBchr=geneBchr, geneBstart=geneBstart, geneBend=geneBend, 
        fusionName=fusionName, outputPath=outputPath, simTranslocationPath=simTranslocationPath
    )
    print(rsvsimCmd)
    # submit job
    print("Submitting RSVSim job")
    jobPath = outputPath
    jobName = fusionName
    job_submitter(threads, jobPath, jobName, rsvsimCmd, outputPath)
    print("RSVSim job submitted!")

    
def subsetGenome(bedtoolsPath, outputPath, bedFile):
    print("Waiting for genome_rearranged.fasta to be produced...")
    # wait for genome_rearranged.fasta to exist
    time_to_wait = 100
    time_counter = 0
    while not os.path.exists(os.path.join(outputPath, "genome_rearranged.fasta")):
        time.sleep(60)
        if time_counter > time_to_wait:break

    print("genome_rearranged.fasta exists!")
    print("Subsetting Genome")
    
    subsetCmd = "{bedtoolsPath}bedtools getfasta -fi {outputPath}genome_rearranged.fasta \
    -bed {bedFile} \
    -fo {outputPath}subset_genome_rearranged.fasta".format(
        bedtoolsPath=bedtoolsPath, outputPath=outputPath, bedFile=bedFile)                    
    print(subsetCmd)
    os.system(subsetCmd)                    
        
    
def run_InSilicoSeq(issPath, nReads, threads, outputPath, fusionName):
    '''
    run InsilicoSeq to simulate Illumina NovaSeq paired-end reads for genome_rearranged.fasta
    Inputs: 
        genome_rearranged.fasta
        nReads
    Outputs:
        reads_R1.fastq
        reads_R2.fastq
    '''
    # wait for subset_genome_rearranged.fasta to exist
    time_to_wait = 100
    time_counter = 0
    while not os.path.exists(os.path.join(outputPath, "subset_genome_rearranged.fasta")):
        time.sleep(60)
        if time_counter > time_to_wait:break
        
    print("subset_genome_arranged.fasta exists!")
        
    issCmd = "nohup {issPath}iss generate \
    --cpus {threads} \
    --genomes {outputPath}subset_genome_rearranged.fasta \
    --model novaseq \
    --abundance uniform \
    --n_reads {nReads} \
    --output {outputPath}{fusionName}_{nReads}".format(
        issPath = issPath, nReads=nReads, threads=threads, outputPath=outputPath, fusionName=fusionName
    )
    
    print(issCmd)

    # submit job
    print("Submitting InSilicoSeq job...")
    jobPath = outputPath
    jobName = fusionName
    job_submitter(threads, jobPath, jobName, issCmd, outputPath)
    
    
def gzipFastqs(outputPath, fusionName, nReads):
    '''
    gzip forward and reverse read fastqs
    '''
    
    print("Waiting for ", outputPath + fusionName + '_' + nReads + "_R2.fastq")
    # wait for reads to exist (e.g.EML4-ALK_100k_R2.fastq)
    time_to_wait = 100
    time_counter = 0
    while not os.path.exists(outputPath + fusionName + '_' + nReads + "_R2.fastq"):
        time.sleep(180)
        if time_counter > time_to_wait:break    
    
    print("Gzipping ....")
    gzipCmd = "gzip {outputPath}{fusionName}_{nReads}_R*".format(
        outputPath=outputPath, fusionName=fusionName, nReads=nReads)
    print(gzipCmd)
    os.system(gzipCmd)    
   

def clean_up(outputPath):
    cleanCmd = "rm -rf {outputPath}*iss.tmp*".format(outputPath=outputPath)
    print(cleanCmd)
    os.system(cleanCmd)
    
    
def main():
    ''' 
    This code simulates translocations into hg19 genome and generates Illumina NovaSeq paired-end reads
    '''
    
    ###### inputs ######
    print('Parsing Args')
    args = parse_args()
    geneAchr = args.geneAchr
    geneAstart =  args.geneAstart
    geneAend = args.geneAend
    geneBchr =  args.geneBchr    
    geneBstart = args.geneBstart
    geneBend =  args.geneBend    
    fusionName =  args.fusionName
    nReads =  args.nReads
    bedFile = args.bedFile
    threads = args.threads
    outputPath = args.outputPath
    
    ###### RSVSim ######
    print("Generating fusion positive reference")
    run_RSVSim(geneAchr, geneAstart, geneAend, geneBchr, geneBstart, geneBend, fusionName, outputPath, simTranslocationPath, threads)
    
    ###### Subset genome_rearranged ######
    subsetGenome(bedtoolsPath, outputPath, bedFile)
    print("subset_genome_rearranged.fasta has been generated!")
    
    ###### InSilicoSeq ######
    run_InSilicoSeq(issPath, nReads, threads, outputPath, fusionName)
    
    ###### gzip fastqs ######
    gzipFastqs(outputPath, fusionName, nReads)
    print("Reads have been generated!")
    
    ###### Cleanup Temp Files ######
    print("Cleaning up temp files...")
    clean_up(outputPath)
    print("Clean up finished!")
    print("Simulation Pipeline Complete!")
    
    
if __name__=="__main__":
#     logging.basicConfig(level=logging.DEBUG, filename="logfile", filemode="a+",
#                         format="%(asctime)-15s %(levelname)-8s %(message)s")
#     logging.info("hello")
    main()
