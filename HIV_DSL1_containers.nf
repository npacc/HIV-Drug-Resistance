#!/usr/bin/env nextflow

//Second DSL1 HIV Pipeline written using containers outside of the earlier conda environment

//PUSHING TO GIT
//git add test_trim.nf
//git commit -m 'message'
//git push origin master

// INPUT RUN FOLDER
FastqDir = "/home/centos/nicole_HIV_fastq"

// INPUT READS
FastqFiles = "${FastqDir}/*_R{1,2}_001.fastq.gz"

// CHANNEL
PairFilesChannel = Channel.fromFilePairs( "${FastqFiles}", flat: true )

//TRIM READS
process TrimmedGalore {
    container "/home/centos/nextflow/singularity/singularity-files/test_Matt.sif"

    cpus 2
    
    publishDir "${FastqDir}/TrimmedGalore/TrimmedReads", pattern: '*_val_{1,2}.fq.gz'
    publishDir "${FastqDir}/TrimmedGalore/htmlfiles", pattern: '*_fastqc.{zip,html}'
    publishDir "${FastqDir}/TrimmedGalore/txtfiles" , pattern: '*_trimming_report.txt'

    input:
    tuple dataset_id, file(forward), file(reverse) from PairFilesChannel

    output:
    tuple dataset_id, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into TrimmedReads
    file("*trimming_report.txt") into txtfiles
    file("*_fastqc.{zip,html}") into htmlfiles

    script:
    """
    trim_galore --fastqc --paired ${forward} ${reverse}
    """
}