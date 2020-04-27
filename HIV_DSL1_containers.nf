#!/usr/bin/env nextflow

//Second DSL1 HIV Pipeline written using containers outside of the earlier conda environment

//PUSHING TO GIT
//git add HIV_DSL1_containers.nf
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
    container "/home/centos/nextflow/Def-Files/singularity-files/trimgalore-test.sif"

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

//MULTI QC REPORT
process MultiQC {
    container "/home/centos/nextflow/Def-Files/singularity-files/multiqc-test.sif"

    publishDir "${FastqDir}/TrimmedGalore/MultiQC"

    input:
    file("*") from txtfiles.collect().combine(htmlfiles.collect())

    output:
    file "multiqc_report.html"

    script:
     """
      multiqc -m cutadapt -m fastqc -n multiqc_report.html .
     """
}

/*
//SETUP REFERENCES 
Ch_HIVComp = Channel.fromPath( "${FastqDir}/HIVRefs/hivcomp2015.fasta" ) 
Ch_HIVHXB2 = Channel.fromPath( "${FastqDir}/HIVRefs/HXB2.fasta" )

//CLEAN HIV READS
process CleanHIVReads {
    container "/home/centos/nextflow/Def-Files/singularity-files/minimap2-test.sif"

    cpus 4

    input:
    tuple dataset_id, file(forward), file(reverse), file(ref) from TrimmedReads.combine(Ch_HIVComp)

    output:
    tuple dataset_id, file("${dataset_id}.clean_1.fq.gz"), file("${dataset_id}.clean_2.fq.gz") into HIVCleanReadsAssembly, HIVCleanReadsPolishing, HIVCleanReadsVariantCalling

    script:
    """
    minimap2 -x sr -a $ref $forward $reverse | samtools view -F 4 -b > clean.bam
    picard SamToFastq VALIDATION_STRINGENCY=LENIENT I=clean.bam F=${dataset_id}.clean_1.fq.gz F2=${dataset_id}.clean_2.fq.gz
    """
}
*/