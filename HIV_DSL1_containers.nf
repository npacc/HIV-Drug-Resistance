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

/*
//ASSEMBLE HIV READS
process AssembleHIVReads {
    container "/home/centos/nextflow/Def-Files/singularity-files/iva-test.sif"

    cpus 4

    input:
    tuple dataset_id, file(forward), file(reverse) from HIVCleanReadsAssembly

    output:
    tuple dataset_id, file("${dataset_id}.iva.fa") optional true into HIVIVAAssembly
    file("*.log")
    script:
    """
    if iva -v -f $forward -r $reverse iva_assembly 2>&1 > ${dataset_id}.iva.log ; then
      mv iva_assembly/contigs.fasta ${dataset_id}.iva.fa
    else
      mv ${dataset_id}.iva.log ${dataset_id}.iva.fail.log
    fi
    """
}

//HIV SHIVER

// Setup init_dir for shiver
ShiverInit = Channel.fromPath( "${FastqDir}/nextflow_pipelines/config/hiv/shiver_init_HIV/shiver_init_HIV" )
ShiverConf = Channel.fromPath( "${FastqDir}/nextflow_pipelines/config/hiv/shiver_init_HIV/config.sh" )

process HIVShiver {
    container "/home/centos/nextflow/Def-Files/singularity-files/shiver-test.sif"

    publishDir "${FastqDir}/Shiver", pattern: "${dataset_id}.shiver.fa", mode: 'copy'
    publishDir "${FastqDir}/Shiver", pattern: "${dataset_id}.shiverlog.txt", mode: 'copy'   

    input:
    tuple dataset_id, file(assembly), file(forward), file(reverse), file(shiverconf), file(shiverinit) from HIVIVAAssembly.join(HIVCleanReadsPolishing, by: [0]).combine(ShiverConf).combine(ShiverInit)

    output:
    tuple dataset_id, file("${dataset_id}.shiver.fa") optional true into HIVAssemblyBAM, HIVAssemblyVariants
    file("${dataset_id}.shiver.txt") optional true
    file("${dataset_id}.shiverlog.txt") optional true

    script:
     """
    if shiver_align_contigs.sh ${shiverinit} ${shiverconf} ${assembly} ${dataset_id}; then
        if [ -f ${dataset_id}_cut_wRefs.fasta ]; then
            shiver_map_reads.sh ${shiverinit} ${shiverconf} ${assembly} ${dataset_id} ${dataset_id}.blast ${dataset_id}_cut_wRefs.fasta ${forward} ${reverse}
        else
            shiver_map_reads.sh ${shiverinit} ${shiverconf} ${assembly} ${dataset_id} ${dataset_id}.blast ${dataset_id}_raw_wRefs.fasta ${forward} ${reverse}
        fi
        seqtk seq -l0 ${dataset_id}_remap_consensus_MinCov_15_30.fasta | head -n2 | sed '/>/!s/-//g' | sed 's/\\?/N/g' | sed 's/_remap_consensus//g' | seqtk seq -l80 > ${dataset_id}.shiver.fa
    else
        echo "No HIV contigs found. This sample is likely to be purely contamination" > ${dataset_id}.shiverlog.txt
    fi
    """
}
*/