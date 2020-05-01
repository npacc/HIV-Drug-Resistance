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

//SHIVER INIT
params.shiverinitdir = "/home/centos/nextflow/HIV-Drug-Resistance/config/hiv/shiver_init_HIV/"
params.shiverconf = "/home/centos/nextflow/HIV-Drug-Resistance/config/hiv/shiver_init_HIV/config/config.sh"
ShiverInitDir = Channel.fromPath( "${params.shiverinitdir}", type: 'dir')
Channel.fromPath( "${params.shiverconf}").set{ShiverConf1}
Channel.fromPath( "${params.shiverconf}").set{ShiverConf2}

process ShiverInit {
    container "/home/centos/nextflow/Def-Files/singularity-files/shiver-test.sif"
    
    input:
    file(dir) from ShiverInitDir
    file(conf) from ShiverConf1

    output:
    file("shiver_init_HIV/") into ShiverInit

    script:
    """
    shiver_init.sh MyInitDir /home/centos/nextflow/HIV-Drug-Resistance/config/hiv/shiver_init_HIV/config/config.sh \
    /home/centos/nextflow/HIV-Drug-Resistance/config/hiv/shiver_init_HIV/ExistingRefAlignment.fasta \
    /home/centos/nextflow/HIV-Drug-Resistance/config/hiv/shiver_init_HIV/adapters.fasta \
    /home/centos/nextflow/HIV-Drug-Resistance/config/hiv/shiver_init_HIV/primers.fasta
    """
}

//SHIVER
process HIVShiver {
    container "/home/centos/nextflow/Def-Files/singularity-files/shiver-test.sif"

    publishDir "${FastqDir}/Shiver", pattern: "${dataset_id}.shiver.fa", mode: 'copy'
    publishDir "${FastqDir}/Shiver", pattern: "${dataset_id}.shiverlog.txt", mode: 'copy'   

    input:
    tuple dataset_id, file(assembly), file(forward), file(reverse), file(shiverconf), file(shiverinit) from HIVIVAAssembly.join(HIVCleanReadsPolishing, by: [0]).combine(ShiverConf2).combine(ShiverInit)

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

/*
//HIV MAPPING VARIANT CALLING
process HIVMappingVariantCalling {
    cpus 4

    container "/home/centos/nextflow/Def-Files/singularity-files/minimap2-test.sif"

    input:
    tuple dataset_id, file(forward), file(reverse), file(assembly) from HIVCleanReadsVariantCalling.join(HIVAssemblyBAM, by: [0])

    output:
    tuple dataset_id, file("${dataset_id}.variants.bam") into HIVMappingNoDupsBAM

    script:
    """
    samtools faidx $assembly
    minimap2 -x sr -a $assembly $forward $reverse | samtools view -@ 2 -b | samtools sort -@ 2 -o ${dataset_id}.variants.bam
    """
}

//SETUP MINOR VARIANT FREQUENCY LIST
params.minvarfreq = ['0.2', '0.1', '0.01']
Channel.from(params.minvarfreq).set{ MinVarFreq}

//Set default variant strategy (VarScan, BCFtools or LoFreq)
params.variantstrategy = 'VarScan'

//HIV VARIANT CALLING VAR SCAN
process HIVVariantCallingVarScan {
    container "/home/centos/nextflow/Def-Files/singularity-files/variantcalling-test.sif"

    publishDir "${FastqDir}/CallVariants/fasta/minor_variants", pattern: "${dataset_id}.${minvarfreq}.minor.fa", mode: 'copy'
    publishDir "${FastqDir}/CallVariants/fasta/IUPAC", pattern: "${dataset_id}.${minvarfreq}.iupac.consensus.fa", mode: 'copy'
    publishDir "${FastqDir}/CallVariants/vcf", pattern: "${dataset_id}.${minvarfreq}.consensus.vcf", mode: 'copy'

    input:
    tuple dataset_id, file(bam), file(assembly), minvarfreq from HIVMappingNoDupsBAM.join(HIVAssemblyVariants, by: [0]).combine(MinVarFreq)

    output:
    tuple dataset_id, minvarfreq, file("*.${minvarfreq}.iupac.consensus.fa") into HIVAssemblyWithVariants
    file "${dataset_id}.${minvarfreq}.minor.fa"
    file "${dataset_id}.${minvarfreq}.consensus.vcf"

    script:
    """
    samtools mpileup --max-depth 10000000 --redo-BAQ --min-MQ 17 --min-BQ 20 --output ${dataset_id}.mpileup --fasta-ref ${assembly} ${bam}
    java -Xmx17G -jar /home/centos/miniconda3/envs/trialrun/share/varscan-2.4.4-0/VarScan.jar mpileup2cns ${dataset_id}.mpileup --min-var-freq ${minvarfreq} --p-value 95e-02 --min-coverage 100 --output-vcf 1 > ${dataset_id}.varscan.cns.vcf
    bgzip ${dataset_id}.varscan.cns.vcf
    tabix -p vcf ${dataset_id}.varscan.cns.vcf.gz
    bcftools view -i'FILTER="PASS"' -Oz -o ${dataset_id}.varscan.cns.filtered.vcf.gz ${dataset_id}.varscan.cns.vcf.gz
    zcat ${dataset_id}.varscan.cns.filtered.vcf.gz > ${dataset_id}.${minvarfreq}.consensus.vcf
    tabix -p vcf ${dataset_id}.varscan.cns.filtered.vcf.gz
    bcftools consensus -f $assembly ${dataset_id}.varscan.cns.filtered.vcf.gz --output ${dataset_id}.${minvarfreq}.minor.fa
    bcftools consensus -I -f $assembly ${dataset_id}.varscan.cns.filtered.vcf.gz --output ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i 's/polished/consensus-minor/g' ${dataset_id}.${minvarfreq}.minor.fa
    sed -i 's/polished/consensus-iupac/g' ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i '/^>/ s/\$/ [Variant caller: ${params.variantstrategy}] [Minor variant bases IUPAC] [Variant frequency: ${minvarfreq}] /' ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i '/^>/ s/\$/ [Variant caller: ${params.variantstrategy}] [Minor variants bases ONLY] [Variant frequency: ${minvarfreq}] /' ${dataset_id}.${minvarfreq}.minor.fa
    """
}

//SET MINOR VARIANT FREQUENCY AT WHICH TO CALL RESISTANCE (must be one of the values in params.minvarfreq)
params.selectedminvarfreq = '0.2'
HIVAssemblyWithSelectedMinVarFreq = HIVAssemblyWithVariants.filter{ it[1] == params.selectedminvarfreq }

//MAKE HIV RESISTANCE REPORT
process HIVMakeResistanceReport {
    container "/home/centos/nextflow/Def-Files/singularity-files/sierappy-test.sif"

    publishDir "${FastqDir}/Resistance/analysis/04-call_resistance", pattern: "${dataset_id}.json", mode: 'copy'
    publishDir "${FastqDir}/Resistance/analysis/05-generate_report", pattern: "${dataset_id}.rtf", mode: 'copy'
    publishDir "${FastqDir}/Resistance/reports", pattern: "${dataset_id}.rtf", mode: 'copy'

    input:
    tuple dataset_id, minvarfreq, file(variantassembly) from HIVAssemblyWithSelectedMinVarFreq

    output:
    file("${dataset_id}.json")
    file("${dataset_id}.rtf")

    script:
    episodenumber = dataset_id.split('-').last()
    """
    buildreport.pl -i ${variantassembly} -n ${dataset_id} -l PHW_Cardiff
    """
} 
*/