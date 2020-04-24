#!/usr/bin/env nextflow

//Original DSL1 HIV piepline made in the following environment: conda activate trialrun

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
//Need to specify the 'fromFilePairs' as we have R1 and R2. 
// $ means vairable
// Flat: true removes the nested sqr brackets
//PairFilesChannel.println() //This can print the contents of the channel but you must first */ /* command out the process below

//TRIM READS
process TrimmedGalore {
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
//You can check the contents of the channels produced using the following:
    //TrimmedReads.println()
    //txtfiles.println()
    //htmlfiles.println()
    //txtfiles.collect().combine(htmlfiles.collect()).println() 
//publishDir pattern (specify a different glob string using the option pattern to store into each directory only the files that match the provided pattern)


//MULTI QC REPORT
process MultiQC {
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
//Multiqc aggregates results from multiple samples to form a single report
    //-m (use only this module)
    //-n (for the fastqc report filename)
    //. (the dot at the end says to process the files found in the current directory)
//To view the MultiQC Report, copy over the html file to local drive, do this on Ubuntu: 
//scp centos@nicole-dev:/home/centos/nicole_HIV_fastq/TrimmedGalore/MultiQC/*.html /mnt/c/Users/ni122520/Documents/hiv_trial */

//SETUP REFERENCES 
Ch_HIVComp = Channel.fromPath( "${FastqDir}/HIVRefs/hivcomp2015.fasta" ) 
Ch_HIVHXB2 = Channel.fromPath( "${FastqDir}/HIVRefs/HXB2.fasta" ) 

//In the PenGU pipeline, the clean reads are compared with the refs to get rid of any non-HIV reads, this is not necessary here as they are all HIV
//Ch_HIVComp.println() 
//Ch_HIVHXB2.println() 


//CLEAN HIV READS
process CleanHIVReads {
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
// Check the channel contents using: TrimmedReads.combine(Ch_HIVComp).println()
// Minimap2 is the tool which actually finds the overlaps between the read and the reference genome specified 
// Minimap2 then outputs this into SAM format and samtools allows us to view the file by converting it into a .BAM file
    // -x (preset string) sr (genomic short-read mapping) -a (output as SAM format)
    //Then the reference (called target in documentation), then any query sequences (.fa/.fastq files)
//Samtools converts SAM to BAM:
    //view (prints all the alignments in the specified input file to stdout)
    //-F 4 (only print reads which do not have a #4 flag - and were succesfully mapped, therefore removing contaminants) 
    //-@ 2 (Number of BAM compression threads to use in addition to main thread)
    //-b (output in the BAM format) > clean.bam
//Picard converts the Sam to Fastq:
    //Validation Stringency Lenient (Emit warnings but keep going if possible)
    //I (input SAM/BAM file to extract reads)
    //F (fastq file output) and F2 (second end fastq file)


//SAMPLE HIV READS
//This step is no longer required as the computing power has increased

//ASSEMBLE HIV READS
process AssembleHIVReads {
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
//IVA is a de novo assembler designed to assemble virus genomes that have no repeat sequences from mixed populations at extremely high and variable depth.
//When it comes off the sequencer, you collapse a single pair of reads into a consensue by comparing the two and creating one final merged contig. 
//This is limited to the size of the original physical fragment because those reads don't have any information outside of themselves.
    //-v prints the messages to the stdout
    //-f & -r are the forward and reverse reads
    //output is 'iva_assembly'
    //2>&1 > puts the stderror into the stdoutput and both of these go to the iva.log file
    //> always points to a file, and means, any stdout from the preceeding process, put into a file. Its called a re-direct 
    //| joins processes together, it takes the stdout from the process on its left, and shoves it into the stdin of the process on the right
    //Then, rename the contigs.fasta file to include the dataset_id 
    //If the first command fails then rename the log to iva.fail.log so we know something went wrong

//ORDER HIV CONTIGS & GAP FILL HIV CONTIGS
//These two steps are redundant as Shiver does them below

//HIV SHIVER

// Setup init_dir for shiver
ShiverInit = Channel.fromPath( "${FastqDir}/nextflow_pipelines/config/hiv/shiver_init_HIV/shiver_init_HIV" )
ShiverConf = Channel.fromPath( "${FastqDir}/nextflow_pipelines/config/hiv/shiver_init_HIV/config.sh" )

//In order to work with Shiver you need to initialise it (only once)
//Check that the version of shiver you have has the correct config file - it changes with every version
    //Once you know your version, find the correct config file on the shiver releases GitHub page. Copy the .tar.gz file and use wget to save it to the folder
        //tar.gz is the linux equivalent to a zip file so we need to unpack it using:
            //tar -xvzf v1.3.5.tar.gz (x expands, v verbose, z uses gzip, f the next thing is a file)
        //Then, check your old config file and move across any missing things to the new one making sure not to change variable names.
//You also need to have bc installed (Unix Calculator) for shiver to work: sudo yum install bc
//You also need to have the correct shiverinit contents. The files which can stay the same between versions are:  
    //primers.fasta, adapters.fasta, config.sh, hivcomp2015.fasta 
    //Once all these are in the folder do the following command: shiver_init.sh shiver_init_HIV config.sh hivcomp2015.fasta adapters.fasta primers.fasta
//If this does not work because of a "lengths" issue, then you need to create a new MSA using sequences from the hiv.lanl.gov database
    //Shiver wants to compare sites in the genome that correspond to each other, so after alignment, corresponding sites ends up in a signle column
    //To create an MSA from the hicvcomp2015.fasta file you use a tool call mafft 
    //mafft --auto hivcomp2015.fasta >  hivcomp2015.aln.fasta
//Make sure you then point the ShiverInit Channel to the directory with the new config.sh file and the new shiver_init_HIV dir

//HIVIVAAssembly.join(HIVCleanReadsPolishing, by: [0]).combine(ShiverConf).combine(ShiverInit).println()

process HIVShiver {
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
//Shiver produces a nicely polished fasta file: this is the "consensus" and it contains the most common base at each position
//Shiver decides if the contig needs correcting and if so creates a cut_wRefs.fasta file
//Shiver maps the reads to our custom reference, it scaffolds the sequence assembly against the HIV reference genoma
    //For shiver, you need: assembly, clean reads, config file, config directory
//publishDir (makes a copy of the files in the output into the location stated). They must be in the output to copy them in PublishDir.
//input: This is a long section with a number of operators:
    //Some of the files are taken from the ShiverConf and ShiverInit Channel defined at the start. 
    //The other files are taken from the earlier HIVIVAAssembly output channel & HIVCleanReadsPolishing Channel
    //HIVIVAAssembly is joined to HIVCleanReadsPolishing using the join operator which finds a common element and matches on this
        //In this case multiple matching elements defined by [0]: this is the first element, in our case: dataset_id
        //This is then combined with the ShiverConf & ShiverInit Channels
//output: optional true (this section allowd the process to run even if the desired output is not generated, the file may be legit missing)
//script: There are two else/if statements
    //The first if statement means run shiver_align_contigs.sh (which requires the following elements: config directory, config file, assembly file, dataset_id)
        //If this is run sucessfully, go to the next bit of the script, if not, then go to the else echo section and save the log file.
    //The next section: [ -f ${dataset_id}_cut_wRefs.fasta ] checked to see if a FILE with that name exists in the work dir
    //The `raw' alignment of contigs is the alignment of all contigs thought to be HIV, without any modification, to the existing references (this is ALWAYS produced)
    //The `cut' alignment is the alignment of the same references after they have been automatically corrected by \shiv, to the existing references (this MIGHT be produced)
    //shiver_align_contigs (this command will produce a file named ID.blast)
        //This details the blast hits of your contigs to those existing references supplied
        //Assuming at least one contig looks like HIV there will be a file called raw_wRefs.fasta 
    //When you have an alighment of contigs that you are happy with, you do map_reads.sh


//seqtk(tool for processing sequences in the FASTA/Q format)
//shiver outputs two sequences in the fasta file, one is the 3.5kb fragment and the second is the whole genome attached. We are only interested in the fragment.
    //There is a single MinCov_15_30.fasta per sample which contains these two sequences
//seqtk seq -l0 ${dataset_id}_remap_consensus_MinCov_15_30.fasta:
    //Print each fasta sequence on one line, so you have:
    //>IDENTIFIER SEQUENCE 1 (3.5kb fragment only)
    //AGTAGATCTAGAGAG
    //>IDENTIFIER SEQUENCE 2 (whole genome added to either end of fragment)
    //AGTAGATCTAGAGAG
//head -n2 
    //takes the first 2 lines only (3.5kb fragment identifier and sequence)
//sed />/!s/-//g:
    //Read this from right to left: replace with a blank whenever a '-' character is present except for when a '>' is present
//sed s/\\?/N/g
    //Replace with a 'N' every time a '?' is found
    //The two slashes before the ? indicate to ignore it's special metacharacter function
//sed s/_remap_consensus//g 
    //replace with a blank everytime you see '_remap_consensus'
//seqtk seq -l80 
    //Output these as 80 bases per line (rather than all in one line as above) and save as ${dataset_id}.shiver.fa


//HIV MAPPING VARIANT CALLING
process HIVMappingVariantCalling {
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
//Variant calling is the process of mapping sequence reads to a reference emitting a bam file. 
    //This generates pileups and looks for evidence of structural variation. 
    //This section maps back to the scaffold assembly using minimap2
    //Variant calling identifys small mismatches in the alighments which may represent mutations present in the sample.
//Input: this HIVAssemblyBAM consists of the MinCov_X_Y.fasta file which is the pairwise alighment of the consensus genome and the reference used for mapping.
    //This is joined according to dataset_id with the two original gz fastq files (read 1 & read 2)
//Output: a variants.bam file
    //Check channel contents: HIVMappingNoDupsBAM.println()
//Script: faidx is a tool which indexes or queries regions from a fasta file. 
    //Faidx-indexed reference file will enable base alighment quality calculations for all reads aligned to a reference in the file (process below).
//minimap2: -x sr (this is the preset which indicates it is short single-end reads without splicing), -a (generate CIGAR output alignments in SAM format)
//samtools view -@ 2 -b (Number of BAM compression threads to use in addition to main thread and output in the BAM format) 
// samtools sort -@ 2 -o ${dataset_id}.variants.bam (-o says to output alighnments to FILE which is dataset_id.variants.bam)


//SETUP MINOR VARIANT FREQUENCY LIST (proportions at which to call variants e.g. 0.2 == 20%)
//Variant allele frequency (VAF) is the percentgae of sequence reads observed matching a specific DNA variant divided by the overall coverage at that locus
//This is a measure of the proportion of DNA molecules in the original specimen carrying the variant (ours is set as 20% below)
params.minvarfreq = ['0.2', '0.1', '0.01']
Channel.from(params.minvarfreq).set{ MinVarFreq}

//Set default variant strategy (VarScan, BCFtools or LoFreq)
//Once you have the pileup you use VarScan to do some inferences on it
params.variantstrategy = 'VarScan'

//HIV VARIANT CALLING VAR SCAN
process HIVVariantCallingVarScan {
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
//This section take the pileup (samtools) and identifies variants (Varscan)
//Using the varscan output, you then generate a consensus sequence (BCFtools)
//Script: samtools mpileup: this produces a "pileup" textual format from an alignment. Each input file produces a sperate group of pileup columns in the output.
//A pileup summarizes each base call of aligned reads to a ref sequence.
    //max-depth (at a position, read maximally 10,000,000)
    //redo-BAQ (Base Alignment Quality is the Phred-scaled probability of a read base being misaligned)
        //It greatly helps to reduce false SNPs caused by misalignments. This command recalculates BAQ on the fly
    //min-MQ 17: Minimum mapping quality for an alignment to be used
    //min-BQ 20: Minimum base quality for a base to be considered
//The next section takes the vcf from varscan and augments the assembly based on variants that reach the minVarFreq value
//So, for example, every variant that gets a minVarFreq of >20 gets a PASS in the FILTER column
    //-J-Xmx17g says the SQ will use 17GB of RAM.
    //A Jar is a package file format used to aggregate many Java class files and associated metadata and resources (text, images, etc.) into one file for distribution
        //To find the .jar file path use: find /home/centos/miniconda3/envs/trialrun -name "*can.jar"
    //mpileup2cns ${dataset_id}.mpileup: This is multi-sample calling
    //--min-var-freq ${minvarfreq}: Minimum variant allele frequency threshold - these were set earlier as a variable minvarfreq
    //p-value 95e-02 (default p value threshold for calling variants) --min-coverage 100 (Min read depth at a position to make a call) 
//bgzip: Bgzip compresses files in a similar manner to, and compatible with, gzip(1)
        //The file is compressed into a series of small (less than 64K) 'BGZF' blocks. 
        //This allows indexes to be built against the compressed file and used to retrieve portions of the data without having to decompress the entire file.
//tabix -p vcf ${dataset_id}.varscan.cns.vcf.gz: Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file
        //-p: Input format for indexing is vcf
//bcftools view -i'FILTER="PASS"' -Oz -o ${dataset_id}.varscan.cns.filtered.vcf.gz ${dataset_id}.varscan.cns.vcf.gz
//This section makes the 3 filtered vcfs (these now include variants which are present at the various thresholds rather than the consensus)
        //Only include sites for which the FILTER=PASS is true (i.e. the variants are >20,10 or 1% of the reads)
        //oz: Output type is compressed VCF and output file
    //zcat ${dataset_id}.varscan.cns.filtered.vcf.gz > ${dataset_id}.${minvarfreq}.consensus.vcf
        //This uncompresses the filtered.vcf.gz and the stdout is re-directed to another file .consensus.vcf
    //tabix -p vcf ${dataset_id}.varscan.cns.filtered.vcf.gz
        //This dictates the input format for indexing as vcf
//There are 2 bcftools consensus parts: minor.fa and iupac.consensus.fa is produced 
    //minor.fa contains only the variant itself replacing the old consensus base
    //iupac.consensus.fa represents both the consensus and the minor variant using an IUPAC code. 
        //The resistance prediction tool will tell you the results for both the variant and consesnus it finds in that position so we don't lose info about the major variant in that position
//bcftools consensus -f $assembly ${dataset_id}.varscan.cns.filtered.vcf.gz --output ${dataset_id}.${minvarfreq}.minor.fa
        //Using the consensus command you can create a consensus sequence for an invididual where the sequence incorporates variants types for this individual.
        //So you take the reference sequence (the assembly from shiver) and a vcf file and it plugs the variants from the vcf file into the consesnsus at the right location
        //-f: applies the filters and output as minor.fa
//bcftools consensus -I -f $assembly ${dataset_id}.varscan.cns.filtered.vcf.gz --output ${dataset_id}.${minvarfreq}.iupac.consensus.fa
        //Same as above but only include sites for which the above expression (FILTER=PASS) is true. 
        //When more than one nucleotide exceeds the interpretation threshold, the app reports the IUPAC code as the consensus for the position.
//sed -i 's/polished/consensus-minor/g' ${dataset_id}.${minvarfreq}.minor.fa
        //replace polished with consensus-minor for the .minor.fa
//sed -i 's/polished/consensus-iupac/g' ${dataset_id}.${minvarfreq}.iupac.consensus.fa
        //replace polished with consensus-iupac for the iupac.consensus.fa file
//sed -i '/^>/ s/\$/ [Variant caller: ${params.variantstrategy}] [Minor variant bases IUPAC] [Variant frequency: ${minvarfreq}] /' ${dataset_id}.${minvarfreq}.iupac.consensus.fa
        //The ^ symbol anchors the match to the start of the line, so "any line that starts with a >"
        //Is substituted at the end with the squr brackets
        //The single \ escapes the special meaning of the dollar sign in nextflow but keeps its special meaning in sed, which is 'the end of the line'
//sed -i '/^>/ s/\$/ [Variant caller: ${params.variantstrategy}] [Minor variants bases ONLY] [Variant frequency: ${minvarfreq}] /' ${dataset_id}.${minvarfreq}.minor.fa
        //Same as above but ammended to represent minor variants only 


//SET MINOR VARIANT FREQUENCY AT WHICH TO CALL RESISTANCE (must be one of the values in params.minvarfreq)
//When setting the filter, make sure to remmeber it uses 0 indexing, so print the channel and look for element 0, element 1, element 2 etc and put that in sqr bracket
params.selectedminvarfreq = '0.2'
HIVAssemblyWithSelectedMinVarFreq = HIVAssemblyWithVariants.filter{ it[1] == params.selectedminvarfreq }

//MAKE HIV RESISTANCE REPORT
process HIVMakeResistanceReport {

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
//Submit the sequence to Standford Database (sierrapy) and generate report from the Stanford output
// Stanford HIV Database: https://hivdb.stanford.edu/pages/FAQ/FAQ.html

//This section uses the 20% as the filter option. It filters the earlier HIVAssemblyWithVariants Channel and only keeps the frequency at 20%
//Script: 
    //episodenumber = dataset_id.split('-').last()
        //Written in Groovy: take dataset_id, split it into a list of elements every time you encounter a '-', then take the last element of that list
    //buildreport.pl is a perl script found in the singularity recepies online 
    //It needs to be saved in the bin directory, nextflow will always automatically look here for an executable file 
    //The -t, -n and -l are all perl command line options
    //If it is not being saved you may need to make the dir executable 
        //chmod +x buildreportl.pl the +x means ADD EXECUTE permission 
    //In order to write to the local drive where you need to add the phw.jpg picture you need to use the sudo before the wget and also the raw url
//Buildreport.pl: 
    //JSON: The Json is what gets returned from sierrapy.The perl script gives sierappy the fasta assembly and sierrapy sends that to the HIVdb API
        //Then, HIVdbAPI searches the sequences against its database and returns a JSON file to the sierappy client containing all of the info about the sequence
        //The JSON stays with the Bioinformatics team as it has more info than the clinicans needs
    //RTF: This goes to the clinicians and is authorised by Sally 
        //Eventually, we want to create a pipeline which will notify the members of each workstream when their workflow has finished via email

