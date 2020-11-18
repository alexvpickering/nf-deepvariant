//  example:
// nextflow run main.nf --reads "/run/media/alex/Seagate Backup Plus Drive/cincinatti_exome_data/vikram/FastQ/09KTH_CGATGT_L001_R{1,2}_001.fastq.gz" --samplename 09KTH_CGATGT
 
 // Define the default parameters
params.reads = ""
params.outdir = "./results"
params.samplename = ""
params.readgroup = "${params.samplename}"
params.reference = "/root/hg38/Homo_sapiens_assembly38.fasta"

// deepvariant parameters
params.model_type = "WES" // one of the following: [WGS,WES,PACBIO]
params.num_shards = "6"

// snpEff parameters
params.snpgen = "GRCh38.99"

/*
 * Create the `read_pairs_ch` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch } 


process map_and_sort {
	publishDir "${params.outdir}/${params.samplename}", mode: 'copy'

	input:
	tuple val(pair_id), path(reads) from read_pairs_ch

	output:
	file "${params.samplename}.sorted.bam" into bamfile_ch
	file "${params.samplename}.sorted.bai" into bamindex_ch
	
	"""
	bwa mem -t 6 -M -R '@RG\\tID:${params.readgroup}\\tSM:${params.samplename}\\tPL:ILLUMINA' $params.reference $reads \
	| gatk --java-options "-Djava.io.tmpdir=./tmp" SortSam -I /dev/stdin -O ${params.samplename}.sorted.bam --SORT_ORDER=coordinate --CREATE_INDEX=true
    """	
}

process run_deepvariant {
    publishDir "${params.outdir}/${params.samplename}", mode: 'copy'

    input:
    file bamfile from bamfile_ch
    file bamindex from bamindex_ch

    output:
    file "${params.samplename}.vcf.gz" into vcf_ch
    file "${params.samplename}.g.vcf.gz" into gvcf_ch

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=${params.model_type} \
        --ref=${params.reference} \
        --reads=$bamfile \
        --output_vcf=${params.samplename}.vcf.gz \
        --output_gvcf=${params.samplename}.g.vcf.gz \
        --num_shards=${params.num_shards}
    """
}

process annotate_vcf {
    publishDir "${params.outdir}/${params.samplename}", mode: 'copy'

    input:
    file vcfgz from vcf_ch

    output:
    file "${params.samplename}.ann.vcf" into vcf_ann_ch

    script:
    """
    java -Xmx8g -jar /gatk/snpEff/snpEff.jar ${params.snpgen} $vcfgz > ${params.samplename}.ann.vcf
    """
}


process filter_vcf {
    publishDir "${params.outdir}/${params.samplename}", mode: 'copy'

    input:
    file vcf_ann from vcf_ann_ch

    output:
    file "${params.samplename}.ann.pass.high_mod.vcf" into vcf_filt_ch

    script:
    """
    cat $vcf_ann | java -jar /gatk/snpEff/SnpSift.jar filter \
        "FILTER = 'PASS' & ((ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE'))"  \
        > ${params.samplename}.ann.pass.high_mod.vcf
    """
}
