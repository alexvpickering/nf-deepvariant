//  example:
// nextflow run main.nf --reads "/run/media/alex/Seagate Backup Plus Drive/cincinatti_exome_data/FastQ/09KTF_TTAGGC_L002_R{1,2}_001.fastq.gz" --samplename 09KTF_TTAGGC
 
 // Define the default parameters
params.reads = ""
params.outdir = "./results"
params.samplename = ""
params.readgroup = "${params.samplename}"
params.reference = "/root/hg38/Homo_sapiens_assembly38.fasta"

// deepvariant parameters
params.model_type = "WES" // one of the following: [WGS,WES,PACBIO]
params.num_shards = "6"

/*
 * Create the `read_pairs_ch` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch } 


process map_and_sort {
	publishDir "${params.outdir}/${params.samplename}"

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

process unzip_vcf {
    publishDir "${params.outdir}/${params.samplename}", mode: 'copy'

    input:
    file vcfgz from vcf_ch

    output:
    file "${params.samplename}.vcf" into vcf_unzip_ch

    script:
    """
    gunzip -f $vcfgz
    """
}


process annotate_csq {
	publishDir "${params.outdir}/${params.samplename}", mode: 'copy'

	input:
	file vcf_unzip from vcf_unzip_ch
	

	output:
	file "${params.samplename}_csq.bcf" into bcf_csq_ch


	script:
	"""
	bcftools csq -f ${params.reference} \
	 -g ~/hg38/Homo_sapiens.GRCh38.100.gff3.gz \
	  $vcf_unzip -Ob -o ${params.samplename}_csq.bcf
	"""
}

