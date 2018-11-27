params.bed = false
params.fasta = false
params.bams = false
params.outdir = false

if( !params.fasta ) {
    exit 1, "No reference fasta was supplied"
}
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
    faidx = file(fasta + ".fai")
    log.info("Reference fasta: ${fasta}")
}
if ( params.bed ){
    bed_file = file(params.bed)
    if ( !bed_file.exists() ) exit 1, "Bed file specified [${params.bed}], but does not exist"
    log.info("Excluded regions: ${params.bed}")
}
outdir = file(params.outdir)

Channel
    .fromPath(params.bams, checkIfExists: true)
    .map { file -> tuple(file.baseName, file, file + '.bai') }
    .into { call_bams; genotype_bams }

log.info("Alignments: ${params.bams}")

process smoove_call {
    tag "sample: $sample"
    publishDir "$outdir", mode: "copy"
    cpus 1

    input:
    set sample, file(bam), file(bai) from call_bams
    file fasta
    file faidx

    output:
    file("${sample}-smoove.genotyped.vcf.gz") into vcfs
    file("${sample}-smoove.genotyped.vcf.gz.csi") into idxs

    script:
    exclude = params.bed ? "--exclude ${params.bed}" : ''
    """
    smoove call --genotype --name $sample -p ${task.cpus} --fasta $fasta $exclude $bam
    """
}

process smoove_merge {
    publishDir "$outdir", mode: "copy"
    cache 'deep'

    input:
    file vcf from vcfs.collect()
    file idx from idxs.collect()
    file fasta
    file faidx

    output:
    file("merged.sites.vcf.gz") into sites

    script:
    """
    smoove merge --name merged -f $fasta $vcf
    """
}

process smoove_genotype {
    tag "sample: $sample"
    publishDir "$outdir", mode: "copy"
    cpus 1

    input:
    set sample, file(bam), file(bai) from genotype_bams
    file sites
    file fasta
    file faidx

    output:
    set sample, file("${sample}-joint-smoove.genotyped.vcf.gz.csi")
    set sample, file("${sample}-joint-smoove.genotyped.vcf.gz")

    script:
    """
    smoove genotype -d -p ${task.cpus} -o ./ --name ${sample}-joint --fasta $fasta --vcf $sites $bam
    """
}
