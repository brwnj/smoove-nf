// PATH=/scratch/ucgd/lustre/u6000771/Data/gemini_install2/tools/bin:/scratch/ucgd/lustre/u6000771/Data/gemini_install2/data/anaconda/bin/:$PATH:~u6000771/bin
// PATH=~/Projects/src/lumpy-sv/bin:$PATH

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
}
if ( params.bed ){
    bed_file = file(params.bed)
    if ( !bed_file.exists() ) exit 1, "Bed file specified [${params.bed}], but does not exist"
}


outdir = file(params.outdir)

Channel
    .fromPath(params.bams, checkIfExists: true)
    .map { file -> tuple(file.baseName, file) }
    .into { call_bams; genotype_bams }

process smoove_call {
    tag "smoove call: $sample"
    publishDir "$outdir", mode: "copy"

    input:
    set sample, file(bam) from call_bams
    file fasta

    output:
    set sample, file("${sample}-smoove.genotyped.vcf.gz") into vcfs
    set sample, file("${sample}-smoove.genotyped.vcf.gz.csi")

    script:
    exclude = params.bed ? "--exclude ${params.bed}" : ''
    """
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/scratch/ucgd/lustre/work/u6000771/Projects/src/htslib/
    smoove call --genotype --name $sample --fasta $fasta $exclude $bam
    """
}

process smoove_merge {
    publishDir "$outdir", mode: "copy"

    input:
    file vcfs
    file fasta

    output:
    file("merged.sites.vcf.gz") into sites

    script:
    """
    smoove merge --name merged -f $fasta $vcfs
    """
}

process smoove_genotype {
    tag "smoove genotype: $sample"
    publishDir "$outdir", mode: "copy"

    input:
    set sample, file(bam) from genotype_bams
    file sites

    output:
    set sample, file("${sample}-joint-smoove.genotyped.vcf.gz.csi")
    set sample, file("${sample}-joint-smoove.genotyped.vcf.gz")

    script:
    """
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/scratch/ucgd/lustre/work/u6000771/Projects/src/htslib/
    smoove genotype -d -p ${task.cpus} -o ./ --name ${sample}-joint --fasta $fasta --vcf $sites $bam
    """
}
