import java.nio.file.Paths

params.bed = false
params.fasta = false
params.bams = false
params.outdir = false
params.excludechroms = false
params.project = false

gff = 'ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.chr.gff3.gz'
project = params.project ?: 'smoove-project'
outdir = file(Paths(params.outdir, project))

if( !params.fasta ) {
    exit 1, "No reference fasta was supplied"
}
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
    faidx = file(fasta + ".fai")
    if( !faidx.exists() ) exit 1, "Fasta index file not found: ${params.fasta}.fai"
    log.info("Reference fasta: ${fasta}")
}
if ( params.bed ){
    bed_file = file(params.bed)
    if ( !bed_file.exists() ) exit 1, "Bed file not found: ${params.bed}"
    log.info("Excluded regions: ${params.bed}")
}
log.info("Alignments: ${params.bams}")
log.info("Output: ${outdir}")

Channel
    .fromPath(params.bams, checkIfExists: true)
    .map { file -> tuple(file.baseName.split("\\.")[0], file, file + (file.endsWith('.cram') ? '.crai' : '.bai')) }
    .into { call_bams; genotype_bams }

process smoove_call {
    tag "sample: $sample"
    publishDir path: "$outdir/smoove-called", mode: "copy"
    cpus 1

    input:
    set sample, file(bam), file(bai) from call_bams
    file fasta
    file faidx
    file bed_file

    output:
    file("${sample}-smoove.genotyped.vcf.gz") into vcfs
    file("${sample}-smoove.genotyped.vcf.gz.csi") into idxs

    script:
    excludechroms = params.excludechroms ? "--excludechroms \"${params.excludechroms}\"" : ''
    """
    smoove call --genotype --name $sample --processes ${task.cpus} --fasta $fasta --exclude $bed_file $excludechroms $bam
    """
}

process smoove_merge {
    publishDir path: "$outdir/smoove-merged", mode: "copy"
    cache 'deep'

    input:
    file vcf from vcfs.collect()
    file idx from idxs.collect()
    file fasta
    file faidx

    output:
    file("${project}.sites.vcf.gz") into sites

    script:
    """
    smoove merge --name $project --fasta $fasta $vcf
    """
}

process smoove_genotype {
    tag "sample: $sample"
    publishDir path: "$outdir/smoove-genotyped", mode: "copy"
    cpus 1

    input:
    set sample, file(bam), file(bai) from genotype_bams
    file sites
    file fasta
    file faidx

    output:
    file("${sample}-smoove.genotyped.vcf.gz.csi") into genotyped_idxs
    file("${sample}-smoove.genotyped.vcf.gz") into genotyped_vcfs

    script:
    """
    wget -q https://raw.githubusercontent.com/samtools/samtools/develop/misc/seq_cache_populate.pl
    perl seq_cache_populate.pl -root \$(pwd)/cache $fasta 1> /dev/null 2> err || (cat err; exit 2)
    export REF_PATH=\$(pwd)/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
    export REF_CACHE=xx

    samtools quickcheck -v $bam
    smoove genotype --duphold --processes ${task.cpus} --removepr --outdir ./ --name ${sample} --fasta $fasta --vcf $sites $bam
    """
}

process smoove_square {
    publishDir path: "$outdir/smoove-squared", mode: "copy"
    cache 'deep'
    cpus 3

    input:
    file vcf from genotyped_vcfs.collect()
    file idx from genotyped_idxs.collect()

    output:
    file("square.anno.vcf.gz") into square_vcf
    file("square.anno.vcf.gz.csi") into square_idx

    script:
    """
    smoove paste --outdir ./ --name $project $vcf

    wget -q $gff
    smoove annotate --gff ${gff.split("\\/")[-1]} ${project}.smoove.square.vcf.gz | bgzip --threads ${task.cpus} -c > $square_vcf
    bcftools index $square_vcf
    """
}
