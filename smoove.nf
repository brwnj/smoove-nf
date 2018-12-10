params.bed = false
params.fasta = false
params.bams = false
params.outdir = false
params.excludechroms = false
params.project = false
params.gff = false

project = params.project ?: 'sites'
outdir = params.outdir
fasta = file(params.fasta)
faidx = file("${params.fasta}.fai")
bed = file(params.bed)
gff = file(params.gff)
indexes = params.bams + ("${params.bams}".endsWith('.cram') ? '.crai' : '.bai')

log.info("\n")
log.info("Project: ${project}")
log.info("Excluded regions: ${params.bed}")
if( params.excludechroms ) log.info("Excluded chroms: ${params.excludechroms}")
log.info("Reference fasta: ${params.fasta}")
log.info("Alignments: ${params.bams}")
log.info("Annotation GFF: ${params.gff}")
log.info("Output: ${outdir}")
log.info("\n")

Channel
    .fromPath(params.bams, checkIfExists: true)
    .map { file -> tuple(file.baseName.split("\\.")[0], file, file + ("${file}".endsWith('.cram') ? '.crai' : '.bai')) }
    .into { flagstat_bams; call_bams; genotype_bams }

Channel
    .fromPath(indexes, checkIfExists: true)
    .set { index_ch }

process run_flagstat {
    tag "sample: $sample"
    publishDir path: "$outdir/logs", mode: "copy"

    input:
    set sample, file(bam), file(bai) from flagstat_bams

    output:
    file("${sample}-flagstat.txt") into sequence_counts

    script:
    """
    samtools flagstat $bam > $sample-flagstat.txt
    """
}

process smoove_call {
    tag "sample: $sample"
    publishDir path: "$outdir/smoove-called", mode: "copy"
    memory { 8.GB * task.attempt }
    errorStrategy { task.attempt == 1 ? 'retry' : 'ignore' }

    input:
    set sample, file(bam), file(bai) from call_bams
    file fasta
    file faidx
    file bed

    output:
    file("${sample}-smoove.genotyped.vcf.gz") into vcfs
    file("${sample}-smoove.genotyped.vcf.gz.csi") into idxs
    file("${sample}-stats.txt") into variant_counts

    script:
    excludechroms = params.excludechroms ? "--excludechroms \"${params.excludechroms}\"" : ''
    """
    smoove call --genotype --name $sample --processes ${task.cpus} --fasta $fasta --exclude $bed $excludechroms $bam
    bcftools stats ${sample}-smoove.genotyped.vcf.gz > ${sample}-stats.txt
    """
}

process smoove_merge {
    publishDir path: "$outdir/smoove-merged", mode: "copy"
    cache 'deep'
    memory 15.GB

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
    memory 45.GB

    input:
    file vcf from genotyped_vcfs.collect()
    file idx from genotyped_idxs.collect()
    file gff

    output:
    file("${project}.smoove.square.anno.vcf.gz") into square_vcf
    file("${project}.smoove.square.anno.vcf.gz.csi") into square_idx

    script:
    """
    smoove paste --outdir ./ --name $project $vcf

    smoove annotate --gff $gff ${project}.smoove.square.vcf.gz | bgzip --threads ${task.cpus} -c > ${project}.smoove.square.anno.vcf.gz
    bcftools index ${project}.smoove.square.anno.vcf.gz
    """
}

process run_indexcov {
    publishDir path: "$outdir/indexcov", mode: "copy"

    input:
    file idx from index_ch.collect()
    file faidx

    output:
    file("${project}*.png")
    file("${project}*.html")
    file("${project}*.bed.gz")
    file("${project}*.ped") into peds
    file("${project}*.roc")

    script:
    excludepatt = params.excludechroms ? "--excludepatt \"${params.excludechroms}\"" : ''
    """
    curl --location 'https://github.com/brentp/goleft/releases/download/v0.2.0/goleft_linux64' > goleft
    chmod +x goleft
    ./goleft indexcov $excludepatt --directory $project --fai $faidx $idx
    mv $project/* .
    """
}

process build_report {
    publishDir path: "$outdir", mode: "copy"

    input:
    file sequence_count from sequence_counts.collect()
    file variant_count from variant_counts.collect()
    file vcf from square_vcf
    file ped from peds

    output:
    file("smoove-nf.html")

    script:
    template 'smoove-report.py'
}
