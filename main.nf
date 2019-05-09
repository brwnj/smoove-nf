// required arguments
params.bed = false
if( !params.bed ) { exit 1, "--bed is not defined" }
params.fasta = false
if( !params.fasta ) { exit 1, "--fasta is not defined" }
params.bams = false
if( !params.bams ) { exit 1, "--bams is not defined" }
params.outdir = false
if( !params.outdir ) { exit 1, "--outdir is not defined" }
params.gff = false
if( !params.gff ) { exit 1, "--gff is not defined" }

// variables
project = params.project ?: 'sites'
sexchroms = params.sexchroms ?: 'X,Y'
sexchroms = sexchroms.replaceAll(" ", "")
outdir = params.outdir
indexes = params.bams + ("${params.bams}".endsWith('.cram') ? '.crai' : '.bai')

log.info("\n")
log.info("Project            (--project)       : ${project}")
log.info("Excluded regions   (--bed)           : ${params.bed}")
if( params.exclude ) {
log.info("Excluded chroms    (--exclude)       : ${params.exclude}")
}
log.info("Reference fasta    (--fasta)         : ${params.fasta}")
log.info("Sex chromosomes    (--sexchroms)     : ${sexchroms}")
log.info("Alignments         (--bams)          : ${params.bams}")
log.info("Indexes                              : ${indexes}")
log.info("Annotation GFF     (--gff)           : ${params.gff}")
log.info("Output             (--outdir)        : ${outdir}")
log.info("\n")

// instantiate files
fasta = file(params.fasta)
faidx = file("${params.fasta}.fai")
bed = file(params.bed)
gff = file(params.gff)

// check file existence
if( !fasta.exists() ) { exit 1, "Missing reference fasta: ${fasta}" }
if( !faidx.exists() ) { exit 1, "Missing reference fasta index: ${faidx}" }
if( !bed.exists() ) { exit 1, "Missing exclude regions: ${bed}" }
if( !gff.exists() ) { exit 1, "Missing annotations: ${gff}" }


Channel
    .fromPath(params.bams, checkIfExists: true)
    .map { file -> tuple(file.baseName, file, file + ("${file}".endsWith('.cram') ? '.crai' : '.bai')) }
    .into { call_bams; genotype_bams }

Channel
    .fromPath(indexes, checkIfExists: true)
    .set { index_ch }

process smoove_call {
    tag "sample: $sample"
    publishDir path: "$outdir/smoove-called", mode: "copy", pattern: "*.vcf.gz*"
    publishDir path: "$outdir/logs", mode: "copy", pattern: "*-stats.txt"
    publishDir path: "$outdir/logs", mode: "copy", pattern: "*-smoove-call.log"
    memory { 16.GB * task.attempt }
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    cache 'lenient'

    input:
    set sample, file(bam), file(bai) from call_bams
    file fasta
    file faidx
    file bed

    output:
    file("${sample}-smoove.genotyped.vcf.gz") into vcfs
    file("${sample}-smoove.genotyped.vcf.gz.csi") into idxs
    file("${sample}-stats.txt") into variant_counts
    file("${sample}-smoove-call.log") into sequence_counts

    script:
    excludechroms = params.exclude ? "--excludechroms \"${params.exclude}\"" : ''
    """
    smoove call --genotype --name $sample --processes ${task.cpus} --fasta $fasta --exclude $bed $excludechroms $bam 2> ${sample}-smoove-call.log
    bcftools stats ${sample}-smoove.genotyped.vcf.gz > ${sample}-stats.txt
    """
}

process smoove_merge {
    publishDir path: "$outdir/smoove-merged", mode: "copy"
    memory 16.GB
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
    errorStrategy { task.attempt < 3 ? 'retry' : 'terminate' }
    memory 16.GB
    cache 'lenient'

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
    publishDir path: "$outdir/smoove-squared", mode: "copy", pattern: "*.vcf.gz*"
    publishDir path: "$outdir/bpbio", mode: "copy", pattern: "*.html"
    cpus 3
    memory 64.GB
    cache 'deep'

    input:
    file vcf from genotyped_vcfs.collect()
    file idx from genotyped_idxs.collect()
    file gff

    output:
    file("${project}.smoove.square.anno.vcf.gz") into square_vcf
    file("${project}.smoove.square.anno.vcf.gz.csi") into square_idx
    file("svvcf.html") into svvcf

    script:
    """
    smoove paste --outdir ./ --name $project $vcf

    smoove annotate --gff $gff ${project}.smoove.square.vcf.gz | bgzip --threads ${task.cpus} -c > ${project}.smoove.square.anno.vcf.gz
    bcftools index ${project}.smoove.square.anno.vcf.gz
    bpbio plot-sv-vcf ${project}.smoove.square.anno.vcf.gz
    """
}

process run_indexcov {
    publishDir path: "$outdir/indexcov", mode: "copy"
    memory 8.GB
    cache 'deep'

    input:
    file idx from index_ch.collect()
    file faidx

    output:
    file("${project}*.png")
    file("*.html")
    file("${project}*.bed.gz") into coverage_bed
    file("${project}*.ped") into ped
    file("${project}*.roc") into roc

    script:
    excludepatt = params.exclude ? "--excludepatt \"${params.exclude}\"" : ''
    """
    goleft indexcov --sex $sexchroms $excludepatt --directory $project --fai $faidx $idx
    mv $project/* .
    """
}

process build_covviz_report {
    publishDir path: "$outdir", mode: "copy", pattern: "*.html"
    container 'jupyter/scipy-notebook:7d427e7a4dde'
    cache 'lenient'

    input:
    file pedfile from ped
    file rocfile from roc
    file bedfile from coverage_bed
    file gff

    output:
    file("covviz_report.html")

    script:
    template 'parse_indexcov.py'
}

process build_report {
    publishDir path: "$outdir", mode: "copy", pattern: "*.html", overwrite: true
    cache false

    input:
    file sequence_count from sequence_counts.collect()
    file variant_count from variant_counts.collect()
    file vcf from square_vcf
    file pedfile from ped
    file variant_html from svvcf

    output:
    file("smoove-nf.html")

    script:
    template 'smoove-report.py'
}
