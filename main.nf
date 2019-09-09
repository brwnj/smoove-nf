params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------

    smoove-nf: a smoove workflow
    ============================

    Documentation and issues can be found at: https://github.com/brwnj/smoove-nf

    smoove is available at: https://github.com/brentp/smoove

    Required arguments:
    -------------------

    --bams                Aligned sequences in .bam and/or .cram format.
                          Indexes (.bai/.crai) must be present.
    --bed                 Bed of exclude regions for `smoove call`.
    --fasta               Reference FASTA. Index (.fai) must exist in same
                          directory.
    --gff                 Annotation GFF used to annotate variants.

    Options:
    --------

    --outdir              Base results directory for output.
                          Default: '/.results'
    --exclude             Regular expression of chromosomes to skip.
                          Default: "~^HLA,~^hs,~:,~^GL,~M,~EBV,~^NC,~^phix,~decoy,~random\$,~Un,~hap,~_alt\$"
    --project             File prefix for merged and annotated VCF files.
                          Default: 'sites'
    --sexchroms           Comma delimited names of the sex chromosome(s)
                          used to infer sex. Default: 'X,Y'
    --sensitive           Preserves more variants from being filtered
                          throughout the workflow. Default: false

    covviz options:
    ---------------

    --zthreshold          A sample must greater than this many standard
                          deviations in order to be found significant.
                          Default: 3.5
    --distancethreshold   Consecutive significant points must span this
                          distance in order to pass this filter. Default:
                          150000
    --slop                Leading and trailing segments added to
                          significant regions to make them more visible.
                          Default: 500000
    --minsamples          Show all traces when analyzing this few samples;
                          ignores z-threshold, distance-threshold, and
                          slop. Default: 8

    somalier options:
    -----------------

    --knownsites          VCF of known polymorphic sites. Download links
                          can be found at:
                          https://github.com/brentp/somalier/releases
                          Default: false
    --ped                 Sample relationship definitions. Default: false

    -----------------------------------------------------------------------
    """.stripIndent()
    exit 0
}


// required arguments
params.bed = false
if( !params.bed ) { exit 1, "--bed is not defined" }
params.fasta = false
if( !params.fasta ) { exit 1, "--fasta is not defined" }
params.bams = false
if( !params.bams ) { exit 1, "--bams is not defined" }
params.gff = false
if( !params.gff ) { exit 1, "--gff is not defined" }

// somalier
params.knownsites = false
params.ped = false

// variables
params.sensitive = false
project = params.project ?: 'sites'
sexchroms = params.sexchroms ?: 'X,Y'
sexchroms = sexchroms.replaceAll(" ", "")
outdir = params.outdir ?: './results'
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
if (params.knownsites) {
log.info("Known sites        (--knownsites)    : ${params.knownsites}")
}
if (params.ped) {
log.info("Pedigree file      (--ped)           : ${params.ped}")
}
log.info("Sensitive          (--sensitive)     : ${params.sensitive}")
log.info("Output             (--outdir)        : ${outdir}")
log.info("\n")

// instantiate files
fasta = file(params.fasta)
faidx = file("${params.fasta}.fai")
bed = file(params.bed)
gff = file(params.gff)

// look for somalier reference VCF
knownsites_file = false
if (params.knownsites) {
    knownsites_file = file(params.knownsites)
    // check file existence
    if (!knownsites_file.exists()) {
        exit 1, "Missing optional known sites file: ${knownsites_file}"
    }
}
custom_ped = false
if (params.ped) {
    custom_ped = file(params.ped)
    if (!custom_ped.exists()) {
        exit 1, "Missing optional ped file: ${custom_ped}"
    }
}

// check file existence
if( !fasta.exists() ) { exit 1, "Missing reference fasta: ${fasta}" }
if( !faidx.exists() ) { exit 1, "Missing reference fasta index: ${faidx}" }
if( !bed.exists() ) { exit 1, "Missing exclude regions: ${bed}" }
if( !gff.exists() ) { exit 1, "Missing annotations: ${gff}" }


Channel
    .fromPath(params.bams, checkIfExists: true)
    .map { file -> tuple(file.baseName, file, file + ("${file}".endsWith('.cram') ? '.crai' : '.bai')) }
    .into { call_bams; genotype_bams; somalier_bams }

Channel
    .fromPath(indexes, checkIfExists: true)
    .set { index_ch }

Channel
    .fromPath(params.bams)
    .map { params.sensitive ? "KEEP" : "FALSE" }
    .into { sensitive_call_ch; sensitive_genotype_ch }


process smoove_call {
    publishDir path: "$outdir/smoove/called", mode: "copy", pattern: "*.vcf.gz*"
    publishDir path: "$outdir/logs", mode: "copy", pattern: "*-stats.txt"
    publishDir path: "$outdir/logs", mode: "copy", pattern: "*-smoove-call.log"

    input:
    env SMOOVE_KEEP_ALL from sensitive_call_ch
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
    excludechroms = params.exclude ? "--excludechroms \"${params.exclude}\"" : ""
    filters = params.sensitive ? "--noextrafilters" : ""
    """
    smoove call --genotype --name $sample --processes ${task.cpus} \
        --fasta $fasta --exclude $bed $excludechroms $filters \
        $bam 2> ${sample}-smoove-call.log
    bcftools stats ${sample}-smoove.genotyped.vcf.gz > ${sample}-stats.txt
    """
}

process smoove_merge {
    // publishDir path: "$outdir/smoove/merged", mode: "copy"

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
    publishDir path: "$outdir/smoove/genotyped", mode: "copy"

    input:
    env SMOOVE_KEEP_ALL from sensitive_genotype_ch
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
    publishDir path: "$outdir/smoove/annotated", mode: "copy", pattern: "*.vcf.gz*"
    publishDir path: "$outdir/reports/bpbio", mode: "copy", pattern: "*.html"

    input:
    file vcf from genotyped_vcfs.collect()
    file idx from genotyped_idxs.collect()
    file gff

    output:
    file("${project}.smoove.square.anno.vcf.gz") into square_vcf
    file("${project}.smoove.square.anno.vcf.gz.csi") into square_idx
    file("svvcf.html") into svvcf

    script:
    smoovepaste = "smoove paste --outdir ./ --name $project $vcf"
    if( vcf.collect().size() < 2 ) {
        paste = "cp $vcf ${project}.smoove.square.vcf.gz && cp $idx ${project}.smoove.square.vcf.gz.csi"
    }
    """
    $smoovepaste

    smoove annotate --gff $gff ${project}.smoove.square.vcf.gz | bgzip --threads ${task.cpus} -c > ${project}.smoove.square.anno.vcf.gz
    bcftools index ${project}.smoove.square.anno.vcf.gz
    bpbio plot-sv-vcf ${project}.smoove.square.anno.vcf.gz
    """
}

process run_indexcov {
    publishDir path: "$outdir/reports/indexcov", mode: "copy"

    input:
    file idx from index_ch.collect()
    file faidx

    output:
    file("${project}*.png")
    file("*.html")
    file("${project}*.bed.gz") into bed_ch
    file("${project}*.ped") into indexcov_ped_ch
    file("${project}*.roc") into roc_ch

    script:
    excludepatt = params.exclude ? "--excludepatt \"${params.exclude}\"" : ""
    """
    goleft indexcov --sex $sexchroms $excludepatt --directory $project --fai $faidx $idx
    mv $project/* .
    """
}

indexcov_ped_ch.into { ped_ch; report_ped_ch }

// account for optional, custom ped and the need to merge that with indexcov output
(merge_ch, report_ch) = (params.ped ? [ped_ch, Channel.empty()]: [Channel.empty(), ped_ch])

process merge_peds {
    label 'covviz'

    input:
    file ped from merge_ch
    file custom_ped

    output:
    file 'merged.ped' into merged_ch

    script:
    template 'merge_peds.py'
}

process build_covviz_report {
    publishDir path: "$outdir/reports", mode: "copy", pattern: "*.html"
    label 'covviz'
    cache 'lenient'

    input:
    file ped from report_ch.mix(merged_ch).collect()
    file bed from bed_ch
    file gff

    output:
    file("covviz_report.html")

    script:
    """
    covviz --min-samples ${params.minsamples} --sex-chroms ${params.sexchroms} --exclude '${params.exclude}' \
        --skip-norm --z-threshold ${params.zthreshold} --distance-threshold ${params.distancethreshold} \
        --slop ${params.slop} --ped ${ped} --gff ${gff} ${bed}
    """
}

process somalier_extract {
    // requires $sample to match name defined in @RG -- https://github.com/brentp/somalier/blob/master/src/somalier.nim#L18
    label 'somalier'
    publishDir path: "$outdir/somalier/extract", mode: "copy"

    input:
    set sample, file(bam), file(bai) from somalier_bams
    file knownsites_file
    file fasta
    file faidx

    output:
    file("${sample}.somalier") into somalier_counts

    // can be run even if user does not specify a ped file for `somalier relate`
    when: params.knownsites != false

    script:
    """
    somalier extract --out-dir ./ --fasta $fasta --sites $knownsites_file $bam
    """
}

process somalier_relate {
    label 'somalier'
    publishDir path: "$outdir/somalier", mode: "copy"

    input:
    file somalier_count from somalier_counts.collect()
    file custom_ped

    output:
    file("somalier.pairs.tsv")
    file("somalier.samples.tsv")
    file("somalier.html")

    when: (params.knownsites != false && params.ped != false)

    script:
    """
    somalier relate --ped $custom_ped $somalier_count
    """
}

process build_report {
    publishDir path: "$outdir/reports", mode: "copy", pattern: "*.html", overwrite: true
    cache false

    input:
    file sequence_count from sequence_counts.collect()
    file variant_count from variant_counts.collect()
    file vcf from square_vcf
    file pedfile from report_ped_ch
    file variant_html from svvcf

    output:
    file("smoove-nf.html")

    script:
    template 'smoove-report.py'
}
