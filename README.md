# smoove-nf

Nextflow implementation of the [smoove](https://github.com/brentp/smoove) toolset (and some others) focused on reliably calling SVs in your data.

## The workflow

The workflow consists of a number of steps, each generally outputing to unique result directories.

#### Call genotypes

`smoove call` is run on individual bam or cram alignment files. Output is written to `$outdir/smoove-called` and includes `$sample-smoove.genotyped.vcf.gz` and an index.

#### Merge genotypes

Next, we collect all SVs across samples into a single, merged (union) VCF using `smoove merge`. Results are written to `$outdir/smoove-merged` and include the file `$project.sites.vcf.gz`.

#### Genotype all samples

Using the union of SVs across all samples, we genotype each sample at those sites using `smoove genotype` with `duphold` for depth annotations. Output is written to `$outdir/smoove-genotyped/$sample-smoove.genotyped.vcf.gz`.

#### Square and annotate VCF

Take all single sample genotyped VCFs and paste into a single, square, joint-called file using `smoove paste`. Then annotate the variants using the annotation supplied from `--gff` with `smoove annotate`. Results are written to:

+ `$outdir/smoove-squared/$project.smoove.square.anno.vcf.gz`
	+ Annotated and indexed VCF for all SVs across all samples.
+ `$outdir/bpbio/svvcf.html`
	+ A report of SV counts per sample by SV type.

#### Coverage profiling

Using [indexcov](https://github.com/brentp/goleft/tree/master/indexcov), estimate coverage across the genome per sample and perform coverage-based quality control. The full report output of `goleft indexcov` is written to `$outdir/indexcov`. Its report is written to `$outdir/indexcov/index.html`.

#### Workflow report

Logs and output of various steps are aggregated and summarized into one report written to `$outdir/smoove-nf.html`.

Cumulative chromosome coverage is available in `$outdir/covviz_report.html`.

## Usage

A Docker container is maintained in parallel with this workflow (https://hub.docker.com/r/brentp/smoove) and will be pulled by Nextflow before data processing begins. There's no need to download and install dependencies outside of Docker or Singularity and [Nextflow](https://www.nextflow.io/).

```
nextflow run brwnj/smoove-nf -latest [nextflow options] [smoove-nf options]
```

### Required parameters

+ `--bams`
	+ Aligned sequences in .bam and/or .cram format. Indexes (.bai/.crai) must be
present.
	+ Use wildcards like `'SRP1234/alignments/*.cram'` to specify your alignment files.

Note: Nextflow will handle wildcard expansion in this case, so it's critical we quote we the string like:

```
nextflow run brwnj/smoove-nf -latest \
	--bams '~/SRP1234/alignments/*.cram'
```

+ `--bed`
	+ File path to bed of exclude regions for `smoove call`.
	+ Exclude regions for b37 and GRCh38 are made available by the Hall lab under [speedseq](https://github.com/hall-lab/speedseq/tree/master/annotations).
+ `--fasta`
	+ File path to reference fasta. Index (.fai) must be present.
	+ GRCh38 at [1k Genomes](ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome)
+ `--gff`
	+ Annotation GFF used in `smoove annotate`.
	+ GRCh38 reference is available via [Ensembl](ftp://ftp.ensembl.org/pub/release-95/gff3/homo_sapiens/Homo_sapiens.GRCh38.95.chr.gff3.gz)

### Optional parameters

+ `--outdir`
	+ The base results directory for output
	+ **default:** './results'
+ `--exclude`
    + regular expression of chromosomes to skip
	+ You should escape '$', e.g. `"hs37d5,~:,~^GL,~decoy,~random\$,~chrUn,~_alt\$"`
    + **default**: "^GL|^hs|^chrEBV$|M$|MT$|^NC|_random$|Un_|^HLA\\-|_alt$|hap\\d+$"
+ `--project`
	+ Acts as the file prefix for merged and squared sites
	+ **default:** 'sites'
+ `--sexchroms`
	+ Comma delimited names of the sex chromosome(s) used to infer sex, e.g. `--sexchroms 'chrX,chrY'`
	+ **default:** 'X,Y'

#### [covviz](https://github.com/brwnj/covviz) params
+ `--zthreshold`
    + a sample must greater than this many standard deviations in order to be found significant
    + **default:** 3.5
+ `--distancethreshold`
    + consecutive significant points must span this distance in order to pass this filter
    + **default:** 150000
+ `--slop`
    + leading and trailing segments added to significant regions to make them more visible
    + **default:** 500000

#### [somalier](https://github.com/brentp/somalier) params
+ `--knownsites`
	+ optional, but required in order to run [somalier](https://github.com/brentp/somalier) quality control
	+ VCF of known polymorphic sites -- download links can be found at https://github.com/brentp/somalier/releases, but any set of common variants will work
	+ **default:** false
+ `--ped`
	+ optional, but required in order to run `somalier relate` and generate somalier's HTML report
	+ sample relationship definitions
	+ **default:** false


## Updating

To pull changes to made to the workflow and ensure you're running the latest version, use:

```
nextflow pull brwnj/smoove-nf
```

That will either pull any changes or confirm you're at the latest version.
