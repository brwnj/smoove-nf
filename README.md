# smoove-nf

Nextflow implementation of the smoove toolset.

## Workflow parameters

### `--bams`

Aligned sequences in .bam and/or .cram format. Indexes (.bai/.crai) must be
present.

```
nextflow run brwnj/smoove-nf --bams 'SRP1234/alignments/*.cram'
```

### `--project`

Acts as the file prefix for merged and squared sites.

* Default: 'sites'

```
nextflow run brwnj/smoove-nf --project SRP1234
```

### `--sexchroms`

Comma delimited names of the sex chromosome(s) used to infer sex.

* Default: 'X,Y'

```
nextflow run brwnj/smoove-nf --sexchroms 'chrX,chrY'
```

### `--fasta`

File path to reference fasta. Index (.fai) must be present.

```
nextflow run brwnj/smoove-nf --fasta /assets/g1k_v37_decoy/g1k_v37_decoy.fa
```

GRCh38 reference is available at:
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome

### `--bed`

File path to bed of exclude regions for `smoove call`.

```
nextflow run brwnj/smoove-nf --bed /assets/regions.exclude.bed.gz
```

Exclude regions for b37 and GRCh38 are available at:
https://github.com/hall-lab/speedseq/tree/master/annotations

### `--outdir`

The base results directory for output.

```
nextflow run brwnj/smoove-nf --outdir /scratch/projects/SRP1234
```

### `--excludechroms`

Chromosomes to exclude during `smoove call`. You should escape '$'.

```
nextflow run brwnj/smoove-nf --excludechroms "hs37d5,~:,~^GL,~decoy,~random\$,~chrUn,~_alt\$"
```

### `--gff`

Annotation GFF used in `smoove annotate`.

GRCh38 reference is available via Ensembl at:
ftp://ftp.ensembl.org/pub/release-95/gff3/homo_sapiens/Homo_sapiens.GRCh38.95.chr.gff3.gz
