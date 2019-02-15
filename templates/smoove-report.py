#!/usr/bin/env python
from __future__ import print_function

import csv
import gzip
import json
import os

from collections import defaultdict


html = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8" />
    <meta name="author" content="Joe Brown" />
    <title>smoove-nf</title>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.10.19/css/dataTables.bootstrap4.min.css">
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css">
    <script src="https://code.jquery.com/jquery-3.3.1.js"></script>
    <script src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.19/js/dataTables.bootstrap4.min.js"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/js/bootstrap.min.js"></script>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style type="text/css">
    .container {max-width: 90%}
    h1:before {
        height: 70px;
        content: "";
        display:block;
    }
    body {
        position: relative;
    }
	/* fix asc/desc arrow positioning */
	table.dataTable.table-sm .sorting:before, table.dataTable.table-sm .sorting_asc:before, table.dataTable.table-sm .sorting_desc:before {
	    right: 1.15em;
	}
    </style>
</head>
<body>
    <nav id="main_nav" class="navbar navbar-expand-md navbar-dark fixed-top bg-dark">
        <span class="navbar-brand mb-0 h1">smoove-nf report</span>
        <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarCollapse" aria-controls="navbarCollapse" aria-expanded="false" aria-label="Toggle navigation">
          <span class="navbar-toggler-icon"></span>
        </button>
        <div class="collapse navbar-collapse" id="navbarCollapse">
          <ul class="navbar-nav mr-auto">
                <li class="nav-item"><a class="nav-link" href="#processing">Samples</a></li>
                <li class="nav-item"><a class="nav-link" href="#variants">Calls</a></li>
                <li class="nav-item"><a class="nav-link" href="#coverage">Coverage</a></li>
                <li class="nav-item"><a class="nav-link" href="#filtering">Filtering</a></li>
                <li class="nav-item"><a class="nav-link" href="#configuration">Configuration</a></li>
            </ul>
        </div>
    </nav>

    <div data-spy="scroll" data-target="#main_nav" data-offset="70" class="container pb-20">
        <h1 class="border-bottom border-dark" id="processing">Summary table</h1>
        <p>Read and variant counts are obtained as samples are processed by <a href="https://github.com/arq5x/lumpy-sv">LUMPY</a>
           via <a href="https://github.com/brentp/smoove">smoove</a>. Samples where variants fail to be called tend
           to have a higher number of variants suggesting issues with the sample. Failed samples are not included
           in the merged VCF, but are genotyped and included in the annotated VCF.</p>
		   <table id="sample_table" class="table table-hover table-sm" width="100%"></table>
        <script>
        \$('body').scrollspy({ target: '#main_nav' })
        var success = '<span class="badge badge-success">Success</span>'
        var fail = '<span class="badge badge-danger">Fail</span>'
        var dataSet = [SAMPLE_SUMMARY];

        \$(document).ready(function() {
            \$('#sample_table').DataTable( {
                data: dataSet,
                columns: [
                    { title: "Sample" },
                    { title: "Reads", render: \$.fn.dataTable.render.number(",")},
                    { title: "Called Variants", render: \$.fn.dataTable.render.number(",") },
                    { title: "Called" },
                    { title: "Genotyped" },
                ]
            } );

            var selected_chrom = \$('#chrom_selector input:radio:checked').data('name')
            build_coverage_by_percent_plot(selected_chrom)
            build_coverage_by_position_plot(selected_chrom)
        } );
        </script>

        <h1 class="border-bottom border-dark" id="variants">SV call summary by sample</h1>
        <h3>Deletions</h3>
        <div class="row">
          <div id="deletions_p1" class="col-6"></div>
          <div id="deletions_p2" class="col-6"></div>
        </div>
        <h3>Duplications</h3>
        <div class="row">
          <div id="duplications_p1" class="col-6"></div>
          <div id="duplications_p2" class="col-6"></div>
        </div>
        <h3>Inversions</h3>
        <div class="row">
          <div id="inversions_p1" class="col-6"></div>
          <div id="inversions_p2" class="col-6"></div>
        </div>
        <h3>Break ends</h3>
        <div class="row">
          <div id="bnd_p1" class="col-6"></div>
          <div id="bnd_p2" class="col-6"></div>
        </div>
        <script>
        var var_samples = [VARIANT_SAMPLE_LIST]
        var variant_bar_layout = {
            height: 350,
            xaxis: {
                title: "Sample",
                showticklabels: false,
				type: 'category',
            },
            yaxis: {
                title: "",
                autorange: true
            },
            barmode: "stack",
            hovermode: "closest",
            legend: {
                orientation: "h",
                x: 0,
                y: 0,
            },
            margin: {t: 10},
        }
        var variant_hist_layout = {
            height: 350,
            xaxis: {
                title: "",
                showticklabels: false,
                autorange: true
            },
            yaxis: {
                title: "Count of Variants",
                autorange: true,
                log: true,
                type: "log"
            },
            hovermode: "closest",
            margin: {t: 10},
        }
        var deletions_p1_data = [
            {
                x: var_samples,
                y: [SMALL_DELETIONS],
                text: var_samples,
                hoverinfo: "text+x+y+name",
                mode: "lines",
                type: "bar",
                name: "small deletions"
            },{
                x: var_samples,
                y: [LARGE_DELETIONS],
                mode: "lines",
                type: "bar",
                name: "large deletions"
            }
        ]
        var deletions_p2_data = [
            {
                x: [DELETIONS_HIST],
                cumulative: {enabled:false},
                histfunc: "count",
                histnorm:"",
                mode: "lines",
                type: "histogram",
                marker: {color: "#7F7F7F"}
            }
        ]
        variant_bar_layout.yaxis.title = "Deletions"
        variant_hist_layout.xaxis.title = "# of Samples with Deletions"
        Plotly.newPlot('deletions_p1', deletions_p1_data, variant_bar_layout)
        Plotly.newPlot('deletions_p2', deletions_p2_data, variant_hist_layout)

        var duplications_p1_data = [
            {
                x: var_samples,
                y: [SMALL_DUPLICATIONS],
                text: var_samples,
                hoverinfo: "text+x+y+name",
                mode: "lines",
                type: "bar",
                name: "small duplications"
            },{
                x: var_samples,
                y: [LARGE_DUPLICATIONS],
                mode: "lines",
                type: "bar",
                name: "large duplications"
            }
        ]
        var duplications_p2_data = [
            {
                x: [DUPLICATIONS_HIST],
                cumulative: {enabled:false},
                histfunc: "count",
                histnorm:"",
                mode: "lines",
                type: "histogram",
                marker: {color: "#7F7F7F"}
            }
        ]
        variant_bar_layout.yaxis.title = "Duplications"
        variant_hist_layout.xaxis.title = "# of Samples with Duplications"
        Plotly.newPlot('duplications_p1', duplications_p1_data, variant_bar_layout)
        Plotly.newPlot('duplications_p2', duplications_p2_data, variant_hist_layout)

        var inversions_p1_data = [
            {
                x: var_samples,
                y: [SMALL_INVERSIONS],
                text: var_samples,
                hoverinfo: "text+x+y+name",
                mode: "lines",
                type: "bar",
                name: "small inversions"
            },{
                x: var_samples,
                y: [LARGE_INVERSIONS],
                mode: "lines",
                type: "bar",
                name: "large inversions"
            }
        ]
        var inversions_p2_data = [
            {
                x: [INVERSIONS_HIST],
                cumulative: {enabled:false},
                histfunc: "count",
                histnorm:"",
                mode: "lines",
                type: "histogram",
                marker: {color: "#7F7F7F"}
            }
        ]
        variant_bar_layout.yaxis.title = "Inversions"
        variant_hist_layout.xaxis.title = "# of Samples with Inversions"
        Plotly.newPlot('inversions_p1', inversions_p1_data, variant_bar_layout)
        Plotly.newPlot('inversions_p2', inversions_p2_data, variant_hist_layout)

        var bnd_p1_data = [
            {
                x: var_samples,
                y: [SMALL_BNDS],
                mode: "lines",
                type: "bar",
                name: "small BNDs",
                text: var_samples
            },{
                x: var_samples,
                y: [LARGE_BNDS],
                mode: "lines",
                type: "bar",
                name: "large BNDs",
                text: var_samples
            },{
                x: var_samples,
                y: [INTERCHROMOSOMAL_BNDS],
                mode: "lines",
                type: "bar",
                name: "interchromosomal BNDs",
                text: var_samples
            }
        ]
        var bnd_p2_data = [
            {
                x: [BNDS_HIST],
                cumulative: {enabled:false},
                histfunc: "count",
                histnorm:"",
                mode: "lines",
                type: "histogram",
                marker: {color: "#7F7F7F"}
            }
        ]
        variant_bar_layout.yaxis.title = "Break Ends"
        variant_hist_layout.xaxis.title = "# of Samples with Break Ends"
        Plotly.newPlot('bnd_p1', bnd_p1_data, variant_bar_layout)
        Plotly.newPlot('bnd_p2', bnd_p2_data, variant_hist_layout)
        </script>

        <h1 class="border-bottom border-dark" id="coverage">Sample coverage profiles</h1>
        <div class="row pb-5">
            <div class="col-6">
                <h3 id="scp_inferred_sex">Inferred sex</h3>
                By scaling coverage to a median of 1, <a href="https://github.com/brentp/goleft/tree/master/indexcov">indexcov</a> infers sex based on X and Y copy-number states. [<a href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-sex.md">docs</a>]
            </div>
            <div class="col-6">
                <h3 id="scp_bin_counts">Bin counts</h3>
                Ideally, sample depth per bin is near 1, indicating uniform coverage. Samples with a high number of low coverage bins fall to the right of the plot and samples towards the top of the plot have many regions outside of the expected coverage. [<a href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-bin.md">docs</a>]
            </div>
        </div>
        <div class="row">
            <div id="inferred_sex" class="col-6"></div>
            <div id="bin_counts" class="col-6"></div>
            <script>
            var samples = [COVERAGE_SAMPLE_LIST];
            var data1 = [{
                x: [INFERRED_X1],
                y: [INFERRED_Y1],
                mode: 'markers',
                type: 'scatter',
                name: 'Inferred CN for X: 1',
                text: [INFERRED_SAMPLES1],
                hoverinfo: 'text',
                marker: {
                  size: 12,
                  color: 'rgba(31,120,180,0.5)',
                  line: {
                      color: 'rgb(40,40,40)',
                      width: 1,
                  }
                },
            },{
                x: [INFERRED_X2],
                y: [INFERRED_Y2],
                mode: 'markers',
                type: 'scatter',
                name: 'Inferred CN for X: 2',
                text: [INFERRED_SAMPLES2],
                hoverinfo: 'text',
                marker: {
                    size: 12,
                    color: 'rgba(227,26,28,0.5)',
                    line: {
                        color: 'rgb(40,40,40)',
                        width: 1,
                    }
                },
            }]
            var layout = {
                title: "Inferred Sex",
                margin: {t: 35},
                height: 450,
                xaxis: {title: "X Copy Number"},
                yaxis: {title: "Y Copy Number"},
                hovermode: "closest",
                legend: {"orientation": "v"}
            }
            var plot0 = Plotly.newPlot('inferred_sex', data1, layout)

            layout.title = "Problematic low and non-uniform coverage bins"
            layout.xaxis.title = "Proportion of bins with depth < 0.15"
            layout.yaxis.title = "Proportion of bins with depth outside of (0.85, 1.15)"
            data = [{
                x: [BIN_X],
                y: [BIN_Y],
                mode: 'markers',
                type: 'scatter',
                name: 'Bins',
                text: samples,
                hoverinfo: 'text',
                marker: {
                    size: 12,
                    color: 'rgba(31,120,180,0.5)',
                    line: {
                        color: 'rgb(40,40,40)',
                        width: 1,
                    }
                },
            }]
            var plot1 = Plotly.newPlot('bin_counts', data, layout)
            </script>
        </div>

        [PCA_DIV]

        <h3>Cumulative chromosome coverage</h3>
        <p>Per chromosome coverage (ROC) plots show how much of the genome is covered at a given (scaled) depth. Shorter curves are indicative of increased low coverage regions. [<a href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-depth.md">docs</a>]</p>
        <div class="row">
            <div class="col-12 pt-3 pb-4">
                 <div id="chrom_selector" class="btn-group btn-group-toggle d-flex justify-content-center flex-wrap" data-toggle="buttons">
                      <div class="input-group-prepend">
                          <span class="input-group-text" id="basic-addon1">Chromosome:</span>
                      </div>
                      [CHROM_BUTTONS]
                </div>
            </div>
            <div class="col-12" id="cov_plot"></div>
            <div class="col-12" id="scaled_cov_plot"</div>
        </div>

        <script>
        var cov_color = 'rgba(108,117,125,0.3)'
        var cov_layout = {
            title: "",
            margin: {t: 10, b: 40},
            height: 350,
            xaxis: {title: "Scaled Coverage", showgrid: false, range: [0, 1.5]},
            yaxis: {title: "Proportion of Regions Covered", range: [0, 1.]},
            hovermode: "closest",
            showlegend: false,
        }
        var scaled_cov_layout = {
            title: "",
            margin: {t: 10, b: 40},
            height: 450,
            xaxis: {
                title: "Position",
                rangeslider: {},
                showgrid: false
            },
            yaxis: {title: "Scaled Coverage", fixedrange: true},
            hovermode: "closest",
            showlegend: false,
        }
        var cov_data = [CHROM_DATA]
        var cov_plot_data = cov_data[\$('#chrom_selector input:radio:checked').data('name')]
        var scaled_cov_plot_data;

        function build_coverage_by_percent_plot(chrom) {
            cov_plot_data = cov_data[chrom]
            var cov_plot = document.getElementById("cov_plot")
            Plotly.newPlot(cov_plot, cov_plot_data, cov_layout)
            cov_plot.on("plotly_click", coverage_plot_click)
        }

        function build_coverage_by_position_plot(chrom) {
            \$.getJSON("report_data/" + selected_chrom + ".json", function(json) {
                scaled_cov_plot_data = json
                var scaled_cov_plot = document.getElementById("scaled_cov_plot")
                Plotly.newPlot(scaled_cov_plot, scaled_cov_plot_data, scaled_cov_layout)
                scaled_cov_plot.on("plotly_click", coverage_plot_click)
            })
        }

        jQuery('#chrom_selector').on("change", function() {
            selected_chrom = \$('#chrom_selector input:radio:checked').data('name')
            build_coverage_by_percent_plot(selected_chrom)
            build_coverage_by_position_plot(selected_chrom)
        })

        function coverage_plot_click(click_data) {
            var local_plot_data = click_data.points[0].data
            var sample_id = local_plot_data.text

            for (var i = 0; i < scaled_cov_plot_data.length; i++) {
                // de-select traces
                scaled_cov_plot_data[i].marker.color = cov_color
                if (scaled_cov_plot_data[i].text == sample_id) {
                    scaled_cov_plot_data[i].marker.color = 'rgb(255,99,71)'
                }
            }

            for (var i = 0; i < cov_plot_data.length; i++) {
                // de-select traces
                cov_plot_data[i].marker.color = cov_color
                if (cov_plot_data[i].text == sample_id) {
                    // set color based on trace index
                    cov_plot_data[i].marker.color = 'rgb(255,99,71)'
                }
            }

            Plotly.redraw("cov_plot")
            Plotly.redraw("scaled_cov_plot")
        }
        </script>

        <h3>Full indexcov Results</h3>
        <p>INDEXCOV_RESULT</p>

        <h1 class="border-bottom border-dark" id="filtering">Alignment Filtering</h1>
        <p>Split and discordant read filtering results as reported by <code>smoove call</code>.</p>

        <div class="row">
            <div id="plot_before" class="col-6"></div>
            <div id="plot_after" class="col-6"></div>
            <script>
            var filter_samples = [FILTERED_SAMPLES]
            var data = [{
                x: [SPLIT_BEFORE],
                y: [DISCORDANT_BEFORE],
                mode: 'markers',
                type: 'scatter',
                name: 'Before Filtering',
                text: filter_samples,
                hoverinfo: 'text',
                marker: {
                    size: 12,
                    color: 'rgba(31,120,180,0.5)',
                    line: {
                        color: 'rgb(40,40,40)',
                        width: 1,
                    }
                },
            }]
            var layout = {title: "Before", margin:{t: 35}, height: 450, xaxis:{title:"Split Reads"}, yaxis:{"title":"Discordant Reads"},"hovermode":"closest"}
            var plot0 = Plotly.newPlot('plot_before', data, layout)

            data = [{
                x: [SPLIT_AFTER],
                y: [DISCORDANT_AFTER],
                mode: 'markers',
                type: 'scatter',
                name: 'After Filtering',
                text: filter_samples,
                hoverinfo: 'text',
                marker: {
                    size: 12,
                    color: 'rgba(31,120,180,0.5)',
                    line: {
                        color: 'rgb(40,40,40)',
                        width: 1,
                    }
                },
            }]
            layout.title = "After"
            var plot1 = Plotly.newPlot('plot_after', data, layout)
            </script>
        </div>

        <h1 class="border-bottom border-dark" id="configuration">Configuration</h1>
        <h3>Parameters</h3>
        <dl class="row small">
            <dt class="col-sm-3"><code>--bams</code></dt>
            <dd class="col-sm-9">${params.bams}</dd>
            <dt class="col-sm-3"><code>--outdir</code></dt>
            <dd class="col-sm-9">${params.outdir}</dd>
            <dt class="col-sm-3"><code>--fasta</code></dt>
            <dd class="col-sm-9">${params.fasta}</dd>
            <dt class="col-sm-3"><code>--bed</code></dt>
            <dd class="col-sm-9">${params.bed}</dd>
            <dt class="col-sm-3"><code>--excludechroms</code></dt>
            <dd class="col-sm-9">${params.excludechroms}</dd>
            <dt class="col-sm-3"><code>--project</code></dt>
            <dd class="col-sm-9">${params.project}</dd>
            <dt class="col-sm-3"><code>--gff</code></dt>
            <dd class="col-sm-9">${params.gff}</dd>
        </dl>

        <h3>Workflow</h3>
        <dl class="row small">
            <dt class="col-sm-3">Repository</dt>
            <dd class="col-sm-9">$workflow.repository</dd>
            <dt class="col-sm-3">Revision</dt>
            <dd class="col-sm-9">$workflow.revision</dd>
            <dt class="col-sm-3">Launch dir</dt>
            <dd class="col-sm-9">$workflow.launchDir</dd>
            <dt class="col-sm-3">Work dir</dt>
            <dd class="col-sm-9">$workflow.workDir</dd>
            <dt class="col-sm-3">Config files</dt>
            <dd class="col-sm-9">$workflow.configFiles</dd>
            <dt class="col-sm-3">Container</dt>
            <dd class="col-sm-9">$workflow.container</dd>
            <dt class="col-sm-3">Container engine</dt>
            <dd class="col-sm-9">$workflow.containerEngine</dd>
            <dt class="col-sm-3">Command line</dt>
            <dd class="col-sm-9">$workflow.commandLine</dd>
        </dl>

        <h1 class="border-bottom border-dark" id="software">Software</h1>
        <dl class="row small">
            <dt class="col-sm-3">smoove</dt>
            <dd class="col-sm-9">https://github.com/brentp/smoove</dd>
            <dt class="col-sm-3">indexcov</dt>
            <dd class="col-sm-9">https://github.com/brentp/goleft/tree/master/indexcov</dd>
            <dt class="col-sm-3">bpbio</dt>
            <dd class="col-sm-9">https://github.com/brentp/bpbio</dd>
            <dt class="col-sm-3">bcftools</dt>
            <dd class="col-sm-9">https://samtools.github.io/bcftools/</dd>
            <dt class="col-sm-3">samtools</dt>
            <dd class="col-sm-9">https://github.com/samtools/samtools</dd>
            <dt class="col-sm-3">nextflow</dt>
            <dd class="col-sm-9">https://www.nextflow.io/</dd>
        </dl>
    </div>
    <footer class="footer">
        <div class="container">
            <p class="text-muted text-right">Generated by <a href="https://github.com/brwnj/smoove-nf">smoove-nf</a>.</p>
        </div>
    </footer>
</body>
</html>
"""

pca_div = """
        <h3 id="scp_principal_components">Principle components</h3>
        <p>PCA is used to summarize samples by their non-sex chromosome bins. Groups of samples may indicate major batch effects in the underlying data. [<a href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-pca.md">docs</a>]</p>
        <div class="row">
            <div id="pca1_2" class="col-6"></div>
            <div id="pca1_3" class="col-6"></div>
            <script>
            layout.title = "PCA: 1 vs 2"
            layout.xaxis.title = "PC1"
            layout.yaxis.title = "PC2"
            pca_x_series = [PCA_x1]
            data = [{
                x: pca_x_series,
                y: [PCA_y1],
                mode: 'markers',
                type: 'scatter',
                name: 'Bins',
                text: samples,
                hoverinfo: 'text',
                marker: {
                    size: 12,
                    color: 'rgba(31,120,180,0.5)',
                    line: {
                        color: 'rgb(40,40,40)',
                        width: 1,
                    }
                },
            }]
            var plot2 = Plotly.newPlot('pca1_2', data, layout)

            layout.title = "PCA: 1 vs 3"
            layout.xaxis.title = "PC1"
            layout.yaxis.title = "PC3"
            data = [{
                x: pca_x_series,
                y: [PCA_y2],
                mode: 'markers',
                type: 'scatter',
                name: 'Bins',
                text: samples,
                hoverinfo: 'text',
                marker: {
                    size: 12,
                    color: 'rgba(31,120,180,0.5)',
                    line: {
                        color: 'rgb(40,40,40)',
                        width: 1,
                    }
                },
            }]
            var plot3 = Plotly.newPlot('pca1_3', data, layout)
            </script>
        </div>
"""

sequence_count_files = "$sequence_count".split(" ")
variant_count_files = "$variant_count".split(" ")
square_vcf_file = "$vcf"
ped_file = "$pedfile"
indexcov_roc_file = "$rocfile"
svvcf_html_file = "$variant_html"
bed_file = "$bedfile"
sex_chroms = "$sexchroms".split(",")

# building the sample summary table
## parse counts
sample_counts = defaultdict(dict)
for count_file in sequence_count_files:
    sample = os.path.basename(count_file).partition("-smoove-call")[0]
    with open(count_file) as fh:
        # [smoove]: ([E]lumpy-filter) 2019/01/16 21:45:01 [lumpy_filter] extracted splits and discordants from 701835557 total aligned reads
        for line in fh:
            if "total aligned reads" in line:
                count = int(line.partition("from ")[-1].partition(" total")[0])
                sample_counts[sample]["mapped"] = count
                sample_counts[sample]["variants"] = -1
                # initialize the status messages
                sample_counts[sample]["called"] = "fail"
                sample_counts[sample]["genotyped"] = "fail"
                break
## parse called
for count_file in variant_count_files:
    sample = os.path.basename(count_file).partition("-stats")[0]
    with open(count_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if "number of records" in line:
                count = int(line.split("\\t")[-1])
                sample_counts[sample]["variants"] = count
                sample_counts[sample]["called"] = "success"
                break
## parse genotyped from annotated vcf
filtering_counts = defaultdict(list)
with gzip.open(square_vcf_file, "rt") as fh:
    for line in fh:
        if not line.startswith("##"):
            break
        if line.startswith("##SAMPLE"):
            sample = line.strip().partition("ID=")[-1].strip(">")
            sample_counts[sample]["genotyped"] = "success"
        # building the read filtering plots
        elif line.startswith("##smoove_count_stats"):
            stats = line.strip().partition("=")[-1]
            sample = stats.partition(":")[0]
            split_before, disc_before, split_after, disc_after = stats.partition(":")[-1].split(",")
            filtering_counts["samples"].append(sample)
            filtering_counts["split_before"].append(split_before)
            filtering_counts["disc_before"].append(disc_before)
            filtering_counts["split_after"].append(split_after)
            filtering_counts["disc_after"].append(disc_after)
sample_list = sorted(sample_counts.keys())
sample_summary_str = ""
for sample in sample_list:
    counts = sample_counts[sample]
    sample_summary_str += "['{sample}', {mapped}, {variants}, {called}, {genotyped}],".format(sample=sample, mapped=counts["mapped"], variants=counts["variants"], called=counts["called"], genotyped=counts["genotyped"])
html = html.replace("SAMPLE_SUMMARY", sample_summary_str)
html = html.replace("FILTERED_SAMPLES", ",".join(["'{sample}'".format(sample=i) for i in filtering_counts["samples"]]))
html = html.replace("SPLIT_BEFORE", ",".join(filtering_counts["split_before"]))
html = html.replace("SPLIT_AFTER", ",".join(filtering_counts["split_after"]))
html = html.replace("DISCORDANT_BEFORE", ",".join(filtering_counts["disc_before"]))
html = html.replace("DISCORDANT_AFTER", ",".join(filtering_counts["disc_after"]))

# building sequence summary plots
ped_data = defaultdict(list)
pca = True
with open(ped_file) as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    for row in reader:
        # inferred sex
        if row["sex"] == "1":
            ped_data["inferred_samples1"].append(row["sample_id"])
            ped_data["inferred_x1"].append(row["CN%s" % sex_chroms[0]])
            try:
                ped_data["inferred_y1"].append(row["CN%s" % sex_chroms[1]])
            except IndexError:
                ped_data["inferred_y1"].append(0)
        else:
            ped_data["inferred_samples2"].append(row["sample_id"])
            ped_data["inferred_x2"].append(row["CN%s" % sex_chroms[0]])
            try:
                ped_data["inferred_y2"].append(row["CN%s" % sex_chroms[1]])
            except IndexError:
                ped_data["inferred_y1"].append(0)
        # bin plot
        total = float(row["bins.in"]) + float(row["bins.out"])
        ped_data["samples"].append(row["sample_id"])
        ped_data["bin_x"].append("%f" % (float(row["bins.lo"]) / total))
        ped_data["bin_y"].append("%f" % (float(row["bins.out"]) / total))
        # PCAs
        try:
            ped_data["pca_1"].append(row["PC1"])
            ped_data["pca_2"].append(row["PC2"])
            ped_data["pca_3"].append(row["PC3"])
        except KeyError:
            pca = False
            pass
html = html.replace("COVERAGE_SAMPLE_LIST", ",".join(["'{sample}'".format(sample=i) for i in ped_data["samples"]]))
html = html.replace("INFERRED_X1", ",".join(ped_data["inferred_x1"]))
html = html.replace("INFERRED_Y1", ",".join(ped_data["inferred_y1"]))
html = html.replace("INFERRED_SAMPLES1", ",".join(["'{sample}'".format(sample=i) for i in ped_data["inferred_samples1"]]))
html = html.replace("INFERRED_X2", ",".join(ped_data["inferred_x2"]))
html = html.replace("INFERRED_Y2", ",".join(ped_data["inferred_y2"]))
html = html.replace("INFERRED_SAMPLES2", ",".join(["'{sample}'".format(sample=i) for i in ped_data["inferred_samples2"]]))
html = html.replace("BIN_X", ",".join(ped_data["bin_x"]))
html = html.replace("BIN_Y", ",".join(ped_data["bin_y"]))

if pca:
    pca_div = pca_div.replace("PCA_x1", ",".join(ped_data["pca_1"]))
    pca_div = pca_div.replace("PCA_y1", ",".join(ped_data["pca_2"]))
    pca_div = pca_div.replace("PCA_y2", ",".join(ped_data["pca_3"]))
else:
    pca_div = ""
html = html.replace("[PCA_DIV]", pca_div)

# fixing the file link to indexcov results
index_cov_output = "$outdir/indexcov/index.html".replace("s3://", "https://s3.amazonaws.com/")
html = html.replace("INDEXCOV_RESULT", '<a href="{path}">{path}</a>'.format(path=index_cov_output))

# build the chromosome coverage plots for 'by percentage'
allowable = ["{}".format(i) for i in list(range(1,23))] + ["X", "Y"]
with open(indexcov_roc_file) as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    data = defaultdict(lambda: defaultdict(list))
    for row in reader:
        # header row is repeated at chromosome breaks
        if row["#chrom"] == "#chrom":
            continue
        chrom = row["#chrom"].strip("chr")
        if chrom not in allowable:
            continue
        data[chrom]["x"].append(row["cov"])
        for sample in sample_list:
            # KeyError here would be interesting
            data[chrom][sample].append(row[sample])
buttons = ""
cov_data = "{"
for chrom in allowable:
    if chrom in data:
        if not buttons:
            buttons += '<label class="btn btn-secondary active"><input type="radio" name="options" data-name="{chrom}" autocomplete="off" checked> {chrom} </label>'.format(chrom=chrom)
        else:
            buttons += '<label class="btn btn-secondary"><input type="radio" name="options" data-name="{chrom}" autocomplete="off"> {chrom} </label>'.format(chrom=chrom)
        trace_data = "["
        for sample in sample_list:
            trace_data += """{{x:[{x}],y:[{y}],hoverinfo:'text',mode:'lines',text:'{sample}',marker:{{ color:cov_color }}}},""".format(x=",".join(data[chrom]["x"]), y=",".join(data[chrom][sample]), sample=sample)
        trace_data += "]"
    cov_data += "'{chrom}': {trace},".format(chrom=chrom, trace=trace_data)
cov_data += "}"
html = html.replace("[CHROM_BUTTONS]", buttons)
html = html.replace("[CHROM_DATA]", cov_data)

# build the chromosome coverage plots for 'by position'
with gzip.open(bed_file, "rt") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    data = defaultdict(lambda: defaultdict(list))
    for row in reader:
        chrom = row["#chrom"].strip("chr")
        if chrom not in allowable:
            continue
        data[chrom]["x"].append(row["start"])
        for sample in sample_list:
            data[chrom][sample].append(row[sample])
for chrom in allowable:
    with open("%s.json" % chrom, "w") as fh:
        traces = list()
        for i, sample in enumerate(sample_list):
            trace = dict(
                x=data[chrom]["x"],
                y=data[chrom][sample],
                hoverinfo="text",
                type="scattergl",
                mode="lines",
                text=sample,
                marker={"color":"rgba(108,117,125,0.3)"}
            )
            if i < 5:
                # ensures lines are in the range selection preview
                trace["type"] = "scatter"
            traces.append(trace)
        print(json.dumps(traces), file=fh)

# build the variant summary plots
var_samples = []
plot_position = 0
histograms = {1: "deletions", 3: "duplications", 5: "inversions", 7: "bnds"}
plot_data = defaultdict(list)
with open(svvcf_html_file) as fh:
    for line in fh:
        line = line.strip()
        if line.startswith("var pdata"):
            false = False
            true = True
            data = eval(line.partition("=")[-1].strip(" ;"))
            # bar
            if len(data) > 1:
                for trace in data:
                    if not var_samples:
                        var_samples = [i.replace("sample:", "") for i in trace["x"]]
                    plot_data[trace["name"]] = trace["y"]
            # hist
            else:
                plot_data[histograms[plot_position]] = data[0]["x"]
            plot_position += 1
html = html.replace("VARIANT_SAMPLE_LIST", ",".join(["'{}'".format(i) for i in var_samples]))
html = html.replace("SMALL_DELETIONS", ",".join(map(str, plot_data["small deletions"])))
html = html.replace("LARGE_DELETIONS", ",".join(map(str, plot_data["large deletions"])))
html = html.replace("DELETIONS_HIST", ",".join(map(str, plot_data["deletions"])))
html = html.replace("SMALL_DUPLICATIONS", ",".join(map(str, plot_data["small duplications"])))
html = html.replace("LARGE_DUPLICATIONS", ",".join(map(str, plot_data["large duplications"])))
html = html.replace("DUPLICATIONS_HIST", ",".join(map(str, plot_data["duplications"])))
html = html.replace("SMALL_INVERSIONS", ",".join(map(str, plot_data["small inversion"])))
html = html.replace("LARGE_INVERSIONS", ",".join(map(str, plot_data["large inversion"])))
html = html.replace("INVERSIONS_HIST", ",".join(map(str, plot_data["inversions"])))
html = html.replace("SMALL_BNDS", ",".join(map(str, plot_data["small BNDs"])))
html = html.replace("LARGE_BNDS", ",".join(map(str, plot_data["large BNDs"])))
try:
    html = html.replace("INTERCHROMOSOMAL_BNDS", ",".join(map(str, plot_data["interchromosomal BNDs"])))
except KeyError:
    pass
html = html.replace("BNDS_HIST", ",".join(map(str, plot_data["bnds"])))

with open("smoove-nf.html", "w") as fh:
    print(html, file=fh)
