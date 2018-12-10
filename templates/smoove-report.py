#!/usr/bin/env python
from __future__ import print_function

import csv
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
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
    <nav class="navbar navbar-expand-md navbar-dark bg-dark">
        <a class="navbar-brand" href="#">smoove-nf report</a>
        <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarCollapse" aria-controls="navbarCollapse" aria-expanded="false" aria-label="Toggle navigation">
            <span class="navbar-toggler-icon"></span>
        </button>
        <div class="collapse navbar-collapse" id="navbarCollapse">
            <ul class="navbar-nav mr-auto">
                <li class="nav-item"><a class="nav-link" href="#">Configuration</a></li>
                <li class="nav-item"><a class="nav-link" href="#summary">Summary</a></li>
                <li class="nav-item"><a class="nav-link" href="#filtering">Filtering</a></li>
                <li class="nav-item"><a class="nav-link" href="#software">Software</a></li>
            </ul>
            <span class="navbar-text">
                [<a href="https://github.com/brwnj/smoove-nf" target="_blank">code</a>, <a href="https://github.com/brwnj/smoove-nf/issues" target="_blank">issues</a>]
            </span>
        </div>
    </nav>

    <div class="jumbotron mb-0 pt-3 pb-3">
        <div class="container">
            <h1 id="configuration" style="padding-top: 20px;">Configuration</h1>
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
        </div>
    </div>

    <div class="container pb-20">
        <h1 id="summary" style="padding-top: 20px;">Sample Summary</h1>
        <table id="sample_table" class="table table-striped" width="100%"></table>
        <script>
        var dataSet = [SAMPLE_SUMMARY];

        \$(document).ready(function() {
            \$('#sample_table').DataTable( {
                data: dataSet,
                columns: [
                    { title: "Sample" },
                    { title: "Reads" },
                    { title: "Called Variants" },
                    { title: "Called" },
                    { title: "Genotyped" },
                ]
            } );
        } );
        </script>

        <h1 id="sequence_summary" style="padding-top: 20px;">Sequence Summary</h1>
        <p><a href="https://github.com/brentp/goleft/tree/master/indexcov">indexcov</a> is used to do a basic QC of the data, screen for sex chromosome anomalies, analyze coverage across samples, and highlight any obvious deletions or duplications.</p>
        <p><a href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-sex.md">Inferred sex</a> and <a href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-bin.md">bin depth</a> analysis plots:</p>
        <div class="row">
            <div id="inferred_sex" class="col-6"></div>
            <div id="bin_counts" class="col-6"></div>
            <script>
            var samples = [SAMPLE_LIST];
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

        <p><a href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-pca.md">Principal component analysis</a> to highlight potential major batch effects:</p>
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
            var plot1 = Plotly.newPlot('pca1_2', data, layout)

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
            var plot1 = Plotly.newPlot('pca1_3', data, layout)
            </script>
        </div>
        <h3>indexcov Results</h3>
        <p>INDEXCOV_RESULT</p>

        <h1 id="filtering" style="padding-top: 20px;">Read Filtering</h1>
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

        <h1 id="software" style="padding-top: 20px;">Software</h1>
        <dl class="row small">
            <dt class="col-sm-3">smoove</dt>
            <dd class="col-sm-9">https://github.com/brentp/smoove</dd>
            <dt class="col-sm-3">indexcov</dt>
            <dd class="col-sm-9">https://github.com/brentp/goleft/tree/master/indexcov</dd>
            <dt class="col-sm-3">bcftools</dt>
            <dd class="col-sm-9">https://samtools.github.io/bcftools/</dd>
            <dt class="col-sm-3">samtools</dt>
            <dd class="col-sm-9">https://github.com/samtools/samtools</dd>
            <dt class="col-sm-3">nextflow</dt>
            <dd class="col-sm-9">https://www.nextflow.io/</dd>
        </dl>
    </div>
    <footer class="footer" style="padding-top: 20px;padding-bottom: 20px;">
        <div class="container">
            <p class="text-muted text-right">Generated by <a href="https://github.com/brwnj/smoove-nf">smoove-nf</a>.</p>
        </div>
    </footer>
</body>
</html>
"""

# building the sample summary table
## parse counts
sample_counts = defaultdict(dict)
for count_file in "$sequence_count".split(" "):
    sample = os.path.basename(count_file).partition("-flagstat")[0]
    with open(count_file) as fh:
        for line in fh:
            if "mapped" in line:
                count = int(line.partition(" ")[0])
                sample_counts[sample]["mapped"] = count
                sample_counts[sample]["variants"] = -1
                # initialize the status messages
                sample_counts[sample]["called"] = '<span class="badge badge-danger">Fail</span>'
                sample_counts[sample]["genotyped"] = '<span class="badge badge-danger">Fail</span>'
                break
# could use later for validation, ordering in js
sample_list = sorted(sample_counts.keys())
## parse called
for count_file in "$variant_count".split(" "):
    sample = os.path.basename(count_file).partition("-stats")[0]
    with open(count_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if "number of records" in line:
                count = int(line.split("\t")[-1])
                sample_counts[sample]["variants"] = count
                sample_counts[sample]["called"] = '<span class="badge badge-success">Success</span>'
                break
## parse genotyped from annotated vcf
with open("$vcf") as fh:
    for line in fh:
        if not line.startswith("##"):
            break
        if line.startswith("##SAMPLE"):
            ##SAMPLE=<ID=1014>
            toks = line.split("=")
            sample = toks[-1].strip(">")
            sample_counts[sample]["genotyped"] = '<span class="badge badge-success">Success</span>'

sample_summary_str = ""
for sample in sample_list:
    counts = sample_counts[sample]
    sample_summary_str += "['{sample}', {mapped}, {variants}, '{called}', '{genotyped}'],".format(sample=sample, mapped=counts["mapped"], variants=counts["variants"], called=counts["called"], genotyped=counts["genotyped"])
html = html.replace("SAMPLE_SUMMARY", sample_summary_str)

# building sequence summary plots
ped_data = defaultdict(list)
with open("$ped") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        # inferred sex
        if row["sex"] == "1":
            ped_data["inferred_samples1"].append(row["sample_id"])
            ped_data["inferred_x1"].append(row["CNX"])
            ped_data["inferred_y1"].append(row["CNY"])
        else:
            ped_data["inferred_samples2"].append(row["sample_id"])
            ped_data["inferred_x2"].append(row["CNX"])
            ped_data["inferred_y2"].append(row["CNY"])
        # bin plot
        total = float(row["bins.in"]) + float(row["bins.out"])
        ped_data["samples"].append(row["sample_id"])
        ped_data["bin_x"].append(float(row["bins.lo"]) / total)
        ped_data["bin_y"].append(float(row["bins.out"]) / total)
        # PCAs
        ped_data["pca_1"].append(row["PC1"])
        ped_data["pca_2"].append(row["PC2"])
        ped_data["pca_3"].append(row["PC3"])
html = html.replace("SAMPLE_LIST", ",".join(["'{sample}'".format(sample=i) for i in ped_data["samples"]]))
html = html.replace("INFERRED_X1", ",".join(ped_data["inferred_x1"]))
html = html.replace("INFERRED_Y1", ",".join(ped_data["inferred_y1"]))
html = html.replace("INFERRED_SAMPLES1", ",".join(["'{sample}'".format(sample=i) for i in ped_data["inferred_samples1"]]))
html = html.replace("INFERRED_X2", ",".join(ped_data["inferred_x2"]))
html = html.replace("INFERRED_Y2", ",".join(ped_data["inferred_y2"]))
html = html.replace("INFERRED_SAMPLES2", ",".join(["'{sample}'".format(sample=i) for i in ped_data["inferred_samples2"]]))
html = html.replace("BIN_X", ",".join(ped_data["bin_x"]))
html = html.replace("BIN_Y", ",".join(ped_data["bin_y"]))
html = html.replace("PCA_x1", ",".join(ped_data["pca_1"]))
html = html.replace("PCA_y1", ",".join(ped_data["pca_2"]))
html = html.replace("PCA_y2", ",".join(ped_data["pca_3"]))

# fixing the file link to indexcov results
index_cov_output = "$outdir/indexcov/index.html".replace("s3://", "https://s3.amazonaws.com/")
html = html.replace("INDEXCOV_RESULT", '<a href="{path}">{path}</a>'.format(path=index_cov_output))

# building the read filtering plots
filtering_counts = defaultdict(list)
with open("$variant_count") as fh:
    for line in fh:
        if not line.startswith("##"):
            break
        if line.startswith("##smoove_count_stats"):
            stats = line.strip().partition("=")[-1]
            sample = stats.partition(":")[0]
            split_before, disc_before, split_after, disc_after = stats.partition(":")[-1].split(",")
            filtering_counts["samples"].append(sample)
            filtering_counts["split_before"].append(split_before)
            filtering_counts["disc_before"].append(disc_before)
            filtering_counts["split_after"].append(split_after)
            filtering_counts["disc_after"].append(disc_after)
html = html.replace("FILTERED_SAMPLES", ",".join(["'{sample}'".format(sample=i) for i in filtering_counts["samples"]]))
html = html.replace("SPLIT_BEFORE", ",".join(filtering_counts["split_before"]))
html = html.replace("SPLIT_AFTER", ",".join(filtering_counts["split_after"]))
html = html.replace("DISCORDANT_BEFORE", ",".join(filtering_counts["disc_before"]))
html = html.replace("DISCORDANT_AFTER", ",".join(filtering_counts["disc_after"]))

with open("smoove-nf.html", "w") as fh:
    print(html, file=fh)
