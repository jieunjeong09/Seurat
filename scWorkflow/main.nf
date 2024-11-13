#!/usr/bin/env nextflow

params {
    input_dir = "./cellranger_outputs"  // CellRanger output
    output_dir = "./results"
}

process preprocess {
    container 'scworkflow:r-latest'  // Docker
    input:
    path input_dir from params.input_dir
    output:
    path "preprocessed_data.rds"
    script:
    """
    Rscript bin/preprocess.R ${input_dir} preprocessed_data.rds
    """
}

process analysis {
    container 'scworkflow:r-latest'
    input:
    path data from preprocess.out
    output:
    path params.output_dir
    script:
    """
    Rscript bin/analysis.R ${data} ${params.output_dir}
    """
}

process generate_rmd {
    container 'scworkflow:r-latest'
    input:
    path params.output_dir from analysis.out
    output:
    path "${params.output.dir}/report.html"
    script:
    """
    R -e rmarkdown::render('Rmarkdown/report_template.Rmd',
       output_file="${params.output_dir}/report.html")
    """
}

workflow {
    preprocess()
    analysis()
    generate_rmd()
}
