process DADA2_TAXONOMY {
    tag "${fasta},${database}"
    label 'process_high'
    publishDir "DADA2_TAXONOMY", mode: 'copy'

    conda "/g/korbel2/weber/miniconda3/envs/idtaxa_env"
    // conda "bioconda::bioconductor-dada2=1.22.0 conda-forge::r-digest=0.6.30"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
    //     'biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    path(fasta)
    path(database)
    val(outfile)
    val(taxlevels_input)

    output:
    path(outfile), emit: tsv
    path( "ASV_tax.rds" ), emit: rds
    path "versions.yml"  , emit: versions
    path "*.args.txt"    , emit: args

    when:
    task.ext.when == null || task.ext.when

    // script:
    // def args = task.ext.args ?: ''
    // def taxlevels = taxlevels_input ?
    //     'c("' + taxlevels_input.split(",").join('","') + '")' :
    //     'c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")'
    // def seed = task.ext.seed ?: '100'
    """
    Rscript /g/korbel2/weber/workspace/ampliseq/bin/idtaxa.R $fasta $task.cpus $outfile
    """
}
