process DADA2_TAXONOMY {
    tag "${fasta},${database}"
    label 'process_high'
    publishDir "DADA2_TAXONOMY", mode: 'copy'

    // conda "bioconda::bioconductor-decipher=2.28.0 conda-forge::r-data.table=1.14.8 bioconda::bioconductor-biostrings=2.68.1"
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
    #!/usr/bin/env Rscript


    library(data.table)
    library(DECIPHER)
    library(Biostrings)

    # args <- commandArgs(TRUE)

    # input <- args[1]
    # input <- $fasta
    # threads <- args[2] |> as.numeric()
    threads <- $task.cpus |> as.numeric()
    trainingSet.path <- "/g/korbel/weber/TMP/VINCENT_GROUP/pr2_version_4.14.0_SSU.decipher.trained.rds"
    # output <- args[3]
    #output <- $outfile
    # output_versions <- args[4]
    # output_args <- args[5]
    # output_rds <- args[6]

    trainingSet <- readRDS(trainingSet.path)

    otus <- Biostrings::readDNAStringSet(\"$fasta\")
    #otus <- Biostrings::readDNAStringSet(input)

    otus_assign <- DECIPHER::IdTaxa(otus,
            trainingSet,
            strand = "top",
            threshold = 50,
            processors = threads
    )

    taxo <- vapply(
            otus_assign,
            function(X) paste(X\$taxon, collapse = ";"),
            character(1)
    )

    confidence <- vapply(
            otus_assign,
            function(X) paste(round(X\$confidence, digits = 1), collapse = ";"),
            character(1)
    )

    data.table(names(taxo), taxo, confidence) |> fwrite(
            file = \"$outfile\",
            sep = "\t",
            quote = FALSE,
            col.names = FALSE
    )

    saveRDS(data.table(names(taxo), taxo, confidence), "ASV_tax.rds")
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: 1.28.0") ), "versions.yml")
    # writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    writeLines(c("TEST_assignTaxonomy"), "assignTaxonomy.args.txt")

    """
}
