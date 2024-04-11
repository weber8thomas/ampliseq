#!/usr/bin/env Rscript

library(data.table)
library(DECIPHER)
library(Biostrings)

args <- commandArgs(TRUE)

input <- args[1]
threads <- args[2] |> as.numeric()
trainingSet.path <- "/g/korbel/weber/TMP/VINCENT_GROUP/pr2_version_4.14.0_SSU.decipher.trained.rds"
output <- args[3]
# output_versions <- args[4]
# output_args <- args[5]
# output_rds <- args[6]

trainingSet <- readRDS(trainingSet.path)

otus <- Biostrings::readDNAStringSet(input)

otus_assign <- DECIPHER::IdTaxa(otus,
        trainingSet,
        strand = "top",
        threshold = 50,
        processors = threads
)

taxo <- vapply(
        otus_assign,
        function(X) paste(X$taxon, collapse = ";"),
        character(1)
)

confidence <- vapply(
        otus_assign,
        function(X) paste(round(X$confidence, digits = 1), collapse = ";"),
        character(1)
)

data.table(names(taxo), taxo, confidence) |> fwrite(
        file = output,
        sep = "\t",
        quote = FALSE,
        col.names = FALSE
)

saveRDS(data.table(names(taxo), taxo, confidence), "ASV_tax.rds")
writeLines(c("\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
writeLines(c("TEST_assignTaxonomy"), "assignTaxonomy.args.txt")
# writeLines(c("TEST_assignTaxonomy"), "versions.yml")
write.table("assignTaxonomy\t$args\ntaxlevels\t$taxlevels\nseed\t$seed", file = "assignTaxonomy.args.txt")
