#!/usr/bin/env Rscript
###############################################################################
# GTDB Dataset creation tool
#   * Selects taxonomic subset from GTDB metadata
#   * Applies quality filters
#   * Samples genomes with different strategies
#   * Downloads FASTA files via NCBI datasets CLI
###############################################################################

suppressMessages(
        for (package in c("data.table", "R.utils", "cluster")) {
                if (!require(package, character.only = T, quietly = T)) {
                        install.packages(package)
                        library(package, character.only = T)
                }
        }
)

###############################################################################
# 1. Command-line interface
###############################################################################

### Parameters
args <- commandArgs(trailingOnly = TRUE)

### Defaults
name_selection <- "Pseudomonas aeruginosa"
taxrank <- "species"
representative <- "f"
dryrun <- FALSE
output_dir <- "output"
domain <- "bacteria"
database_file <- NULL
## Sampling specific
n_samples <- 20
quality <- "good"
sampling_method <- "custom"
sampling_proportion <- c("A:0.5, O:0.25, R:0.25")
## Genome Size bands
min_size <- 0
max_size <- Inf

## Logic for default sampling_proportion
average_samples <- as.numeric(gsub("A:(.*?),.*", "\\1", sampling_proportion)) * n_samples
outliers_samples <- round(as.numeric(gsub(".*O:([0-9.]+)[, ].*", "\\1", sampling_proportion)) * n_samples)
random_samples <- round(as.numeric(gsub(".*R:([0-9.]+).*", "\\1", sampling_proportion)) * n_samples)

# Parse named arguments
i <- 1
while (i <= length(args)) {
        if (args[i] %in% c("--name_selection", "-n")) {
                name_selection <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--taxrank", "-t")) {
                taxrank <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--n_samples", "-ns")) {
                n_samples <- args[i + 1]
                n_samples <- as.numeric(n_samples)
                i <- i + 2
        } else if (args[i] %in% c("--sampling", "-ss")) {
                sampling_proportion <- args[i + 1]
                average_samples <- as.numeric(gsub("A:(.*?),.*", "\\1", sampling_proportion)) * n_samples
                outliers_samples <- round(as.numeric(gsub(".*O:([0-9.]+)[, ].*", "\\1", sampling_proportion)) * n_samples)
                random_samples <- round(as.numeric(gsub(".*R:([0-9.]+).*", "\\1", sampling_proportion)) * n_samples)
                i <- i + 2
        } else if (args[i] %in% c("--quality", "-q")) {
                quality <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--min_size")) {
                min_size <- as.numeric(args[i + 1])
                i <- i + 2
        } else if (args[i] %in% c("--max_size")) {
                max_size <- as.numeric(args[i + 1])
                i <- i + 2
        } else if (args[i] %in% c("--representative_check", "-r")) {
                representative_check <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--dryrun", "-d")) {
                dryrun <- TRUE
                i <- i + 1
        } else if (args[i] %in% c("--output_dir", "-o")) {
                output_dir <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--domain", "-m")) {
                domain <- tolower(args[i + 1])
                if (!domain %in% c("bacteria", "archaea")) {
                        stop("Domain must be either 'bacteria' or 'archaea'")
                }
                i <- i + 2
        } else if (args[i] %in% c("--database", "-db")) {
                database_file <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--sampling_method")) {
                sampling_method <- tolower(args[i + 1])
                i <- i + 2
        } else if (args[i] %in% c("--help", "-h")) {
                cat("Usage: gtdb_downloader.R [options]\n\n")
                cat("Required:\n")
                cat("  -n, --name_selection <str>  Target taxon name\n")
                cat("  -ns, --n_samples <int>        Number of genomes to sample\n")
                cat("       --sampling_method <custom|kmedoids|maxmin|stratified|hybrid>\n")
                cat("Optional:\n")
                cat("  -ss, --sampling <str>       Proportion schema (A:avg,O:out,R:rand)\n\n")
                cat("  -t,  --taxrank <rank>       [species|genus|family]  (default: species)\n")
                cat("  -q,  --quality <level>     [good|bad|mixed]        (default: good)\n")
                cat("       --min_size <num>                               (default: 0)\n")
                cat("       --max_size <num>                               (default: Inf)\n")
                cat("  -r,  --representative <t|f>  Only representatives  (default: f)\n")
                cat("  -m,  --domain <bacteria|archaea>                 (default: bacteria)\n")
                cat("  -o,  --output_dir <dir>                           (default: output)\n")
                cat("  -db, --database <file>  Custom GTDB TSV           (optional)\n")
                cat("       --sampling_method <custom|kmedoids|maxmin|stratified|hybrid>\n")
                cat("  -d,  --dryrun                                   (default: FALSE)\n")
                cat("  -h,  --help                                     Show this text\n")
                quit(status = 0)
        } else {
                # Fallback to positional arguments
                if (is.null(name_selection)) {
                        name_selection <- args[i]
                } else if (i == 2) {
                        taxrank <- args[i]
                } else if (i == 3) {
                        representative_check <- args[i]
                } else if (i == 4) dryrun <- as.logical(args[i])
                i <- i + 1
        }
}

# Check if name_selection was givven
if (is.null(name_selection)) stop("--name_selection is required")
if (!taxrank %in% c("species", "genus", "family")) {
        stop("taxrank must be species, genus, or family")
}
if (!quality %in% c("good", "bad", "mixed")) {
        stop("quality must be good, bad, or mixed")
}
if (!representative %in% c("t", "f")) {
        stop("representative must be t or f")
}
if (!domain %in% c("bacteria", "archaea")) {
        stop("domain must be bacteria or archaea")
}

###############################################################################
# 4. Console header
###############################################################################


cat("\n\033[1mGTDB Dataset Builder Parameters\033[0m\n")
cat(rep("-", 60), "\n", sep = "")

# Core parameters
cat("\033[1;34mCore Settings:\033[0m\n")
cat(sprintf("  %-25s: %s", "Taxonomic name", name_selection), "\n")
cat(sprintf("  %-25s: %s", "Taxonomic rank", taxrank), "\n")
cat(sprintf("  %-25s: %s", "Domain", domain), "\n")
cat(sprintf("  %-25s: %s", "Output directory", output_dir), "\n")
cat(sprintf(
        "  %-25s: %s", "Representative genomes",
        ifelse(representative == "t", "Yes", "No")
), "\n")

# Sampling parameters
cat("\n\033[1;34mSampling Strategy:\033[0m\n")
cat(sprintf("  %-25s: %d", "Total samples", n_samples), "\n")
cat(sprintf("  %-25s: %s", "Quality level", quality), "\n")
cat(sprintf("  %-25s: %s", "Sampling schema", sampling_proportion), "\n")
cat(sprintf("  %-25s: %s", "Sampling method", sampling_method), "\n")
cat(sprintf(
        "  %-25s: %s", "Average genomes",
        average_samples
), "\n")
cat(sprintf(
        "  %-25s: %s", "Outlier genomes",
        outliers_samples
), "\n")
cat(sprintf(
        "  %-25s: %s", "Random genomes",
        random_samples
), "\n")

# Operational parameters
cat("\n\033[1;34mOperational Settings:\033[0m\n")
cat(sprintf(
        "  %-25s: %s", "Dry run mode",
        ifelse(dryrun, "Yes", "No")
), "\n")
cat(sprintf(
        "  %-25s: %s", "Database",
        ifelse(is.null(database_file), "Default", database_file)
), "\n")

cat(rep("-", 60), "\n\n", sep = "")

############################
## NOTE: 2. Tool download ##
############################

### Database
if (!is.null(database_file)) {
        # Custom database
        if (!file.exists(database_file)) {
                stop(paste("Specified database file not found:", database_file))
        }
        cat("Loading custom database from:", database_file, "\n")
        data <- fread(database_file)
} else {
        # NOTE: DATABASE UPDATE MODIFY JUST THIS LINES
        # Default GTDB database
        if (domain == "bacteria") {
                default_files <- c("bac120_metadata_r226.tsv.gz", "bac120_metadata_r226.tsv")
                download_url <- "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/bac120_metadata_r226.tsv.gz"
        } else if (domain == "archaea") {
                default_files <- c("ar53_metadata_r226.tsv.gz", "ar53_metadata_r226.tsv")
                download_url <- "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/ar53_metadata_r226.tsv.gz"
        }
        ## File check
        file_detection <- list.files()
        file_check <- default_files %in% file_detection
        if (!any(file_check)) {
                cat("Downloading", domain, "metadata from GTDB...\n")
                data <- fread(download_url, header = TRUE)
                local_file <- default_files[2] # Save as uncompressed TSV
                fwrite(data, local_file, sep = "\t")
                cat("Saved database to:", local_file, "\n")
        } else {
                existing_file <- default_files[which(file_check)[1]]
                cat("GTDB", domain, "database found, loading:", existing_file, "\n")
                data <- fread(existing_file)
        }
}

## NCBI dataset
if (!file.exists("datasets")) {
        cat("Downloading NCBI datasets CLI tool...\n")
        download.file(
                "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets",
                "datasets",
                quiet = TRUE
        )
        system("chmod +x datasets", ignore.stdout = TRUE, ignore.stderr = TRUE)
        cat("NCBI datasets CLI tool installed and made executable\n")
} else {
        cat("NCBI datasets CLI tool found\n")
}

###############################################################################
# 6. Metadata processing
###############################################################################

tax_split <- function(dt) {
        ## tax split
        tax_levels <- c("phylium", "class", "order", "family", "genus", "species")
        dt[, (tax_levels) := tstrsplit(gtdb_taxonomy, ";")[2:7]]
        ## Pref remove
        for (col in tax_levels) {
                dt[, (col) := sub("^[a-z]__", "", get(col))]
        }
        return(dt)
}

processed_data <- tax_split(data)
###############################################################################
# 7. Subset selection
###############################################################################

selection <- processed_data[
        get(taxrank) == name_selection &
                gtdb_representative == representative &
                genome_size >= min_size &
                genome_size <= max_size,
]

if (nrow(selection) == 0) {
        print("No genomes match the specified criteria.")
        quit(status = 150)
}
if (n_samples > nrow(selection)) {
        warning("Requested samples > available genomes. Using all ", nrow(selection), " entries.")
        n_samples <- nrow(selection)
}


###############################################################################
# 8. Quality filtering
###############################################################################

# Define thresholds based on external 'quality' variable
if (quality == "good") {
        contamination_thresh <- as.numeric(quantile(selection$checkm2_contamination, 0.25)) # Lower = better
        completeness_thresh <- as.numeric(quantile(selection$checkm2_completeness, 0.75)) # Higher = better
        contig_count_thresh <- as.numeric(quantile(selection$contig_count, 0.25)) # Lower = better
} else if (quality == "bad") {
        contamination_thresh <- as.numeric(quantile(selection$checkm2_contamination, 0.75)) # Higher = worse
        completeness_thresh <- as.numeric(quantile(selection$checkm2_completeness, 0.25)) # Lower = worse
        contig_count_thresh <- as.numeric(quantile(selection$contig_count, 0.75)) # Higher = worse
} else if (quality == "mixed") {
        # Skip filtering entirely
        selection_quality <- selection
        message("Quality 'mixed': No thresholds applied. Using all ", nrow(selection), " samples.")
} else {
        stop("Invalid 'quality' value. Must be 'good', 'bad', or 'mixed'.")
}


# Apply thresholds (skip if quality = "mixed")
if (quality == "good") {
        selection_quality <- selection[
                checkm2_contamination <= contamination_thresh &
                        checkm2_completeness >= completeness_thresh &
                        contig_count <= contig_count_thresh,
        ]
        relax_increment <- 0.05
        relax_steps <- 0
        while (nrow(selection_quality) < n_samples & relax_steps < 4 & quality != "mixed") {
                relax_steps <- relax_steps + 1
                contamination_thresh <- contamination_thresh * (1 + relax_increment)
                completeness_thresh <- completeness_thresh * (1 - relax_increment)
                contig_count_thresh <- contig_count_thresh * (1 + relax_increment)
                selection_quality <- selection[
                        selection$checkm2_contamination <= contamination_thresh &
                                selection$checkm2_completeness >= completeness_thresh &
                                selection$contig_count <= contig_count_thresh,
                ]
                message(
                        "Relaxing thresholds (step ", relax_steps, "): ",
                        "Contamination ≤ ", round(contamination_thresh, 2), ", ",
                        "Completeness ≥ ", round(completeness_thresh, 2), ", ",
                        "Contigs ≤ ", round(contig_count_thresh, 2), ". ",
                        "Now ", nrow(selection_quality), " samples available."
                )
        }
} else if (quality == "bad") {
        selection_quality <- selection[
                checkm2_contamination >= contamination_thresh &
                        checkm2_completeness <= completeness_thresh
        ]
        relax_increment <- 0.05
        relax_steps <- 0
        while (nrow(selection_quality) < n_samples & relax_steps < 4 & quality != "mixed") {
                relax_steps <- relax_steps + 1
                contamination_thresh <- contamination_thresh * (1 + relax_increment)
                completeness_thresh <- completeness_thresh * (1 - relax_increment)

                selection_quality <- selection[
                        selection$checkm2_contamination <= contamination_thresh &
                                selection$checkm2_completeness >= completeness_thresh
                ]
                message(
                        "Relaxing thresholds (step ", relax_steps, "): ",
                        "Contamination ≤ ", round(contamination_thresh, 2), ", ",
                        "Completeness ≥ ", round(completeness_thresh, 2), ", ",
                        "Now ", nrow(selection_quality), " samples available."
                )
        }
}
# Final sampling
if (nrow(selection_quality) == 0) {
        stop("No samples met even relaxed quality thresholds. Adjust criteria or data.")
} else if (nrow(selection_quality) < n_samples) {
        warning(
                "Only ", nrow(selection_quality), " high-quality samples available. ",
                "Returning all (fewer than ", n_samples, ")."
        )
        selection_final <- selection_quality
} else {
        selection_final <- selection_quality
}

rnames <- selection_final[, c("ncbi_genbank_assembly_accession")][["ncbi_genbank_assembly_accession"]]
to_num_cols <- c("lsu_23s_contig_len", "lsu_23s_length", "lsu_5s_contig_len", "lsu_5s_length", "lsu_silva_23s_blast_align_len", "lsu_silva_23s_blast_bitscore", "lsu_silva_23s_blast_perc_identity", "ncbi_contig_count", "ncbi_contig_n50", "ncbi_ncrna_count", "ncbi_protein_count", "ncbi_rrna_count", "ncbi_scaffold_count", "ncbi_scaffold_l50", "ncbi_scaffold_n50", "ncbi_scaffold_n75", "ncbi_scaffold_n90", "ncbi_ssu_count", "ncbi_trna_count", "ssu_contig_len", "ssu_gg_blast_align_len", "ssu_gg_blast_bitscore", "ssu_gg_blast_evalue", "ssu_gg_blast_perc_identity", "ssu_length", "ssu_silva_blast_align_len", "ssu_silva_blast_bitscore", "ssu_silva_blast_perc_identity")

distance_variables <- c(
        # Category: Identification and Accession
        "accession",
        # Category: Taxonomy
        "gtdb_taxonomy",
        "ncbi_taxonomy",
        "phylium",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "ncbi_translation_table",
        ## Category: Basic Stats
        "gc_percentage",
        "genome_size",
        # Category: Genome Quality
        "checkm2_completeness",
        "checkm2_contamination",
        "checkm_strain_heterogeneity",
        # Category: Assembly Statistics
        "n50_contigs",
        "contig_count",
        "longest_contig",
        "ncbi_contig_n50",
        "total_gap_length",
        "ncbi_total_gap_length",
        "ncbi_contig_count",
        # Category: Gene and RNA Counts
        "ssu_count",
        "lsu_23s_count",
        "protein_count",
        "ncbi_ncrna_count",
        "ncbi_rrna_count",
        "ncbi_trna_count",
        "trna_aa_count",
        "trna_selenocysteine_count",
        # Category: Environmental and Geodata
        "ncbi_genome_category",
        "ncbi_bioproject",
        "ncbi_country",
        "ncbi_genome_representation",
        "ncbi_isolation_source",
        "ncbi_lat_lon"
)

# Convert all character columns to factors in-place
char_cols <- names(selection_final)[sapply(selection_final, is.character)]
selection_final[, (char_cols) := lapply(.SD, as.factor), .SDcols = char_cols]
selection_final[, (to_num_cols) := lapply(.SD, as.numeric), .SDcols = to_num_cols]


###############################################################################
# 9. Sampling
###############################################################################

## NOTE: Added more robust sampling
need_dissimilarity <- sampling_method %in% c("kmedoids", "maxmin", "hybrid") ||
        (sampling_method == "custom" && average_samples > 0)

dissimilarity_matrix <- NULL
row_similarity <- NULL
if (need_dissimilarity) {
        gower_dist <- daisy(selection_final[, ..distance_variables], metric = "gower")
        dissimilarity_matrix <- as.matrix(gower_dist) # 0..1 dissimilarity
        similarity_matrix <- 1 - dissimilarity_matrix
        row_similarity <- rowMeans(similarity_matrix, na.rm = TRUE)
        names(row_similarity) <- selection_final$ncbi_genbank_assembly_accession
}


### To print
# global_similarity <- mean(distance)
#
# print(paste0("Average sample similarity: ", round(global_similarity, digits = 2)))

sample_custom <- function(selection_final, row_similarity, counts) {
        # counts: named vector c(average=..., outliers=..., random=...)
        if (is.null(row_similarity)) stop("row_similarity required for custom sampling")
        ids <- names(row_similarity)
        # Most similar (average)
        n_avg <- counts["average"]
        n_out <- counts["outliers"]
        n_rand <- counts["random"]
        avg_ids <- if (n_avg > 0) names(sort(row_similarity, decreasing = TRUE))[seq_len(min(n_avg, length(row_similarity)))] else character(0)
        remaining <- setdiff(ids, avg_ids)
        out_ids <- if (n_out > 0 && length(remaining) > 0) names(sort(row_similarity[remaining], decreasing = FALSE))[seq_len(min(n_out, length(remaining)))] else character(0)
        remaining2 <- setdiff(remaining, out_ids)
        rand_ids <- if (n_rand > 0 && length(remaining2) > 0) sample(remaining2, size = min(n_rand, length(remaining2)), replace = FALSE) else character(0)
        selected <- unique(c(avg_ids, out_ids, rand_ids))
        list(selected_ids = selected, groups = list(average = avg_ids, outliers = out_ids, random = rand_ids))
}
sample_kmedoids <- function(selection_final, dissimilarity_matrix, k) {
        if (is.null(dissimilarity_matrix)) stop("dissimilarity_matrix required for kmedoids")
        # cluster::pam accepts a dissimilarity object or matrix (diss=TRUE)
        pam_res <- pam(as.dist(dissimilarity_matrix), k = min(k, nrow(selection_final)), diss = TRUE)
        medoid_idx <- pam_res$id.med
        ids <- selection_final$ncbi_genbank_assembly_accession[medoid_idx]
        list(selected_ids = ids, groups = list(medoids = ids, pam = pam_res))
}
sample_maxmin <- function(selection_final, dissimilarity_matrix, k, seed = 1256) {
        set.seed(seed)
        ids <- selection_final$ncbi_genbank_assembly_accession
        n <- length(ids)
        k <- min(k, n)
        # start with a random seed or most central (choose central for reproducibility)
        central_idx <- which.min(rowMeans(dissimilarity_matrix, na.rm = TRUE))
        selected_idx <- central_idx
        while (length(selected_idx) < k) {
                remaining_idx <- setdiff(seq_len(n), selected_idx)
                # distance of each remaining to nearest selected
                nearest_dist <- sapply(remaining_idx, function(i) min(dissimilarity_matrix[i, selected_idx], na.rm = TRUE))
                next_idx <- remaining_idx[which.max(nearest_dist)]
                selected_idx <- c(selected_idx, next_idx)
        }
        selected_ids <- ids[selected_idx]
        list(selected_ids = selected_ids, groups = list(maxmin = selected_ids))
}
sample_stratified <- function(selection_final, strata_cols = "ncbi_genome_category", total_n, seed = 1256) {
        set.seed(seed)
        if (!all(strata_cols %in% names(selection_final))) stop("strata columns not present")
        sel <- copy(selection_final)
        # create a single stratum key by pasting columns
        sel[, .stratum := do.call(paste, c(.SD, sep = "__")), .SDcols = strata_cols]
        strata_sizes <- sel[, .N, by = .stratum]
        # allocate proportional counts
        strata_sizes[, alloc := floor(N / sum(N) * total_n)]
        remainder <- total_n - sum(strata_sizes$alloc)
        if (remainder > 0) {
                # give extra to largest strata
                ord <- order(strata_sizes$N, decreasing = TRUE)
                strata_sizes$alloc[ord[seq_len(remainder)]] <- strata_sizes$alloc[ord[seq_len(remainder)]] + 1
        }
        selected_ids <- character(0)
        for (i in seq_len(nrow(strata_sizes))) {
                st <- strata_sizes$.stratum[i]
                k <- strata_sizes$alloc[i]
                rows <- sel[.stratum == st]
                if (k > 0 && nrow(rows) > 0) {
                        choose_n <- min(k, nrow(rows))
                        selected_ids <- c(selected_ids, sample(rows$ncbi_genbank_assembly_accession, choose_n))
                }
        }
        list(selected_ids = unique(selected_ids), groups = list(stratified = selected_ids))
}
sample_hybrid <- function(selection_final, dissimilarity_matrix, strata_cols, per_stratum_k, seed = 1256) {
        set.seed(seed)
        sel <- copy(selection_final)
        sel[, .stratum := do.call(paste, c(.SD, sep = "__")), .SDcols = strata_cols]
        selected_ids <- character(0)
        groups <- list()
        for (st in unique(sel$.stratum)) {
                subset_idx <- which(sel$.stratum == st)
                sub_n <- length(subset_idx)
                k <- min(per_stratum_k, sub_n)
                if (k <= 0) next
                sub_diss <- dissimilarity_matrix[subset_idx, subset_idx, drop = FALSE]
                res <- sample_maxmin(sel[subset_idx], sub_diss, k = k, seed = seed)
                selected_ids <- c(selected_ids, res$selected_ids)
                groups[[st]] <- res$selected_ids
        }
        list(selected_ids = unique(selected_ids), groups = groups)
}
# counts parsed earlier: average_samples,outliers_samples,random_samples
counts <- c(average = average_samples, outliers = outliers_samples, random = random_samples)

sampling_result <- switch(sampling_method,
        "custom" = sample_custom(selection_final, row_similarity, counts),
        "kmedoids" = {
                k <- n_samples
                sample_kmedoids(selection_final, dissimilarity_matrix, k)
        },
        "maxmin" = {
                k <- n_samples
                sample_maxmin(selection_final, dissimilarity_matrix, k, seed = 1256)
        },
        "stratified" = {
                # TODO choose strata col(s); could add CLI option e.g. --strata "ncbi_country"
                strata_cols <- "ncbi_genome_category"
                sample_stratified(selection_final, strata_cols, n_samples, seed = 1256)
        },
        "hybrid" = {
                strata_cols <- "ncbi_genome_category"
                per_stratum_k <- max(1, floor(n_samples / length(unique(selection_final[[strata_cols]]))))
                sample_hybrid(selection_final, dissimilarity_matrix, strata_cols, per_stratum_k, seed = 1256)
        },
        stop("Unknown sampling_method: ", sampling_method)
)

sample_list <- sampling_result$selected_ids
dataset_list <- subset(selection_final, ncbi_genbank_assembly_accession %in% sample_list)

## Final list
selection <- dataset_list[, c("ncbi_genbank_assembly_accession")]

## Accession write (for datasets tool)
fwrite(dataset_list, "sample_metadata.csv", col.names = TRUE)
fwrite(selection, "download_accession_list.txt", col.names = FALSE)


###############################################################################
# 11. Genome download
###############################################################################
if (dryrun == TRUE) {
        dataset_list
        cat("Dry run mode quitting")
        quit(status = 0)
}
system("./datasets download genome accession --inputfile download_accession_list.txt --dehydrated --include genome")
unzip("ncbi_dataset.zip")
system("./datasets rehydrate --directory .")
dir.create(output_dir)
suppressMessages(file.rename(
        from = list.files(path = "ncbi_dataset/data/", pattern = "\\.fna$", recursive = TRUE, full.names = TRUE),
        to = file.path(paste0("./", output_dir), basename(list.files("ncbi_dataset/data/", "\\.fna$", recursive = TRUE)))
))
###############################################################################
# 12. Clean-up
###############################################################################
rm_list <- c("download_accession_list.txt", "ncbi_dataset.zip", "md5sum.txt", "README.md")
suppressMessages(file.remove(rm_list))
suppressMessages(system("rm -r ncbi_dataset"))

print("Download completed! \n", "output directory:", output_dir)
