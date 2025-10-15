#!/usr/bin/env Rscript
###############################################################################
# GTDB Dataset creation tool 
#   * Selects taxonomic subset from GTDB metadata
#   * Applies quality filters with a relaxation strategy
#   * Samples genomes with multiple advanced strategies
#   * Downloads FASTA files via NCBI datasets CLI
###############################################################################

# ##############################################################################
# 1. Setup and Dependencies
# ##############################################################################
gtdb_release <- "226" # NOTE: Update this for future GTDB releases

suppressMessages({
    packages <- c("data.table", "R.utils", "cluster", "optparse")
    for (package in packages) {
        if (!require(package, character.only = TRUE, quietly = TRUE)) {
            install.packages(package)
            library(package, character.only = TRUE)
        }
    }
})

## Verification if they are loaded
all_loaded <- all(sapply(packages, require, character.only = TRUE, quietly = TRUE))
if (all_loaded) {
    message("All required packages are loaded and ready.")
} else {
    # Find which ones failed to load
    failed_packages <- packages[!sapply(required_packages, require, character.only = TRUE, quietly = TRUE)]
    stop(paste("Error: Failed to load the following packages:", 
               paste(failed_packages, collapse = ", "),
               "Please check your R environment and internet connection."), 
         call. = FALSE)
}


# ##############################################################################
# 2. Command-Line Interface Definition (using optparse)
# ##############################################################################
option_list <- list(
    # Core required options
    make_option(c("-n", "--name"), type = "character", 
        help = "Required: Target taxon name (e.g., 'Pseudomonas aeruginosa')"),
    make_option(c("--n_samples"), type = "integer", default = 20, 
        help = "Number of genomes to sample [default: %default]"),
    make_option(c("--sampling_method"), type = "character", default = "custom", 
        help = "Sampling method: custom|kmedoids|maxmin|stratified|hybrid [default: %default]"),

    # Filtering and selection options
    make_option(c("-t", "--taxrank"), type = "character", default = "species", 
        help = "Taxonomic rank: species|genus|family|class|phylum|domain [default: %default]"),
    make_option(c("-q", "--quality"), type = "character", default = "good", 
        help = "Quality level: good|bad|mixed [default: %default]"),
    make_option(c("-r", "--representative"), type = "character", default = "f", 
        help = "Only include representative genomes (t/f) [default: %default]"),
    make_option(c("--min_size"), type = "numeric", default = 0, 
        help = "Minimum genome size in base pairs [default: %default]"),
    make_option(c("--max_size"), type = "numeric", default = Inf, 
        help = "Maximum genome size in base pairs [default: %default]"),

    # Custom sampling and data source
    make_option(c("-s", "--sampling_proportion"), type = "character", default = "A:0.5,O:0.25,R:0.25", 
        help = "Proportion schema for 'custom' method (A:average, O:outlier, R:random) [default: '%default']"),
    make_option(c("--strata_cols"), type = "character", default = "ncbi_genome_category", 
        help = "Column(s) for stratified/hybrid sampling, comma-separated [default: '%default']"),
    make_option(c("--database"), type = "character", default = NULL, 
        help = "Path to a custom GTDB metadata TSV file (optional)"),
    make_option(c("-m", "--domain"), type = "character", default = "bacteria", 
        help = "Domain: bacteria|archaea [default: %default]"),

    # Operational options
    make_option(c("-o", "--output_dir"), type = "character", default = "output", 
        help = "Output directory for FASTA files [default: %default]"),
    make_option(c("-d", "--dryrun"), action = "store_true", default = FALSE, 
        help = "Run script without downloading genomes")
)

parser <- OptionParser(option_list = option_list, description = "GTDB Dataset Creation Tool")
args <- parse_args(parser)

# --- Argument Validation ---
#if (is.null(args$name)) {
#    print_help(parser)
#    stop("Error: --name is a required argument.", call. = FALSE)
#}
#valid_ranks <- c("species", "genus", "family", "class", "phylum", "domain")
#if (!tolower(args$taxrank) %in% valid_ranks) stop("Invalid --taxrank. Must be one of: ", paste(valid_ranks, collapse = ", "))
#if (!tolower(args$quality) %in% c("good", "bad", "mixed")) stop("Invalid --quality. Must be 'good', 'bad', or 'mixed'.")
#if (!tolower(args$representative) %in% c("t", "f")) stop("Invalid --representative. Must be 't' or 'f'.")
#if (!tolower(args$domain) %in% c("bacteria", "archaea")) stop("Invalid --domain. Must be 'bacteria' or 'archaea'.")

# Sampling proprotion
tryCatch({
    avg_prop <- as.numeric(gsub(".*A:([0-9.]+).*", "\\1", args$sampling_proportion))
    out_prop <- as.numeric(gsub(".*O:([0-9.]+).*", "\\1", args$sampling_proportion))
    rand_prop <- as.numeric(gsub(".*R:([0-9.]+).*", "\\1", args$sampling_proportion))
    if(sum(avg_prop, out_prop, rand_prop) != 1) warning("Sampling proportions do not sum to 1.")
    
    average_samples <- floor(avg_prop * args$n_samples)
    outliers_samples <- floor(out_prop * args$n_samples)
    random_samples <- args$n_samples - (average_samples + outliers_samples) # Remainder goes to random
}, error = function(e) {
    stop("Could not parse --sampling_proportion. Ensure it follows the 'A:x,O:y,R:z' format.", call. = FALSE)
})


# ##############################################################################
# 3. Helper Functions
# ##############################################################################

#' Splits GTDB taxonomy string into separate columns.
#' @param dt A data.table object.
#' @return The modified data.table with new taxonomy columns.
tax_split <- function(dt) {
    tax_levels <- c("domain", "phylum", "class", "order", "family", "genus", "species") # BUG FIX: Corrected "phylium" to "phylum"
    dt[, (tax_levels) := tstrsplit(gtdb_taxonomy, ";", fixed = TRUE)[1:7]]
    for (col in tax_levels) {
        dt[, (col) := sub("^[dpcofgs]__", "", get(col))]
    }
    return(dt)
}

#' Filters genomes based on quality and iteratively relaxes thresholds if needed.
#' @param data The input data.table.
#' @param quality_mode The quality setting ('good' or 'bad').
#' @param n_target The target number of samples.
#' @return A data.table of quality-filtered genomes.
filter_and_relax <- function(data, quality_mode, n_target) {
    # REFACTOR: Centralized quality filtering logic
    if (quality_mode == "good") {
        # Lower contamination/contigs and higher completeness are better
        contam_q <- 0.25; compl_q <- 0.75; contig_q <- 0.25
        contam_op <- `<=`; compl_op <- `>=`; contig_op <- `<=`
        contam_relax <- 1.05; compl_relax <- 0.95; contig_relax <- 1.05
    } else if (quality_mode == "bad") {
        # Higher contamination/contigs and lower completeness are worse
        contam_q <- 0.75; compl_q <- 0.25; contig_q <- 0.75
        contam_op <- `>=`; compl_op <- `<=`; contig_op <- `>=`
        contam_relax <- 0.95; compl_relax <- 1.05; contig_relax <- 0.95 # Relax towards the median
    }

    # Initial thresholds
    contam_thresh <- as.numeric(quantile(data$checkm2_contamination, contam_q, na.rm = TRUE))
    compl_thresh <- as.numeric(quantile(data$checkm2_completeness, compl_q, na.rm = TRUE))
    contig_thresh <- as.numeric(quantile(data$contig_count, contig_q, na.rm = TRUE))
    
    filtered_data <- data[
        contam_op(checkm2_contamination, contam_thresh) &
        compl_op(checkm2_completeness, compl_thresh) &
        contig_op(contig_count, contig_thresh)
    ]
    
    # Relaxation loop
    relax_steps <- 0
    while (nrow(filtered_data) < n_target && relax_steps < 5) {
        relax_steps <- relax_steps + 1
        message(sprintf("Relaxing thresholds (step %d)...", relax_steps))
        contam_thresh <- contam_thresh * contam_relax
        compl_thresh <- compl_thresh * compl_relax
        contig_thresh <- contig_thresh * contig_relax
        
        filtered_data <- data[
            contam_op(checkm2_contamination, contam_thresh) &
            compl_op(checkm2_completeness, compl_thresh) &
            contig_op(contig_count, contig_thresh)
        ]
        message(sprintf(
            "   Contamination %s %.2f, Completeness %s %.2f, Contigs %s %.0f. Found %d samples.",
            deparse(substitute(contam_op)), contam_thresh,
            deparse(substitute(compl_op)), compl_thresh,
            deparse(substitute(contig_op)), contig_thresh,
            nrow(filtered_data)
        ))
    }
    return(filtered_data)
}


# 4. Console Header
# ==============================================================================
cat("\n\033[1mGTDB Dataset Builder Parameters\033[0m\n")
cat(rep("-", 60), "\n", sep = "")
cat("\033[1;34mCore Settings:\033[0m\n")
cat(sprintf("  %-25s: %s\n", "Taxonomic name", args$name))
cat(sprintf("  %-25s: %s\n", "Taxonomic rank", args$taxrank))
cat(sprintf("  %-25s: %s\n", "Domain", args$domain))
cat(sprintf("  %-25s: %s\n", "Output directory", args$output_dir))
cat(sprintf("  %-25s: %s\n", "Representative genomes", ifelse(args$representative == "t", "Yes", "No")))
cat("\n\033[1;34mSampling Strategy:\033[0m\n")
cat(sprintf("  %-25s: %d\n", "Total samples", args$n_samples))
cat(sprintf("  %-25s: %s\n", "Quality level", args$quality))
cat(sprintf("  %-25s: %s\n", "Sampling method", args$sampling_method))
if (args$sampling_method == "custom") {
    cat(sprintf("  %-25s: %s\n", "Sampling schema", args$sampling_proportion))
    cat(sprintf("  %-25s: %d\n", " -> Average genomes", average_samples))
    cat(sprintf("  %-25s: %d\n", " -> Outlier genomes", outliers_samples))
    cat(sprintf("  %-25s: %d\n", " -> Random genomes", random_samples))
}
cat(rep("-", 60), "\n\n", sep = "")

# ##############################################################################
# 5. Data Download 
# ##############################################################################
## NOTE: Download Database (gtdb_release at the top)
if (!is.null(args$database)) {
    if (!file.exists(args$database)) stop("Specified database file not found: ", args$database)
    cat("Loading custom database from:", args$database, "\n")
    gtdb_data <- fread(args$database)
} else {
    if (args$domain == "bacteria") {
        db_file <- paste0("bac120_metadata_r", sub("\\.", "", gtdb_release), ".tsv")
        db_url <- paste0("https://data.ace.uq.edu.au/public/gtdb/data/releases/release", sub("\\.0", "", gtdb_release), "/", gtdb_release, "/", db_file, ".gz")
    } else {
        db_file <- paste0("ar53_metadata_r", sub("\\.", "", gtdb_release), ".tsv")
        db_url <- paste0("https://data.ace.uq.edu.au/public/gtdb/data/releases/release", sub("\\.0", "", gtdb_release), "/", gtdb_release, "/", db_file, ".gz")
    }

    if (!file.exists(db_file)) {
        cat("Downloading", args$domain, "metadata from GTDB release", gtdb_release, "...\n")
        gtdb_data <- fread(db_url, header = TRUE)
        fwrite(gtdb_data, db_file, sep = "\t")
        cat("Saved database to:", db_file, "\n")
    } else {
        cat("GTDB", args$domain, "database found, loading:", db_file, "\n")
        gtdb_data <- fread(db_file)
    }
}

## NOTE: Download tool
if (!file.exists("datasets")) {
    cat("Downloading NCBI datasets CLI tool...\n")
    download.file(
        "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets", "datasets", quiet = TRUE
    )
    system("chmod +x datasets")
    cat("NCBI datasets CLI tool installed.\n")
} else {
    cat("NCBI datasets CLI tool found.\n")
}

# ##############################################################################
# 6. Metadata Processing and Filtering
# ##############################################################################
cat("Processing database, this may take a moment...\n")
processed_data <- tax_split(gtdb_data)

selection <- processed_data[
    get(args$taxrank) == args$name &
    gtdb_representative == args$representative &
    genome_size >= args$min_size &
    genome_size <= args$max_size
]

if (nrow(selection) == 0) {
    stop("No genomes match the specified criteria. Exiting.")
}
if (args$n_samples > nrow(selection)) {
    warning("Requested samples > available genomes. Using all ", nrow(selection), " entries.")
    args$n_samples <- nrow(selection)
}

## Quality filtering
if (args$quality == "mixed") {
    selection_quality <- selection
    message("Quality 'mixed': No filtering applied. Using all ", nrow(selection), " samples.")
} else {
    selection_quality <- filter_and_relax(selection, args$quality, args$n_samples)
}

if (nrow(selection_quality) == 0) {
    stop("No samples met even the most relaxed quality thresholds. Adjust criteria or data.")
}

if (nrow(selection_quality) < args$n_samples) {
    warning("Only ", nrow(selection_quality), " quality samples available. Returning all.")
    selection_final <- selection_quality
    args$n_samples <- nrow(selection_quality) # Adjust n_samples to what's available
} else {
    selection_final <- selection_quality
}

# ##############################################################################
# 7. Dissimilarity Matrix Calculation
# ##############################################################################
distance_variables <- c(
    "gc_percentage", "genome_size", "checkm2_completeness", "checkm2_contamination",
    "checkm_strain_heterogeneity", "contig_count", "n50_contigs", "protein_count",
    "trna_aa_count", "ncbi_genome_category", "ncbi_isolation_source"
)
distance_variables <- intersect(distance_variables, names(selection_final))
char_cols <- names(selection_final)[sapply(selection_final[, ..distance_variables], is.character)]
if(length(char_cols) > 0) selection_final[, (char_cols) := lapply(.SD, as.factor), .SDcols = char_cols]

dissimilarity_matrix <- NULL
row_similarity <- NULL
need_dissimilarity <- args$sampling_method %in% c("kmedoids", "maxmin", "hybrid") || 
                     (args$sampling_method == "custom" && (average_samples > 0 || outliers_samples > 0))

if (need_dissimilarity) {
    cat("Calculating Gower dissimilarity matrix for sampling...\n")
    gower_dist <- daisy(selection_final[, ..distance_variables], metric = "gower")
    dissimilarity_matrix <- as.matrix(gower_dist)
    similarity_matrix <- 1 - dissimilarity_matrix
    row_similarity <- rowMeans(similarity_matrix, na.rm = TRUE)
    names(row_similarity) <- selection_final$ncbi_genbank_assembly_accession
}

# ##############################################################################
# 8. Sampling
# ##############################################################################

## Sampling Functions 
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

## Execute sampling
counts <- c(average = average_samples, outliers = outliers_samples, random = random_samples)
strata_cols <- trimws(strsplit(args$strata_cols, ",")[[1]])

sampling_result <- switch(args$sampling_method,
    "custom" = sample_custom(selection_final, row_similarity, counts),
    "kmedoids" = sample_kmedoids(selection_final, dissimilarity_matrix, args$n_samples),
    "maxmin" = sample_maxmin(selection_final, dissimilarity_matrix, args$n_samples),
    "stratified" = sample_stratified(selection_final, strata_cols, args$n_samples),
    "hybrid" = {
        per_stratum_k <- max(1, floor(args$n_samples / uniqueN(selection_final, by = strata_cols)))
        sample_hybrid(selection_final, dissimilarity_matrix, strata_cols, per_stratum_k)
    },
    stop("Unknown sampling_method: ", args$sampling_method)
)

final_accessions <- sampling_result$selected_ids
dataset_metadata <- subset(selection_final, ncbi_genbank_assembly_accession %in% final_accessions)

# ##############################################################################
# 9. Output and Download
# ##############################################################################
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("\nSelected %d genomes for download.\n", length(final_accessions)))
fwrite(dataset_metadata,       paste0("./",args$output_dir,"/sample_metadata.csv"        ), col.names = TRUE)
fwrite(list(final_accessions), paste0("./",args$output_dir,"/download_accession_list.txt"), col.names = FALSE)

if (args$dryrun) {
    cat("Dry run mode active. Metadata saved to 'sample_metadata.csv'. No genomes will be downloaded.\n")
    quit(status = 0)
}

## Genome Download
cat("Downloading genomes using NCBI datasets...\n")
system("./datasets download genome accession --inputfile download_accession_list.txt --dehydrated --include genome")
unzip("ncbi_dataset.zip")
system("./datasets rehydrate --directory .")

# Move fastas to output dir
fasta_files <- list.files(path = "ncbi_dataset/data/", pattern = "\\.fna$", recursive = TRUE, full.names = TRUE)
file.rename(from = fasta_files, to = file.path(args$output_dir, basename(fasta_files)))

# ##############################################################################
# 10. Clean-up
# ##############################################################################
cat("Cleaning up temporary files...\n")
files_to_remove <- c("download_accession_list.txt", "ncbi_dataset.zip", "md5checksums.txt", "README.md")
suppressMessages(file.remove(files_to_remove[file.exists(files_to_remove)]))
suppressMessages(unlink("ncbi_dataset", recursive = TRUE))

cat(paste("\n Download complete!\n"))
cat(paste0("FASTA files are located in the '", args$output_dir, "' directory.\n"))
