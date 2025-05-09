#!/usr/bin/env Rscript

suppressMessages(require(data.table))
suppressMessages(require(R.utils))

### Parameters
args <- commandArgs(trailingOnly = TRUE)

### Defaults 
name_selection <- NULL
taxrank <- "species"
representative <- "t"
dryrun <- FALSE
output_dir <- "output"
domain <- "bacteria"  
database_file <- NULL 

# Parse named arguments
i <- 1
while (i <= length(args)) {
  if (args[i] %in% c("--name_selection", "-n")) {
    name_selection <- args[i+1]
    i <- i + 2
  } else if (args[i] %in% c("--taxrank", "-t")) {
    taxrank <- args[i+1]
    i <- i + 2
  } else if (args[i] %in% c("--representative_check", "-r")) {
    representative_check <- args[i+1]
    i <- i + 2
  } else if (args[i] %in% c("--dryrun", "-d")) {
    dryrun <- TRUE
    i <- i + 1
  } else if (args[i] %in% c("--output_dir", "-o")) {
    output_dir <- args[i+1]
    i <- i + 2
  } else if (args[i] %in% c("--domain", "-m")) {  
    domain <- tolower(args[i+1])
    if (!domain %in% c("bacteria", "archaea")) {
      stop("Domain must be either 'bacteria' or 'archaea'")
    }
    i <- i + 2
  } else if (args[i] %in% c("--database", "-db")) {  
    database_file <- args[i+1]
    i <- i + 2
  } else if (args[i] %in% c("--help", "-h")) {
    cat("Usage: fasta_extractor [options]\n")
    cat("Options:\n")
    cat("  -n, --name_selection <name>   Name of the taxon to select (required)\n")
    cat("  -t, --taxrank <rank>         Taxonomic rank [default: species]\n")
    cat("  -r, --representative_check <t/f> Check for representative genomes [default: t]\n")
    cat("  -d, --dryrun                 Dry run (no actual processing) [default: FALSE]\n")
    cat("  -o, --output_dir <dir>       Output directory [default: output]\n")
    cat("  -m, --domain <domain>        Domain (bacteria or archaea) [default: bacteria]\n")
    cat("  -db, --database <file>       Custom database file path\n")
    cat("  -h, --help                   Show this help message\n")
    quit()
  } else {
    # Fallback to positional arguments
    if (is.null(name_selection)) name_selection <- args[i]
    else if (i == 2) taxrank <- args[i]
    else if (i == 3) representative_check <- args[i]
    else if (i == 4) dryrun <- as.logical(args[i])
    i <- i + 1
  }
}

# Check if name_selection was givven
if (is.null(name_selection)) {
  stop("Error: --name_selection argument is required\nUse --help for usage information", call. = FALSE)
}
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
    local_file <- default_files[2]  # Save as uncompressed TSV
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

#################
## NOTE: 2. SRC##
#################

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

########################
## NOTE: 3. File input##
########################
## List prep
processed_data <- tax_split(data)

selection <- processed_data[
        get(taxrank) == name_selection & gtdb_representative == representative, 
]

#######################
## NOTE: 3. Sampling ##
#######################

###FIX change to parameters
n_samples<-20

if(n_samples>nrow(selection)){
print("Selected number of samples is higher than selection!")
        quit(status = 0)
}

### Highquality genomes

n_samples<-20
contamination_thresh<-as.numeric(quantile(selection$checkm2_contamination,0.25))
completness_thresh<-as.numeric(quantile(selection$checkm2_completeness,0.75))
contig_count_thresh<-as.numeric(quantile(selection$contig_count,0.25))

selection_quality<- selection[selection$checkm2_contamination<=contamination&
                              selection$checkm2_completeness>=completness&
                              selection$contig_count<=contig_count
                             ,]


relax_increment <- 0.05  # How much to relax thresholds if too few samples pass
relax_steps <- 0

while (nrow(selection_quality) < n_samples & relax_steps < 4) {  
  relax_steps <- relax_steps + 1
  contamination_thresh <- contamination_thresh * (1 + relax_increment)
  completenes_thresh <- completenes_thresh * (1 - relax_increment)
  contig_count_thresh <- contig_count_thresh * (1 + relax_increment)
  selection_quality <- selection[
    selection$checkm2_contamination <= contamination_thresh &
    selection$checkm2_completeness >= completenes_thresh &
    selection$contig_count <= contig_count_thresh,
  ]
  n_available <- nrow(selection_quality)
  message(
    "Relaxing thresholds (step ", relax_steps, "): ",
    "Contamination ≤ ", round(contamination_thresh, 2), ", ",
    "Completeness ≥ ", round(completenes_thresh, 2), ", ",
    "Contigs ≤ ", round(contig_count_thresh, 2), ". ",
    "Now ", nrow(selection_quality), " samples available."
  )
}

# Final sampling
if (nrow(selection_quality) == 0) {
  stop("No samples met even relaxed quality thresholds. Adjust criteria or data.")
} else if (nrow(selection_quality) < n_samples) {
  if (allow_replacement) {
    warning(
      "Only ", nrow(selection_quality), " high-quality samples available. ",
      "Sampling with replacement to reach ", n_samples, "."
    )
    selection_final <- selection_quality[sample(n_available, n_samples, replace = TRUE), ]
  } else {
    warning(
      "Only ", nrow(selection_quality), " high-quality samples available. ",
      "Returning all (fewer than ", n_samples, ")."
    )
    selection_final <- selection_quality
  }
} else {
  selection_final <- selection_quality[sample(nrow(selection_quality), n_samples, replace = FALSE), ]
}

require(cluster)

col_omitt<-c(1, 5, 9, 18, 19, 20, 21, 22, 23, 28, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 78, 79, 80, 82, 83, 86, 87, 88, 91, 94, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 114, 115, 116, 117, 118, 119)
selection_final<-selection_final[, .SD, .SDcols = -col_omitt]

gower_dist <- daisy(selection_final, 
                    metric = "gower"
                    )
distance<- 1 -  as.matrix(gower_dist)
row_similarity <- rowMeans(distance)
### To print
global_similarity<-mean(distance)

## User confirmation and sample count check
sample_sel<-order(row_similarity, decreasing = TRUE)[1:n_samples]











## User confirmation and sample count check
n_records <- nrow(selection)
unique_species <- processed_data[, .N, by = get(taxrank)][order(-N)]
colnames(unique_species) <- c("Name", "N")
unique_species[, distance := adist(Name, name_selection)]
setorder(unique_species, distance)
top_similar <- unique_species[Name != name_selection][1:min(10, .N)]

if (dryrun == TRUE) {
        cat(paste0("\n", name_selection, " ", n_records, " records found.\n"))
        if (nrow(top_similar) > 0) {
                cat("\nSimilar records:\n(Not selected for download)\n")
                print(top_similar[, .(Name = Name, Records = N)])
                quit(status = 0)
        } else {
                cat("\nNo similar species found\n")
        }
}

cat(paste0("\n", name_selection, " ", n_records, " records found.\n"))
if (nrow(top_similar) > 0) {
        cat("\nSimilar records:\n(Not selected for download)\n")
        print(top_similar[, .(Name = Name, Records = N)])
} else {
        cat("\nNo similar species found\n")
}
cat("Would you like to proceed with the download? [y/n]: ")
response <- readLines("stdin", n = 1)
if (!tolower(response) %in% c("y", "yes")) {
        cat("Aborting at user request.\n")
        quit(status = 0)
}

## Accession write (for datasets tool)
fwrite(selection, "download_accession_list.txt", col.names = FALSE)

## Genomes download
system("./datasets download genome accession --inputfile download_accession_list.txt --dehydrated --include genome")
unzip("ncbi_dataset.zip")
system("./datasets rehydrate --directory .")
dir.create(output_dir)
suppressMessages(file.rename(
        from = list.files(path = "ncbi_dataset/data/", pattern = "\\.fna$", recursive = TRUE, full.names = TRUE),
        to = file.path(paste0("./", output_dir), basename(list.files("ncbi_dataset/data/", "\\.fna$", recursive = TRUE)))
))

## Cleanup
rm_list <- c("download_accession_list.txt", "ncbi_dataset.zip", "md5sum.txt", "README.md")
suppressMessages(file.remove(rm_list))
suppressMessages(system("rm -r ncbi_dataset"))

print("Download completed")
