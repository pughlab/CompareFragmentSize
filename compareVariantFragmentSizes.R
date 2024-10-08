### compareVariantFragmentSizes.R ##################################################################
# Wrapper script to compare variant fragment size on a single sample.
# This function loads a processed BAM file, subsets reads to target positions, annotated each read
# 	as REF or ALT, performs summary statisitcs for each variant and outputs visualizations.

# Example:
# Rscript compareVariantFragmentSizes.R -b /path/to/sample.bam -o /path/to/output/directory -s Sample1 -t /path/to/targets.maf -r hg38

### PREPARE SESSION ################################################################################
# import libraries
library(CompareFragmentSize);
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-b', '--bam', type = 'character', help = 'path to input file (BAM)');
parser$add_argument('-o', '--output_dir', type = 'character', help = 'path to output directory');
parser$add_argument('-s', '--sample', type = 'character', default = NULL,
	help = "sample ID (should match value in 'Sample' or 'Tumor_Sample_Barcode' fields of 'targets')");
parser$add_argument('-t', '--targets', type = 'character', help = 'path to target regions');
parser$add_argument('-r', '--ref', type = 'character', help = 'fragment-length reference (Vessies et al.)',
	default = system.file('extdata', 'vessies_reference_set.txt', package = 'CompareFragmentSize'));
parser$add_argument('-g', '--genome', type = 'character', help = 'reference genome (hg38 or hg19)',
	default = 'hg38');

arguments <- parser$parse_args();

# check for required arguments
if (is.null(arguments$targets)) {
        stop('No target mutations provided; please specify path to target mutations file using --targets.');
        }

if (is.null(arguments$bam)) {
	stop('No BAM file provided; please specify path to input BAM file using --bam.');
	}

if (is.null(arguments$output_dir)) {
	stop('No output directory provided; please specify path to output directory using --output_dir.');
	}

if (is.null(arguments$sample)) {
	warning('No argument to --sample provided; we will assume that every variant is to be checked.');
	}

ignore.insertions <- FALSE;
if (!arguments$genome %in% c('hg19','hg38')) {
	warning("Argument --genome must be one of 'hg19' or 'hg38' in order for insertions to be checked.");
	ignore.insertions <- TRUE;
	}

if (is.null(arguments$ref)) {
	warning("No argument to --ref provided; will not provide fragment scores.");
	}

# make and set an output directory
output.dir <- if (!is.null(arguments$sample)) {
	paste0(arguments$output_dir, '/', arguments$sample);
	} else {
	arguments$output_dir;
	}

if (!dir.exists(output.dir)) {
	dir.create(output.dir, recursive = TRUE);
	}

### MAIN ###########################################################################################
# read in target positions (+ any annotations; ie, Sample or Gene names)
mutations <- formatTargets(arguments$targets, id = arguments$sample);

if (ignore.insertions) {
	mutations <- mutations[which(mutations$Variant_Type != 'INS'),];
	}

if (nrow(mutations) == 0) {
	stop('No appropriate mutations available; halting program.');
	}

# load in fragment score reference (if provided)
vessies.ref <- NULL;
if (!is.null(arguments$ref)) {
	vessies.ref <- as.numeric(read.delim(arguments$ref, header = FALSE)$V1);
	}

# move to output directory
setwd(output.dir);

# initiate list to hold variant sizes
metrics.data <- list();

# loop over each variant
for (idx in 1:nrow(mutations)) {

	# format fragment data
	target.reads <- getReadsFromBAM(
		filePath = arguments$bam,
		targets = mutations[idx,],
		dedup = TRUE
		);

	if (length(target.reads$raw) == 0) { next; }

	target.reads <- annotateReads(
		ga = target.reads,
		targets = mutations[idx,],
		refGenome = arguments$genome,
		fs_ref = vessies.ref
		);

	metrics.data[[idx]] <- collectMetrics(
		fragment.data = target.reads,
		target = mutations[idx,],
		verbose = FALSE
		);

	if (is.na(metrics.data[[idx]]$KS.p)) { next; }

	plotFragmentSizes(
		fragment.data = target.reads,
		id = arguments$sample,
		target = metrics.data[[idx]]
		);
	}

summary.metrics <- do.call(rbind, metrics.data);

# format for output
write.table(
	summary.metrics,
	file = paste0(Sys.Date(), '_fragment_size_summary.tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

