### getReadsFromBAM.R ##############################################################################
# Read in BAM and subset reads to those overlapping target region.

### getReadsFromBAM ################################################################################
#' Load BAM and extract reads
#'
#' This function loads a processed BAM file and subsets reads to thos overlapping target positions.
#'
#' @param filePath path to BAM file
#' @param targets data.frame containing target mutations (see formatTargets)
#' @param mapQthreshold threshold to filter low quality reads (default = 30)
#' @param dedup logical indicating whether to remove duplicates or not (default = FALSE)
#' @return a genomic alignments object
#' @examples
#' test.maf.file <- system.file('extdata', 'test_snp_data.maf', package = 'CompareFragmentSize'); 
#' test.data.maf <- formatTargets(input = test.maf.file, id = 'Sample1');
#' test.bam.file <- system.file('extdata', 'sample1_example.bam', package = 'CompareFragmentSize'); 
#' target.reads <- getReadsFromBAM(filePath = test.bam.file, targets = test.data.maf[1,]);
#' @export
getReadsFromBAM <- function(filePath,  mapQthreshold = 30, targets = NULL, dedup = FALSE) {

	# dedup should be a logical (yes/no to removing duplicates)
	if (!is.logical(dedup)) {
		warning('dedup should be a logical (TRUE/FALSE); setting to default (FALSE)');
		dedup <- FALSE;
		}

	# for scanBamFlag isDuplicate:
	# 	TRUE = keep only duplicate reads
	#	FALSE = keep only un-duplicated reads
	#	NA = keep any read
	remove.duplicates <- ifelse(is.na(dedup) | !dedup, NA, FALSE);

	# ensure target is provided
	if (is.null(targets) || nrow(targets) == 0) {
		stop('No targets provided.');
		}

	# find the BAM index
	indexed.bam <- gsub(".bam$", ".bai", filePath);
	if (!file.exists(indexed.bam)) { 
		indexed.bam <- gsub("$", ".bai", filePath);
		}

	if (!file.exists(indexed.bam)) {
		Rsamtools::indexBam(filePath);
		}

	# indicate parameters to use
	param <- Rsamtools::ScanBamParam(
		mapqFilter = mapQthreshold,
		which = GenomicRanges::makeGRangesFromDataFrame(targets, seqnames.field = 'Chromosome',
			start.field = 'Start', end.field = 'End'),
		what = c('seq','isize'),
		flag = Rsamtools::scanBamFlag(
			isPaired = TRUE,
			isProperPair = TRUE,
			isDuplicate = remove.duplicates,
			isUnmappedQuery = FALSE
			)
		);

	# create genomic alignments object
	ga <- list(
		param = param,
		raw = GenomicAlignments::readGAlignments(filePath, param = param, index = indexed.bam)
		);

	return(ga);
	}
