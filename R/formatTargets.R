### formatTargets.R ################################################################################
# Format mutations - this may be SNV/INDELs (MAF format) or a table containing variant positions
#	with variant information (chr/start/end/ref/alt)

### formatTargets ##################################################################################
#' Load and format the mutations
#'
#' This function loads a SNV/INDEL file (MAF or tab-separated table).
#' Required MAF fields are: Chromosome, Start_Position, End_Position, Reference_Allele,
#' Tumor_Seq_Allele2 (and, for multi-sample files, Tumor_Sample_Barcode).
#' For non-MAF input, required fields are: Chromosome, Start, End, REF, ALT (and,
#' for multi-sample files, Sample).
#' Additional columns are allowed.
#'
#' @param input path to the input SNV/INDEL file (MAF or tab-separated format)
#' @param id sample ID used to select rows from input 
#' @return A data.frame of formatted input
#' @examples
#' test.maf.file <- system.file('extdata', 'test_snp_data.maf', package = 'CompareFragmentSize'); 
#' test.tsv.file <- system.file('extdata', 'test_snp_data.tsv', package = 'CompareFragmentSize'); 
#' test.data.maf <- formatTargets(input = test.maf.file, id = 'Sample1');
#' test.data.tsv <- formatTargets(input = test.tsv.file, id = 'plasma2');
#' @export
formatTargets <- function(input, id = NULL) {

	# read in mutation data
	if (is.null(input)) {
		stop('Please provide path to mutation file.');
		} else {
		input.data <- read.delim(input, stringsAsFactors = FALSE, comment.char = '#');
		}

	# filter mutations to specified sample
	if ('Tumor_Sample_Barcode' %in% colnames(input.data)) {
		colnames(input.data)[which(colnames(input.data) == 'Tumor_Sample_Barcode')] <- 'Sample';
		}

	if (!is.null(id) & ('Sample' %in% colnames(input.data))) {
		input.data <- input.data[which(input.data$Sample == id),];
		}

	# check to make sure there are any mutations left
	if (nrow(input.data) == 0) {
		stop('No target mutations present in "targets" file.');
		}

	# ensure Start/End positions are provided
	# if data is provided as TSV file, these should already be present
	# else, for MAF format
	if ('Start_Position' %in% colnames(input.data)) {
		colnames(input.data)[which(colnames(input.data) == 'Start_Position')] <- 'Start';
		}

	if ('End_Position' %in% colnames(input.data)) {
		colnames(input.data)[which(colnames(input.data) == 'End_Position')] <- 'End';
		}

	# final check to ensure these are present
	if (! 'Start' %in% colnames(input.data)) {
		stop('"Start" must be a column in target mutation file that indicates the start position (1-based coordinates).');
		}

	if (! 'End' %in% colnames(input.data)) {
		stop('"End" must be a column in target mutation file that indicates the end position (1-based coordinates).');
		}

	# ensure REF/ALT allele are provided
	# if data is provided as TSV file, these should already be present
	# else, for MAF format
	if ('Tumor_Seq_Allele2' %in% colnames(input.data)) {
		if (! 'REF' %in% colnames(input.data)) {
			colnames(input.data)[which(colnames(input.data) == 'Reference_Allele')] <- 'REF';
			}
		if (! 'ALT' %in% colnames(input.data)) {
			colnames(input.data)[which(colnames(input.data) == 'Tumor_Seq_Allele2')] <- 'ALT';
			}
		}

	# final check to ensure these are present
	if (! 'REF' %in% colnames(input.data)) {
		stop('"REF" must be a column in target mutation file that indicates the reference allele.');
		}

	if (! 'ALT' %in% colnames(input.data)) {
		stop('"ALT" must be a column in target mutation file that indicates the alternate allele.');
		}

	# indicate the type of variant
	if (! 'Variant_Type' %in% colnames(input.data)) {
		input.data$Variant_Type <- apply(input.data[,c('REF','ALT')],1,function(i) {
			if (nchar(i[1]) > nchar(i[2])) { 'DEL'
			} else if (nchar(i[1]) < nchar(i[2])) { 'INS'
			} else { 'SNP' }
			});
		}

	# indicate fields to return
	keep.fields <- c('Sample','Chromosome','Start','End','REF','ALT','Variant_Type');
	if ('Hugo_Symbol' %in% colnames(input.data)) { keep.fields <- c(keep.fields, 'Hugo_Symbol'); }
	if ('HGVSc' %in% colnames(input.data)) { keep.fields <- c(keep.fields, 'HGVSc'); }

	# return formatted dataframe
	return(input.data[,keep.fields]);
	}

