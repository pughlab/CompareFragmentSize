### annotateReads.R ################################################################################
# Annotate each read as having the expected REF or ALT allele.

### annotateReads ##################################################################################
#' Annotate reads as REF or ALT
#'
#' This function annotates each read in a given GenomicAlignment object as having the expected REF
#' or ALT allele (as provided in targets file).
#'
#' @param ga list containing a GenomicAlignments object (see getReadsFromBAM)
#' @param targets data.frame containing target mutations (see formatTargets)
#' @param refGenome BSgenome to use (required for evaluating insertions; currently supports hg19 and hg38)
#' @param fs_ref vector of fragment scores (ie, from Vessies et al)
#' @return list containing original and annotated GenomicAlignments objects
#' @examples
#' test.maf.file <- system.file('extdata', 'test_snp_data.maf', package = 'CompareFragmentSize'); 
#' test.bam.file <- system.file('extdata', 'sample1_example.bam', package = 'CompareFragmentSize'); 
#' mutations <- formatTargets(input = test.maf.file, id = 'Sample1');
#' target.reads <- getReadsFromBAM(filePath = test.bam.file, targets = mutations[1,]);
#' fragment.data <- annotateReads(ga = target.reads, targets = mutations[1,], refGenome = 'hg19');
#' @export
annotateReads <- function(ga, targets = NULL, fs_ref = NULL, refGenome = "hg38") {

	# ensure target is provided
	if (is.null(targets) || nrow(targets) == 0) {
		stop('No targets provided.');
		}

	if (!all(c('Chromosome','Start','End','REF','ALT','Variant_Type') %in% colnames(targets))) {
		stop("target file must contain fields 'Chromosome', 'Start', 'End', 'REF', 'ALT' and 'Variant_Type'");
		}

	# ensure GA is an output from getReadsFromBAM
	if (!all(c('raw','param') %in% names(ga))) {
		stop('ga must be a GenomicAlignments object generated using getReadsFromBAM');
		}
	if ('GAlignments' != class(ga$raw)) {
		stop('ga must be a GenomicAlignments object generated using getReadsFromBAM');
		}
	
	# if this is an insertion
	if ('INS' == targets[1,]$Variant_Type) {

		alt.size <- nchar(targets[1,]$ALT);

		target.gr <- GenomicRanges::GRanges(
			seqnames = targets[1,]$Chromosome,
			IRanges::IRanges(start = targets[1,]$Start, end = targets[1,]$End)
			);

		targets[1,]$REF <- if ('hg38' == refGenome) {
			as.character(Biostrings::DNAStringSet(IRanges::Views(
			BSgenome.Hsapiens.UCSC.hg38::Hsapiens, target.gr)));
			} else if ('hg19' == refGenome) {
			as.character(Biostrings::DNAStringSet(IRanges::Views(
			BSgenome.Hsapiens.UCSC.hg19::Hsapiens, target.gr)));
			}

	# if this is a deletion
	} else if ('DEL' == targets[1,]$Variant_Type) {

		alt.size <- nchar(targets[1,]$REF);
		targets[1,]$ALT <- paste(rep('-', alt.size), collapse = '');

	# otherwise, it's probably a SNV
	} else {
		alt.size <- 1;
	}

	# extract read information
	sam <- data.frame(ga$raw);
	sam[,c('Group','Allele')] <- NA;

	# determine if reads are REF or ALT
	for (i in 1:nrow(sam)) {

		# ensure full variant is within read
		if (! all(targets$Start[1]:targets$End[1] %in% sam[i,]$start:sam[i,]$end)) {
			next;
			}

		# if the read only contains matches
		if ('150M' == sam[i,]$cigar) {
			sam[i,]$Group <- 'REF';
			next;
			}

		# expand the dna sequence
		dna.seq <- GenomicAlignments::sequenceLayer(
			Biostrings::DNAStringSet(sam[i,]$seq),
			sam[i,]$cigar)[[1]];

		if (sam[i,]$strand == '+') {
			dna.seq <- as.character(dna.seq);
			} else {
			dna.seq <- as.character(Biostrings::complement(dna.seq));
			}

		# how far is it from the beginning of the read and does the allele match?
		# for deletions:
		if ('DEL' == targets[1,]$Variant_Type) {
			pos1 <- targets$Start[1] - sam[i,]$start + 1;
			pos2 <- targets$Start[1] - sam[i,]$start + alt.size;
			allele <- substr(dna.seq,pos1,pos2);

			sam[i,]$Group <- if (allele == targets[1,]$REF) {
				'REF' } else if (allele == targets[1,]$ALT) {
				'ALT' } else { NA }

			# for SNVs:
			} else if (alt.size == 1) {
			pos1 <- targets$Start[1] - sam[i,]$start + 1;
			allele <- substr(dna.seq, pos1, pos1);

			sam[i,]$Group <- if (allele == targets[1,]$REF) {
				'REF' } else if (allele == targets[1,]$ALT) {
				'ALT' } else { NA }

			# for insertions:
			} else {
			dna.seq <- as.character(GenomicAlignments::sequenceLayer(
				Biostrings::DNAStringSet(sam[i,]$seq),
				sam[i,]$cigar, to = 'pairwise')[[1]]);

			pos1 <- targets$Start[1] - sam[i,]$start + 2;
			pos2 <- targets$Start[1] - sam[i,]$start + alt.size + 1;
			allele <- substr(dna.seq,pos1,pos2);

			if ((allele == targets[1,]$ALT) & (grepl(paste0(alt.size,'I'), sam[i,]$cigar))) {
				sam[i,]$Group <- 'ALT';
				} else {
				pos1 <- targets$Start[1] - sam[i,]$start + 1;
				pos2 <- targets$End[1] - sam[i,]$start + 1;
				allele <- substr(dna.seq,pos1,pos2);
				if (allele == targets[1,]$REF) { sam[i,]$Group <- 'REF'; }
				}
			}

		sam[i,]$Allele <- allele;
		}

	# add Group, Allele and insert size annotations to fragment data
	sam$isize <- abs(sam$isize);

	if (!is.null(fs_ref)) {
		sam$FS <- fs_ref[match(sam$isize, 1:length(fs_ref))];
		}

	ga$annotated <- GenomicRanges::makeGRangesFromDataFrame(
		sam[,c('seqnames','strand','start','end','isize','FS','Group','Allele')],
		keep.extra.columns = TRUE
		);

	return(ga);
	}
