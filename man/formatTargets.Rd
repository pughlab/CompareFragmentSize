% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatTargets.R
\name{formatTargets}
\alias{formatTargets}
\title{Load and format the mutations}
\usage{
formatTargets(input, id = NULL)
}
\arguments{
\item{input}{path to the input SNV/INDEL file (MAF or tab-separated format)}

\item{id}{sample ID used to select rows from input}
}
\value{
A data.frame of formatted input
}
\description{
This function loads a SNV/INDEL file (MAF or tab-separated table).
Required MAF fields are: Chromosome, Start_Position, End_Position, Reference_Allele,
Tumor_Seq_Allele2 (and, for multi-sample files, Tumor_Sample_Barcode).
For non-MAF input, required fields are: Chromosome, Start, End, REF, ALT (and,
for multi-sample files, Sample).
Additional columns are allowed.
}
\examples{
test.maf.file <- system.file('extdata', 'test_snp_data.maf', package = 'CompareFragmentSize'); 
test.tsv.file <- system.file('extdata', 'test_snp_data.tsv', package = 'CompareFragmentSize'); 
test.data.maf <- formatTargets(input = test.maf.file, id = 'Sample1');
test.data.tsv <- formatTargets(input = test.tsv.file, id = 'plasma2');
}
