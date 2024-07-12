### format.pvalue.R ################################################################################
# Cleanly format p-value to show on plot 

### format.pvalue ##################################################################################
#' Format p-value to show on plot
#'
#' This function formats a p-value for so that it is easy to read on a plot.
#'
#' @param x value to be formatted
#' @param digits number of digits to show
#' @param lower.value lowest value to show for scientific notation
#' @return expression containing formatted p-value
#' @examples
#' format.pvalue(0.092);
#' format.pvalue(0.000006231);
format.pvalue <- function(x, digits = 1, lower.value = 1e-8) {

	# determine if value of x is below our lower limit
	symbol <- ' = ';
	if (x < lower.value) {
		symbol <- ' < ';
		x <- lower.value;
		}

	# determine if scientific notation is required
	if (x > 0.01) {
		to.return <- as.expression(paste('p', '=', signif(x, 2)));
		} else {
		exponent <- floor(log10(x));
		base <- sprintf(paste("%.", digits, "f", sep = ""), x/10^exponent);

		to.return <- as.expression(substitute(
			expr = paste('p', symbol, base %*% 10^exponent),
			env = list(base = base, exponent = exponent, symbol = symbol)
			));
		}

	return(to.return);
	}
