#' Plot weight series with confidence interval
#' 
#' Plots time series of predicted weights and conf intervals
#' 
#' @param time vector of years
#' @param predictions list with elements fit (estimate) and se (std error) 
#' @param date string for title for date weights are predicted
#' @param ... other arguments for plot
#' @return None
#' @export
#' @author Jeff Laake
plot.weight.series=function (time, predictions, date,...)
{
	plotCI(time, predictions$fit, 1.96 * predictions$se,
			1.96 * predictions$se, xlab = "Cohort",
			ylab = paste("Predicted average weight (kg)",date),
			...)
	lines(time, predictions$fit)
	invisible()
}
