#' Interactive exploration and CMA computation.
#'
#' Interactively plot a given patient's data, allowing the real-time exploration
#' of the various CMAs and their parameters.
#' It can use \code{Rstudio}'s \code{manipulate} library or \code{Shiny}.
#'
#' This is merely a stub for the actual implementation in package
#' \code{AdhereRViz}: it just checks if this package is installed and functional,
#' in which case it calls the actual implementation, otherwise warns the user that
#' \code{AdhereRViz} must be instaled.
#'
#' @seealso Function \code{\link[AdhereR]{plot_interactive_cma}} in package
#' \code{AdhereRViz}.
#'
#' @param ... Parameters to be passed to \code{plot_interactive_cma()} in package
#' \code{AdhereRViz}.
#'
#' @return Nothing
#' @examples
#' \dontrun{
#' plot_interactive_cma(med.events,
#'                      ID.colname="PATIENT_ID",
#'                      event.date.colname="DATE",
#'                      event.duration.colname="DURATION",
#'                      event.daily.dose.colname="PERDAY",
#'                      medication.class.colname="CATEGORY");}
#' @export
plot_interactive_cma <- function(...)
{
  if( requireNamespace("AdhereRViz", quietly=TRUE) )
  {
    # Pass the parameters to AdhereRViz:
    AdhereRViz::plot_interactive_cma(...);
  } else {
    .report.ewms("Package 'AdhereRViz' must be installed for the interactive plotting to work! Please either install it or use the 'normal' plotting functions provided by 'AdhereR'...\n", "error", "plot_interactive_cma", "AdhereR");
    if( interactive() )
    {
      if( menu(c("Yes", "No"), graphics=FALSE, title="Do you want to install 'AdhereRViz' now?") == 1 )
      {
        # Try to install AdhereRViz:
        install.packages("AdhereRViz", dependencies=TRUE);
        if( requireNamespace("AdhereRViz", quietly=TRUE) )
        {
          # Pass the parameters to AdhereRViz:
          AdhereRViz::plot_interactive_cma(...);
        } else
        {
          .report.ewms("Failed to install 'AdhereRViz'!\n", "error", "plot_interactive_cma", "AdhereR");
          return (invisible(NULL));
        }
      } else
      {
        return (invisible(NULL));
      }
    }
  }
}