#' Access the actual CMA estimate from a CMA object.
#'
#' Retrieve the actual CMA estimate(s) encapsulated in a simple, per episode,
#' or sliding window CMA object.
#'
#' @param x a CMA object.
#' @param flatten.medication.groups \emph{Logical}, if \code{TRUE} and there are
#' medication groups defined, then the return value is flattened to a single
#' \code{data.frame} with an extra column containing the medication group (its
#' name is given by \code{medication.groups.colname}).
#' @param medication.groups.colname a \emph{string} (defaults to ".MED_GROUP_ID")
#' giving the name of the column storing the group name when
#' \code{flatten.medication.groups} is \code{TRUE}.
#' @return a \emph{data.frame} containing the CMA estimate(s).
#' @examples
#' cma1 <- CMA1(data=med.events,
#'              ID.colname="PATIENT_ID",
#'              event.date.colname="DATE",
#'              event.duration.colname="DURATION",
#'              followup.window.start=30,
#'              observation.window.start=30,
#'              observation.window.duration=365,
#'              date.format="%m/%d/%Y"
#'             );
#' getCMA(cma1);
#' \dontrun{
#' cmaE <- CMA_per_episode(CMA="CMA1",
#'                         data=med.events,
#'                         ID.colname="PATIENT_ID",
#'                         event.date.colname="DATE",
#'                         event.duration.colname="DURATION",
#'                         event.daily.dose.colname="PERDAY",
#'                         medication.class.colname="CATEGORY",
#'                         carry.only.for.same.medication=FALSE,
#'                         consider.dosage.change=FALSE,
#'                         followup.window.start=0,
#'                         observation.window.start=0,
#'                         observation.window.duration=365,
#'                         date.format="%m/%d/%Y"
#'                        );
#' getCMA(cmaE);}
#' @export
getCMA <- function(x, flatten.medication.groups=FALSE, medication.groups.colname=".MED_GROUP_ID") UseMethod("getCMA")
#' @export
getCMA.CMA0 <- function(x, flatten.medication.groups=FALSE, medication.groups.colname=".MED_GROUP_ID")
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA0") || !("CMA" %in% names(cma)) || is.null(cma$CMA) ) return (NULL);
  if( inherits(cma$CMA, "data.frame") || !flatten.medication.groups )
  {
    return (cma$CMA);
  } else
  {
    # Flatten the medication groups into a single data.frame:
    ret.val <- do.call(rbind, cma$CMA);
    if( is.null(ret.val) || nrow(ret.val) == 0 ) return (NULL);
    ret.val <- cbind(ret.val, unlist(lapply(1:length(cma$CMA), function(i) if(!is.null(cma$CMA[[i]])){rep(names(cma$CMA)[i], nrow(cma$CMA[[i]]))}else{NULL})));
    names(ret.val)[ncol(ret.val)] <- medication.groups.colname; rownames(ret.val) <- NULL;
    return (ret.val);
  }
}