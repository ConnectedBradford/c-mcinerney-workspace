    getMGs <- function(x) UseMethod("getMGs")
#' @export
getMGs.CMA0 <- function(x)
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA0") || is.null(cma$medication.groups) ) return (NULL);
  return (cma$medication.groups);
}


CMA0 <- function(data=NULL, # the data used to compute the CMA on
                 # Important columns in the data
                 ID.colname=NA, # the name of the column containing the unique patient ID (NA = undefined)
                 event.date.colname=NA, # the start date of the event in the date.format format (NA = undefined)
                 event.duration.colname=NA, # the event duration in days (NA = undefined)
                 event.daily.dose.colname=NA, # the prescribed daily dose (NA = undefined)
                 medication.class.colname=NA, # the classes/types/groups of medication (NA = undefined)
                 # Groups of medication classes:
                 medication.groups=NULL, # a named vector of medication group definitions, the name of a column in the data that defines the groups, or NULL
                 flatten.medication.groups=FALSE, medication.groups.colname=".MED_GROUP_ID", # if medication.groups were defined, return CMAs and event.info as single data.frame?
                 # Various types methods of computing gaps:
                 carryover.within.obs.window=NA, # if TRUE consider the carry-over within the observation window (NA = undefined)
                 carryover.into.obs.window=NA, # if TRUE consider the carry-over from before the starting date of the observation window (NA = undefined)
                 carry.only.for.same.medication=NA, # if TRUE the carry-over applies only across medication of same type (NA = undefined)
                 consider.dosage.change=NA, # if TRUE carry-over is adjusted to reflect changes in dosage (NA = undefined)
                 # The follow-up window:
                 followup.window.start=0, # if a number is the earliest event per participant date plus number of units, or a Date object, or a column name in data (NA = undefined)
                 followup.window.start.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!) (NA = undefined)
                 followup.window.start.per.medication.group=FALSE, # if there are medication groups and this is TRUE, then the first event is relative to each medication group separately, otherwise is relative to the patient
                 followup.window.duration=365*2, # the duration of the follow-up window in the time units given below (NA = undefined)
                 followup.window.duration.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!)  (NA = undefined)
                 # The observation window (embedded in the follow-up window):
                 observation.window.start=0, # the number of time units relative to followup.window.start (NA = undefined)
                 observation.window.start.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!) (NA = undefined)
                 observation.window.duration=365*2, # the duration of the observation window in time units (NA = undefined)
                 observation.window.duration.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!) (NA = undefined)
                 # Date format:
                 date.format="%m/%d/%Y", # the format of the dates used in this function (NA = undefined)
                 # Comments and metadata:
                 summary="Base CMA object",
                 # Misc:
                 suppress.warnings=FALSE,
                 suppress.special.argument.checks=FALSE, # used internally to suppress the check that we don't use special argument names
                 arguments.that.should.not.be.defined=NULL, # the list of argument names and values for which a warning should be thrown if passed to the function
                 ...
)
{
  # Preconditions:
  if( !is.null(data) )
  {
    # data's class and dimensions:
    if( inherits(data, "matrix") || inherits(data, "tbl") || inherits(data, "data.table") ) data <- as.data.frame(data); # convert various things to data.frame
    if( !inherits(data, "data.frame") )
    {
      if( !suppress.warnings ) .report.ewms("The 'data' for a CMA must be of type 'data.frame'!\n", "error", "CMA0", "AdhereR");
      return (NULL);
    }
    if( nrow(data) < 1 )
    {
      if( !suppress.warnings ) .report.ewms("The 'data' for a CMA must have at least one row!\n", "error", "CMA0", "AdhereR");
      return (NULL);
    }
    # the column names must exist in data:
    if( !is.na(ID.colname) && !(ID.colname %in% names(data)) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Column ID.colname='",ID.colname,"' must appear in the 'data'!\n"), "error", "CMA0", "AdhereR");
      return (NULL);
    }
    if( !is.na(event.date.colname) && !(event.date.colname %in% names(data)) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Column event.date.colname='",event.date.colname,"' must appear in the 'data'!\n"), "error", "CMA0", "AdhereR");
      return (NULL);
    }
    if( !is.na(event.duration.colname) && !(event.duration.colname %in% names(data)) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Column event.duration.colname='",event.duration.colname,"' must appear in the 'data'!\n"), "error", "CMA0", "AdhereR");
      return (NULL);
    }
    if( !is.na(event.daily.dose.colname) && !(event.daily.dose.colname %in% names(data)) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Column event.daily.dose.colname='",event.daily.dose.colname,"' must appear in the 'data'!\n"), "error", "CMA0", "AdhereR");
      return (NULL);
    }
    if( !is.na(medication.class.colname) && !(medication.class.colname %in% names(data)) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Column medication.class.colname='",medication.class.colname,"' must appear in the 'data'!\n"), "error", "CMA0", "AdhereR");
      return (NULL);
    }
    # carry-over conditions:
    if( !is.na(carryover.within.obs.window) && !is.logical(carryover.within.obs.window) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Parameter 'carryover.within.obs.window' must be logical!\n"), "error", "CMA0", "AdhereR");
      return (NULL);
    }
    if( !is.na(carryover.into.obs.window) && !is.logical(carryover.into.obs.window) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Parameter 'carryover.into.obs.window' must be logical!\n"), "error", "CMA0", "AdhereR");
      return (NULL);
    }
    if( !is.na(carry.only.for.same.medication) && !is.logical(carry.only.for.same.medication) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Parameter 'carry.only.for.same.medication' must be logical!\n"), "error", "CMA0", "AdhereR");
      return (NULL);
    }
    if( !is.na(consider.dosage.change) && !is.logical(consider.dosage.change) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Parameter 'consider.dosage.change' must be logical!\n"), "error", "CMA0", "AdhereR");
      return (NULL);
    }
    if( (!is.na(carryover.within.obs.window) && !is.na(carryover.into.obs.window) && !is.na(carry.only.for.same.medication)) &&
          (!carryover.within.obs.window && !carryover.into.obs.window && carry.only.for.same.medication) )
    {
      if( !suppress.warnings ) .report.ewms("Cannot carry over only for same medication when no carry over at all is considered!\n", "error", "CMA0", "AdhereR")
      return (NULL);
    }
    # the follow-up window:
    if( !is.na(followup.window.start) && !inherits(followup.window.start,"Date") && !is.numeric(followup.window.start) && !(followup.window.start %in% names(data)) )
    {
      # See if it can be forced to a valid Date:
      if( !is.na(followup.window.start) && (is.character(followup.window.start) || is.factor(followup.window.start)) && !is.na(followup.window.start <- as.Date(followup.window.start, format=date.format, optional=TRUE)) )
      {
        # Ok, it was apparently successfully converted to Date: nothing else to do...
      } else
      {
        if( !suppress.warnings ) .report.ewms("The follow-up window start must be either a number, a Date, or a valid column name in 'data'!\n", "error", "CMA0", "AdhereR")
        return (NULL);
      }
    }
    if( !inherits(followup.window.start,"Date") && !is.na(followup.window.start.unit) && !(followup.window.start.unit %in% c("days", "weeks", "months", "years") ) )
    {
      if( !suppress.warnings ) .report.ewms("The follow-up window start unit is not recognized!\n", "error", "CMA0", "AdhereR")
      return (NULL);
    }
    if( !is.na(followup.window.duration) && is.numeric(followup.window.duration) && followup.window.duration < 0 )
    {
      if( !suppress.warnings ) .report.ewms("The follow-up window duration must be a positive number of time units after the first event!\n", "error", "CMA0", "AdhereR")
      return (NULL);
    } else if( !is.na(followup.window.duration) && !is.numeric(followup.window.duration) && !(followup.window.duration %in% names(data)) )
    {
      if( !suppress.warnings ) .report.ewms("The follow-up window duration must be either a positive number or a valid column name in 'data'!\n", "error", "CMA0", "AdhereR")
      return (NULL);
    }
    if( !is.na(followup.window.duration.unit) && !(followup.window.duration.unit %in% c("days", "weeks", "months", "years") ) )
    {
      if( !suppress.warnings ) .report.ewms("The follow-up window duration unit is not recognized!\n", "error", "CMA0", "AdhereR")
      return (NULL);
    }
    # the observation window:
    if( !is.na(observation.window.start) && is.numeric(observation.window.start) && observation.window.start < 0 )
    {
      if( !suppress.warnings ) .report.ewms("The observation window start must be a positive number of time units after the first event!\n", "error", "CMA0", "AdhereR")
      return (NULL);
    } else if( !is.na(observation.window.start) && !inherits(observation.window.start,"Date") && !is.numeric(observation.window.start) && !(observation.window.start %in% names(data)) )
    {
      # See if it can be forced to a valid Date:
      if( !is.na(observation.window.start) && (is.character(observation.window.start) || is.factor(observation.window.start)) && !is.na(observation.window.start <- as.Date(observation.window.start, format=date.format, optional=TRUE)) )
      {
        # Ok, it was apparently successfully converted to Date: nothing else to do...
      } else
      {
        if( !suppress.warnings ) .report.ewms("The observation window start must be either a positive number, a Date, or a valid column name in 'data'!\n", "error", "CMA0", "AdhereR")
        return (NULL);
      }
    }
    if( !inherits(observation.window.start,"Date") && !is.na(observation.window.start.unit) && !(observation.window.start.unit %in% c("days", "weeks", "months", "years") ) )
    {
      if( !suppress.warnings ) .report.ewms("The observation window start unit is not recognized!\n", "error", "CMA0", "AdhereR")
      return (NULL);
    }
    if( !is.na(observation.window.duration) && is.numeric(observation.window.duration) && observation.window.duration < 0 )
    {
      if( !suppress.warnings ) .report.ewms("The observation window duration must be a positive number of time units after the first event!\n", "error", "CMA0", "AdhereR")
      return (NULL);
    } else if( !is.na(observation.window.duration) && !is.numeric(observation.window.duration) && !(observation.window.duration %in% names(data)) )
    {
      if( !suppress.warnings ) .report.ewms("The observation window duration must be either a positive number or a valid column name in 'data'!\n", "error", "CMA0", "AdhereR")
      return (NULL);
    }
    if( !is.na(observation.window.duration.unit) && !(observation.window.duration.unit %in% c("days", "weeks", "months", "years") ) )
    {
      if( !suppress.warnings ) .report.ewms("The observation window duration unit is not recognized!\n", "error", "CMA0", "AdhereR")
      return (NULL);
    }

    # Arguments that should not have been passed:
    if( !suppress.warnings && !is.null(arguments.that.should.not.be.defined) )
    {
      # Get the actual list of arguments (including in the ...); the first is the function's own name:
      args.list <- as.list(match.call(expand.dots = TRUE));
      args.mathing <- (names(arguments.that.should.not.be.defined) %in% names(args.list)[-1]);
      if( any(args.mathing) )
      {
        for( i in which(args.mathing) )
        {
          .report.ewms(paste0("Please note that '",args.list[[1]],"' overrides argument '",names(arguments.that.should.not.be.defined)[i],"' with value '",arguments.that.should.not.be.defined[i],"'!\n"), "warning", "CMA0", "AdhereR");
        }
      }
    }
  } else
  {
    return (NULL);
  }

  # Parse, check and evaluate the medication groups (if any):
  if( !is.null(medication.groups) )
  {
    if( is.null(mg <- .apply.medication.groups(medication.groups=medication.groups, data=data, suppress.warnings=suppress.warnings)) )
    {
      return (NULL);
    }
  } else
  {
    mg <- NULL;
  }

  structure(list("data"=data,
                 "ID.colname"=ID.colname,
                 "event.date.colname"=event.date.colname,
                 "event.duration.colname"=event.duration.colname,
                 "event.daily.dose.colname"=event.daily.dose.colname,
                 "medication.class.colname"=medication.class.colname,
                 "medication.groups"=mg,
                 "flatten.medication.groups"=flatten.medication.groups,
                 "medication.groups.colname"=medication.groups.colname,
                 "carryover.within.obs.window"=carryover.within.obs.window,
                 "carryover.into.obs.window"=carryover.into.obs.window,
                 "carry.only.for.same.medication"=carry.only.for.same.medication,
                 "consider.dosage.change"=consider.dosage.change,
                 "followup.window.start"=followup.window.start,
                 "followup.window.start.unit"=followup.window.start.unit,
                 "followup.window.start.per.medication.group"=followup.window.start.per.medication.group,
                 "followup.window.duration"=followup.window.duration,
                 "followup.window.duration.unit"=followup.window.duration.unit,
                 "observation.window.start"=observation.window.start,
                 "observation.window.start.unit"=observation.window.start.unit,
                 "observation.window.duration"=observation.window.duration,
                 "observation.window.duration.unit"=observation.window.duration.unit,
                 "date.format"=date.format,
                 "summary"=summary),
            class="CMA0");
}


#' Print CMA0 (and derived) objects.
#'
#' Prints and summarizes a basic CMA0, or derived, object.
#'
#' Can produce output for the console (text), R Markdown or LaTeX,
#' showing various types of information.
#'
#' @param x A \emph{\code{CMA0}} or derived object, representing the CMA to
#' print.
#' @param inline \emph{Logical}, should print inside a line of text or as a
#' separate, extended object?
#' @param format A \emph{string}, the type of output: plain text ("text";
#' default), LaTeX ("latex") or R Markdown ("markdown").
#' @param print.params \emph{Logical}, should print the parameters?
#' @param print.data \emph{Logical}, should print a summary of the data?
#' @param exclude.params A vector of \emph{strings}, the names of the object
#' fields to exclude from printing (usually, internal information irrelevant to
#' the end-user).
#' @param skip.header \emph{Logical}, should the header be printed?
#' @param cma.type A \emph{string}, used to override the reported object's class.
#' @param ... other possible parameters
#' @examples
#' cma0 <- CMA0(data=med.events,
#'              ID.colname="PATIENT_ID",
#'              event.date.colname="DATE",
#'              event.duration.colname="DURATION",
#'              event.daily.dose.colname="PERDAY",
#'              medication.class.colname="CATEGORY",
#'              followup.window.start=0,
#'              followup.window.start.unit="days",
#'              followup.window.duration=2*365,
#'              followup.window.duration.unit="days",
#'              observation.window.start=30,
#'              observation.window.start.unit="days",
#'              observation.window.duration=365,
#'              observation.window.duration.unit="days",
#'              date.format="%m/%d/%Y",
#'              summary="Base CMA");
#' cma0;
#' print(cma0, format="markdown");
#' cma1 <- CMA1(data=med.events,
#'              ID.colname="PATIENT_ID",
#'              event.date.colname="DATE",
#'              event.duration.colname="DURATION",
#'              followup.window.start=30,
#'              observation.window.start=30,
#'              observation.window.duration=365,
#'              date.format="%m/%d/%Y"
#'             );
#' cma1;
#' @export
print.CMA0 <- function(x,                                                  # the CMA0 (or derived) object
                       ...,                                                # required for S3 consistency
                       inline=FALSE,                                       # print inside a line of text or as a separate, extended object?
                       format=c("text", "latex", "markdown"),              # the format to print to
                       print.params=TRUE,                                  # show the parameters?
                       print.data=TRUE,                                    # show the summary of the data?
                       exclude.params=c("event.info", "real.obs.windows"), # if so, should I not print some?
                       skip.header=FALSE,                                  # should I print the generic header?
                       cma.type=class(cma)[1]
)
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA0") ) return (invisible(NULL));

  if( format[1] == "text" )
  {
    # Output text:
    if( !inline )
    {
      # Extended print:
      if( !skip.header ) cat(paste0(cma.type,":\n"));
      if( print.params )
      {
        params <- names(cma); params <- params[!(params %in% c("data",exclude.params))]; # exlude the 'data' (and any other requested) params from printing
        if( length(params) > 0 )
        {
          if( "summary" %in% params )
          {
            cat(paste0("  \"",cma$summary,"\"\n"));
            params <- params[!(params %in% "summary")];
          }
          cat("  [\n");
          for( p in params )
          {
            if( p == "CMA" )
            {
              cat(paste0("    ",p," = CMA results for ",nrow(cma[[p]])," patients\n"));
            } else if( p == "medication.groups" )
            {
              if( !is.null(cma[[p]]) )
              {
                cat(paste0("    ", p, " = ", nrow(cma[[p]]$defs), " [", ifelse(nrow(cma[[p]]$defs)<4, paste0("'",cma[[p]]$defs$name,"'", collapse=", "), paste0(paste0("'",cma[[p]]$defs$name[1:4],"'", collapse=", ")," ...")), "]\n"));
              } else
              {
                cat(paste0("    ", p, " = <NONE>\n"));
              }
            } else if( !is.null(cma[[p]]) && length(cma[[p]]) > 0 && !is.na(cma[[p]]) )
            {
              cat(paste0("    ",p," = ",cma[[p]],"\n"));
            }
          }
          cat("  ]\n");
        }
        if( print.data && !is.null(cma$data) )
        {
          # Data summary:
          cat(paste0("  DATA: ",nrow(cma$data)," (rows) x ",ncol(cma$data)," (columns)"," [",length(unique(cma$data[,cma$ID.colname]))," patients]",".\n"));
        }
      }
    } else
    {
      # Inline print:
      cat(paste0(cma$summary,ifelse(print.data && !is.null(cma$data),paste0(" (on ",nrow(cma$data)," rows x ",ncol(cma$data)," columns",", ",length(unique(cma$data[,cma$ID.colname]))," patients",")"),"")));
    }
  } else if( format[1] == "latex" )
  {
    # Output LaTeX: no difference between inline and not inline:
    cat(paste0("\\textbf{",cma$summary,"} (",cma.type,"):",
               ifelse(print.data && !is.null(cma$data),paste0(" (on ",nrow(cma$data)," rows x ",ncol(cma$data)," columns",", ",length(unique(cma$data[,cma$ID.colname]))," patients",")"),"")));
  } else if( format[1] == "markdown" )
  {
    # Output Markdown: no difference between inline and not inline:
    cat(paste0("**",cma$summary,"** (",cma.type,"):",
               ifelse(print.data && !is.null(cma$data),paste0(" (on ",nrow(cma$data)," rows x ",ncol(cma$data)," columns",", ",length(unique(cma$data[,cma$ID.colname]))," patients",")"),"")));
  } else
  {
    .report.ewms("Unknown format for printing!\n", "error", "print.CMA0", "AdhereR");
    return (invisible(NULL));
  }
}

#' Plot CMA0 objects.
#'
#' Plots the events (prescribing or dispensing) data encapsulated in a basic
#' CMA0 object.
#'
#' The x-axis represents time (either in days since the earliest date or as
#' actual dates), with consecutive events represented as ascending on the y-axis.
#'
#' Each event is represented as a segment with style \code{lty.event} and line
#' width \code{lwd.event} starting with a \code{pch.start.event} and ending with
#' a \code{pch.end.event} character, coloured with a unique color as given by
#' \code{col.cats}, extending from its start date until its end date.
#' Consecutive events are thus represented on consecutive levels of the y-axis
#' and are connected by a "continuation" line with \code{col.continuation}
#' colour, \code{lty.continuation} style and \code{lwd.continuation} width;
#' these continuation lines are purely visual guides helping to perceive the
#' sequence of events, and carry no information about the availability of
#' medication in this interval.
#'
#' When several patients are displayed on the same plot, they are organized
#' vertically, and alternating bands (white and gray) help distinguish
#' consecutive patients.
#' Implicitly, all patients contained in the \code{cma} object will be plotted,
#' but the \code{patients.to.plot} parameter allows the selection of a subset
#' of patients.
#'
#' @param x A \emph{\code{CMA0}} or derived object, representing the CMA to
#' plot
#' @param patients.to.plot A vector of \emph{strings} containing the list of
#' patient IDs to plot (a subset of those in the \code{cma} object), or
#' \code{NULL} for all
#' @param duration A \emph{number}, the total duration (in days) of the whole
#' period to plot; in \code{NA} it is automatically determined from the event
#' data such that the whole dataset fits.
#' @param align.all.patients \emph{Logical}, should all patients be aligned
#' (i.e., the actual dates are discarded and all plots are relative to the
#' earliest date)?
#' @param align.first.event.at.zero \emph{Logical}, should the first event be
#' placed at the origin of the time axis (at 0)?
#' @param show.period A \emph{string}, if "dates" show the actual dates at the
#' regular grid intervals, while for "days" (the default) shows the days since
#' the beginning; if \code{align.all.patients == TRUE}, \code{show.period} is
#' taken as "days".
#' @param period.in.days The \emph{number} of days at which the regular grid is
#' drawn (or 0 for no grid).
#' @param show.legend \emph{Logical}, should the legend be drawn?
#' @param legend.x The position of the legend on the x axis; can be "left",
#' "right" (default), or a \emph{numeric} value.
#' @param legend.y The position of the legend on the y axis; can be "bottom"
#' (default), "top", or a \emph{numeric} value.
#' @param legend.bkg.opacity A \emph{number} between 0.0 and 1.0 specifying the
#' opacity of the legend background.
#' @param cex,cex.axis,cex.lab,legend.cex,legend.cex.title \emph{numeric} values
#' specifying the cex of the various types of text.
#' @param xlab Named vector of x-axis labels to show for the two types of periods
#' ("days" and "dates"), or a single value for both, or \code{NULL} for nothing.
#' @param ylab Named vector of y-axis labels to show without and with CMA estimates,
#' or a single value for both, or \code{NULL} for nonthing.
#' @param title Named vector of titles to show for and without alignment, or a
#' single value for both, or \code{NULL} for nonthing.
#' @param col.cats A \emph{color} or a \emph{function} that specifies the single
#' colour or the colour palette used to plot the different medication; by
#' default \code{rainbow}, but we recommend, whenever possible, a
#' colorblind-friendly palette such as \code{viridis} or \code{colorblind_pal}.
#' @param unspecified.category.label A \emph{string} giving the name of the
#' unspecified (generic) medication category.
#' @param medication.groups.to.plot the names of the medication groups to plot or
#' \code{NULL} (the default) for all.
#' @param medication.groups.separator.show a \emph{boolean}, if \code{TRUE} (the
#' default) visually mark the medication groups the belong to the same patient,
#' using horizontal lines and alternating vertical lines.
#' @param medication.groups.separator.lty,medication.groups.separator.lwd,medication.groups.separator.color
#' graphical parameters (line type, line width and colour describing the visual
#' marking og medication groups as beloning to the same patient.
#' @param medication.groups.allother.label a \emph{string} giving the label to
#' use for the implicit \code{__ALL_OTHERS__} medication group (defaults to "*").
#' @param lty.event,lwd.event,pch.start.event,pch.end.event The style of the
#' event (line style, width, and start and end symbols).
#' @param plot.events.vertically.displaced Should consecutive events be plotted
#' on separate rows (i.e., separated vertically, the default) or on the same row?
#' @param print.dose \emph{Logical}, should the daily dose be printed as text?
#' @param cex.dose \emph{Numeric}, if daily dose is printed, what text size
#' to use?
#' @param print.dose.outline.col If \emph{\code{NA}}, don't print dose text with
#' outline, otherwise a color name/code for the outline.
#' @param print.dose.centered \emph{Logical}, print the daily dose centered on
#' the segment or slightly below it?
#' @param plot.dose \emph{Logical}, should the daily dose be indicated through
#' segment width?
#' @param lwd.event.max.dose \emph{Numeric}, the segment width corresponding to
#' the maximum daily dose (must be >= lwd.event but not too big either).
#' @param plot.dose.lwd.across.medication.classes \emph{Logical}, if \code{TRUE},
#' the line width of the even is scaled relative to all medication classes (i.e.,
#' relative to the global minimum and maximum doses), otherwise it is scale
#' relative only to its medication class.
#' @param col.continuation,lty.continuation,lwd.continuation The style of the
#' "continuation" lines connecting consecutive events (colour, line style and
#' width).
#' @param col.na The colour used for missing event data.
#' @param highlight.followup.window \emph{Logical}, should the follow-up window
#' be plotted?
#' @param followup.window.col The follow-up window's colour.
#' @param highlight.observation.window \emph{Logical}, should the observation
#' window be plotted?
#' @param observation.window.col,observation.window.density,observation.window.angle,observation.window.opacity
#' Attributes of the observation window (colour, shading density, angle and
#' opacity).
#' @param alternating.bands.cols The colors of the alternating vertical bands
#' distinguishing the patients; can be \code{NULL} = don't draw the bandes;
#' or a vector of colors.
#' @param bw.plot \emph{Logical}, should the plot use grayscale only (i.e., the
#' \code{\link[grDevices]{gray.colors}} function)?
#' @param rotate.text \emph{Numeric}, the angle by which certain text elements
#' (e.g., axis labels) should be rotated.
#' @param force.draw.text \emph{Logical}, if \code{TRUE}, always draw text even
#' if too big or too small
#' @param min.plot.size.in.characters.horiz,min.plot.size.in.characters.vert
#' \emph{Numeric}, the minimum size of the plotting surface in characters;
#' horizontally (min.plot.size.in.characters.horiz) refers to the the whole
#' duration of the events to plot; vertically (min.plot.size.in.characters.vert)
#' refers to a single event. If the plotting is too small, possible solutions
#' might be: if within \code{RStudio}, try to enlarge the "Plots" panel, or
#' (also valid outside \code{RStudio} but not if using \code{RStudio server}
#' start a new plotting device (e.g., using \code{X11()}, \code{quartz()}
#' or \code{windows()}, depending on OS) or (works always) save to an image
#' (e.g., \code{jpeg(...); ...; dev.off()}) and display it in a viewer.
#' @param suppress.warnings \emph{Logical}: show or hide the warnings?
#' @param max.patients.to.plot \emph{Numeric}, the maximum patients to attempt
#' to plot.
#' @param export.formats a \emph{string} giving the formats to export the figure
#' to (by default \code{NULL}, meaning no exporting); can be any combination of
#' "svg" (just an \code{SVG} file), "html" (\code{SVG} + \code{HTML} + \code{CSS}
#' + \code{JavaScript}, all embedded within one \code{HTML} document), "jpg",
#' "png", "webp", "ps" or "pdf".
#' @param export.formats.fileprefix a \emph{string} giving the file name prefix
#' for the exported formats (defaults to "AdhereR-plot").
#' @param export.formats.height,export.formats.width \emph{numbers} giving the
#' desired dimensions (in pixels) for the exported figure (defaults to sane
#' values if \code{NA}).
#' @param export.formats.save.svg.placeholder a \emph{logical}, if TRUE, save an
#' image placeholder of type given by \code{export.formats.svg.placeholder.type}
#'for the \code{SVG} image.
#' @param export.formats.svg.placeholder.type a \emph{string}, giving the type of
#' placeholder for the \code{SVG} image to save; can be "jpg",
#' "png" (the default) or "webp".
#' @param export.formats.svg.placeholder.embed a \emph{logical}, if \code{TRUE},
#' embed the placeholder image in the HTML document (if any) using \code{base64}
#' encoding, otherwise (the default) leave it as an external image file (works
#' only when an \code{HTML} document is exported and only for \code{JPEG} or
#' \code{PNG} images.
#' @param export.formats.html.template,export.formats.html.javascript,export.formats.html.css
#' \emph{character strings} or \code{NULL} (the default) giving the path to the
#' \code{HTML}, \code{JavaScript} and \code{CSS} templates, respectively, to be
#' used when generating the HTML+CSS semi-interactive plots; when \code{NULL},
#' the default ones included with the package will be used. If you decide to define
#' new templates please use the default ones for inspiration and note that future
#' version are not guaranteed to be backwards compatible!
#' @param export.formats.directory a \emph{string}; if exporting, which directory
#' to export to; if \code{NA} (the default), creates the files in a temporary
#' directory.
#' @param generate.R.plot a \emph{logical}, if \code{TRUE} (the default),
#' generate the standard (base \code{R}) plot for plotting within \code{R}.
#' @param do.not.draw.plot a \emph{logical}, if \code{TRUE} (\emph{not} the default),
#' does not draw the plot itself, but only the legend (if \code{show.legend} is
#' \code{TRUE}) at coordinates (0,0) irrespective of the given legend coordinates.
#' This is intended to allow (together with the \code{get.legend.plotting.area()}
#' function) the separate plotting of the legend.
#' @param ... other possible parameters
#' @examples
#' cma0 <- CMA0(data=med.events,
#'              ID.colname="PATIENT_ID",
#'              event.date.colname="DATE",
#'              event.duration.colname="DURATION",
#'              event.daily.dose.colname="PERDAY",
#'              medication.class.colname="CATEGORY",
#'              followup.window.start=0,
#'              followup.window.start.unit="days",
#'              followup.window.duration=2*365,
#'              followup.window.duration.unit="days",
#'              observation.window.start=30,
#'              observation.window.start.unit="days",
#'              observation.window.duration=365,
#'              observation.window.duration.unit="days",
#'              date.format="%m/%d/%Y",
#'              summary="Base CMA");
#' plot(cma0, patients.to.plot=c("1","2"));
#' @export
plot.CMA0 <- function(x,                                     # the CMA0 (or derived) object
                      ...,                                   # required for S3 consistency
                      patients.to.plot=NULL,                 # list of patient IDs to plot or NULL for all
                      duration=NA,                           # duration to plot in days (if missing, determined from the data)
                      align.all.patients=FALSE, align.first.event.at.zero=FALSE, # should all patients be aligned? and, if so, place the first event as the horizontal 0?
                      show.period=c("dates","days")[2],      # draw vertical bars at regular interval as dates or days?
                      period.in.days=90,                     # the interval (in days) at which to draw veritcal lines
                      show.legend=TRUE, legend.x="right", legend.y="bottom", legend.bkg.opacity=0.5, legend.cex=0.75, legend.cex.title=1.0, # legend params and position
                      cex=1.0, cex.axis=0.75, cex.lab=1.0,   # various graphical params
                      xlab=c("dates"="Date", "days"="Days"), # Vector of x labels to show for the two types of periods, or a single value for both, or NULL for nothing
                      ylab=c("withoutCMA"="patient", "withCMA"="patient (& CMA)"), # Vector of y labels to show without and with CMA estimates, or a single value for both, or NULL ofr nonthing
                      title=c("aligned"="Event patterns (all patients aligned)", "notaligned"="Event patterns"), # Vector of titles to show for and without alignment, or a single value for both, or NULL for nonthing
                      col.cats=rainbow,                      # single color or a function mapping the categories to colors
                      unspecified.category.label="drug",     # the label of the unspecified category of medication
                      medication.groups.to.plot=NULL,        # the names of the medication groups to plot (by default, all)
                      medication.groups.separator.show=TRUE, medication.groups.separator.lty="solid", medication.groups.separator.lwd=2, medication.groups.separator.color="blue", # group medication events by patient?
                      medication.groups.allother.label="*",  # the label to use for the __ALL_OTHERS__ medication class (defaults to *)
                      lty.event="solid", lwd.event=2, pch.start.event=15, pch.end.event=16, # event style
                      plot.events.vertically.displaced=TRUE, # display the events on different lines (vertical displacement) or not (defaults to TRUE)?
                      print.dose=FALSE, cex.dose=0.75, print.dose.outline.col="white", print.dose.centered=FALSE, # print daily dose
                      plot.dose=FALSE, lwd.event.max.dose=8, plot.dose.lwd.across.medication.classes=FALSE, # draw daily dose as line width
                      col.continuation="black", lty.continuation="dotted", lwd.continuation=1, # style of the contuniation lines connecting consecutive events
                      col.na="lightgray",                    # color for missing data
                      highlight.followup.window=TRUE, followup.window.col="green",
                      highlight.observation.window=TRUE, observation.window.col="yellow", observation.window.density=35, observation.window.angle=-30, observation.window.opacity=0.3,
                      alternating.bands.cols=c("white", "gray95"), # the colors of the alternating vertical bands across patients (NULL=don't draw any; can be >= 1 color)
                      rotate.text=-60,                       # some text (e.g., axis labels) may be rotated by this much degrees
                      force.draw.text=FALSE,                 # if true, always draw text even if too big or too small
                      bw.plot=FALSE,                         # if TRUE, override all user-given colors and replace them with a scheme suitable for grayscale plotting
                      min.plot.size.in.characters.horiz=0, min.plot.size.in.characters.vert=0, # the minimum plot size (in characters: horizontally, for the whole duration, vertically, per event)
                      suppress.warnings=FALSE,               # suppress warnings?
                      max.patients.to.plot=100,              # maximum number of patients to plot
                      export.formats=NULL,                   # the formats to export the figure to (by default, none); can be any subset of "svg" (just SVG file), "html" (SVG + HTML + CSS + JavaScript all embedded within the HTML document), "jpg", "png", "webp", "ps" and "pdf"
                      export.formats.fileprefix="AdhereR-plot", # the file name prefix for the exported formats
                      export.formats.height=NA, export.formats.width=NA, # desired dimensions (in pixels) for the exported figure (defaults to sane values)
                      export.formats.save.svg.placeholder=TRUE,
                      export.formats.svg.placeholder.type=c("jpg", "png", "webp")[2],
                      export.formats.svg.placeholder.embed=FALSE, # save a placeholder for the SVG image?
                      export.formats.html.template=NULL, export.formats.html.javascript=NULL, export.formats.html.css=NULL, # HTML, JavaScript and CSS templates for exporting HTML+SVG
                      export.formats.directory=NA,           # if exporting, which directory to export to (if not give, creates files in the temporary directory)
                      generate.R.plot=TRUE,                  # generate standard (base R) plot for plotting within R?
                      do.not.draw.plot=FALSE                 # if TRUE, don't draw the actual plot, but only the legend (if required)
)
{
  .plot.CMAs(x,
             patients.to.plot=patients.to.plot,
             duration=duration,
             align.all.patients=align.all.patients,
             align.first.event.at.zero=align.first.event.at.zero,
             show.period=show.period,
             period.in.days=period.in.days,
             show.legend=show.legend,
             legend.x=legend.x,
             legend.y=legend.y,
             legend.bkg.opacity=legend.bkg.opacity,
             legend.cex=legend.cex,
             legend.cex.title=legend.cex.title,
             cex=cex,
             cex.axis=cex.axis,
             cex.lab=cex.lab,
             show.cma=FALSE, # for CMA0, there's no CMA to show by definition
             xlab=xlab,
             ylab=ylab,
             title=title,
             col.cats=col.cats,
             unspecified.category.label=unspecified.category.label,
             medication.groups.to.plot=medication.groups.to.plot,
             medication.groups.separator.show=medication.groups.separator.show,
             medication.groups.separator.lty=medication.groups.separator.lty,
             medication.groups.separator.lwd=medication.groups.separator.lwd,
             medication.groups.separator.color=medication.groups.separator.color,
             medication.groups.allother.label=medication.groups.allother.label,
             lty.event=lty.event,
             lwd.event=lwd.event,
             show.event.intervals=FALSE, # not for CMA0
             plot.events.vertically.displaced=plot.events.vertically.displaced,
             pch.start.event=pch.start.event,
             pch.end.event=pch.end.event,
             print.dose=print.dose,
             cex.dose=cex.dose,
             print.dose.outline.col=print.dose.outline.col,
             print.dose.centered=print.dose.centered,
             plot.dose=plot.dose,
             lwd.event.max.dose=lwd.event.max.dose,
             plot.dose.lwd.across.medication.classes=plot.dose.lwd.across.medication.classes,
             col.na=col.na,
             col.continuation=col.continuation,
             lty.continuation=lty.continuation,
             lwd.continuation=lwd.continuation,
             print.CMA=FALSE, # for CMA0, there's no CMA to show by definition
             plot.CMA=FALSE, # for CMA0, there's no CMA to show by definition
             plot.partial.CMAs.as=NULL, # for CMA0, there's no "partial" CMAs by definition
             highlight.followup.window=highlight.followup.window,
             followup.window.col=followup.window.col,
             highlight.observation.window=highlight.observation.window,
             observation.window.col=observation.window.col,
             observation.window.density=observation.window.density,
             observation.window.angle=observation.window.angle,
             observation.window.opacity=observation.window.opacity,
             alternating.bands.cols=alternating.bands.cols,
             bw.plot=bw.plot,
             rotate.text=rotate.text,
             force.draw.text=force.draw.text,
             min.plot.size.in.characters.horiz=min.plot.size.in.characters.horiz,
             min.plot.size.in.characters.vert=min.plot.size.in.characters.vert,
             max.patients.to.plot=max.patients.to.plot,
             export.formats=export.formats,
             export.formats.fileprefix=export.formats.fileprefix,
             export.formats.height=export.formats.height,
             export.formats.width=export.formats.width,
             export.formats.save.svg.placeholder=export.formats.save.svg.placeholder,
             export.formats.svg.placeholder.type=export.formats.svg.placeholder.type,
             export.formats.svg.placeholder.embed=export.formats.svg.placeholder.embed,
             export.formats.html.template=export.formats.html.template,
             export.formats.html.javascript=export.formats.html.javascript,
             export.formats.html.css=export.formats.html.css,
             export.formats.directory=export.formats.directory,
             generate.R.plot=generate.R.plot,
             do.not.draw.plot=do.not.draw.plot,
             suppress.warnings=suppress.warnings);
}