#' CMA1 and CMA3 constructors.
#'
#' Constructs a CMA (continuous multiple-interval measures of medication
#' availability/gaps) type 1 or type 3 object.
#'
#' \code{CMA1} considers the total number of days with medication supplied in
#' all medication events in the observation window, excluding the last event.
#' \code{CMA3} is identical to \code{CMA1} except that it is capped at 100\%.
#'
#' The formula is
#' \deqn{(number of days supply excluding last) / (first to last event)}
#' Thus, the durations of all events are added up, possibly resulting in an CMA
#' estimate (much) bigger than 1.0 (100\%).
#'
#' \code{\link{CMA2}} and \code{CMA1} differ in the inclusion or not of the last
#' event.
#'
#' @param data A \emph{\code{data.frame}} containing the events used to compute
#' the CMA. Must contain, at a minimum, the patient unique ID, the event date
#' and duration, and might also contain the daily dosage and medication type
#' (the actual column names are defined in the following four parameters).
#' @param ID.colname A \emph{string}, the name of the column in \code{data}
#' containing the unique patient ID; must be present.
#' @param event.date.colname A \emph{string}, the name of the column in
#' \code{data} containing the start date of the event (in the format given in
#' the \code{date.format} parameter); must be present.
#' @param event.duration.colname A \emph{string}, the name of the column in
#' \code{data} containing the event duration (in days); must be present.
#' @param medication.groups A \emph{vector} of characters defining medication
#' groups or the name of a column in \code{data} that defines such groups.
#' The names of the vector are the medication group unique names, while
#' the content defines them as logical expressions. While the names can be any
#' string of characters except "\}", it is recommended to stick to the rules for
#' defining vector names in \code{R}. For example,
#' \code{c("A"="CATEGORY == 'medA'", "AA"="{A} & PERDAY < 4"} defines two
#' medication groups: \emph{A} which selects all events of type "medA", and
#' \emph{B} which selects all events already defined by "A" but with a daily
#' dose lower than 4. If \code{NULL}, no medication groups are defined. If
#' medication groups are defined, there is one CMA estimate for each group;
#' moreover, there is a special group \emph{__ALL_OTHERS__} automatically defined
#' containing all observations \emph{not} covered by any of the explicitly defined
#' groups.
#' @param flatten.medication.groups \emph{Logical}, if \code{FALSE} (the default)
#' then the \code{CMA} and \code{event.info} components of the object are lists
#' with one medication group per element; otherwise, they are \code{data.frame}s
#' with an extra column containing the medication group (its name is given by
#' \code{medication.groups.colname}).
#' @param medication.groups.colname a \emph{string} (defaults to ".MED_GROUP_ID")
#' giving the name of the column storing the group name when
#' \code{flatten.medication.groups} is \code{TRUE}.
#' @param followup.window.start If a \emph{\code{Date}} object, it represents
#' the actual start date of the follow-up window; if a \emph{string} it is the
#' name of the column in \code{data} containing the start date of the follow-up
#' window either as the numbers of \code{followup.window.start.unit} units after
#' the first event (the column must be of type \code{numeric}) or as actual
#' dates (in which case the column must be of type \code{Date} or a string
#' that conforms to the format specified in \code{date.format}); if a
#' \emph{number} it is the number of time units defined in the
#' \code{followup.window.start.unit} parameter after the begin of the
#' participant's first event; or \code{NA} if not defined.
#' @param followup.window.start.unit can be either \emph{"days"},
#' \emph{"weeks"}, \emph{"months"} or \emph{"years"}, and represents the time
#' units that \code{followup.window.start} refers to (when a number), or
#' \code{NA} if not defined.
#' @param followup.window.start.per.medication.group a \emph{logical}: if there are
#' medication groups defined and this is \code{TRUE}, then the first event
#' considered for the follow-up window start is relative to each medication group
#' separately, otherwise (the default) it is relative to the patient.
#' @param followup.window.duration either a \emph{number} representing the
#' duration of the follow-up window in the time units given in
#' \code{followup.window.duration.unit}, or a \emph{string} giving the column
#' containing these numbers. Should represent a period for which relevant
#' medication events are recorded accurately (e.g. not extend after end of
#' relevant treatment, loss-to-follow-up or change to a health care provider
#' not covered by the database).
#' @param followup.window.duration.unit can be either \emph{"days"},
#' \emph{"weeks"}, \emph{"months"} or \emph{"years"}, and represents the time
#' units that \code{followup.window.duration} refers to, or \code{NA} if not
#' defined.
#' @param observation.window.start,observation.window.start.unit,observation.window.duration,observation.window.duration.unit the definition of the observation window
#' (see the follow-up window parameters above for details).
#' @param date.format A \emph{string} giving the format of the dates used in the
#' \code{data} and the other parameters; see the \code{format} parameters of the
#' \code{\link[base]{as.Date}} function for details (NB, this concerns only the
#' dates given as strings and not as \code{Date} objects).
#' @param summary Metadata as a \emph{string}, briefly describing this CMA.
#' @param event.interval.colname A \emph{string}, the name of a newly-created
#' column storing the number of days between the start of the current event and
#' the start of the next one; the default value "event.interval" should be
#' changed only if there is a naming conflict with a pre-existing
#' "event.interval" column in \code{event.info}.
#' @param gap.days.colname A \emph{string}, the name of a newly-created column
#' storing the number of days when medication was not available (i.e., the
#' "gap days"); the default value "gap.days" should be changed only if there is
#' a naming conflict with a pre-existing "gap.days" column in \code{event.info}.
#' @param force.NA.CMA.for.failed.patients \emph{Logical} describing how the
#' patients for which the CMA estimation fails are treated: if \code{TRUE}
#' they are returned with an \code{NA} CMA estimate, while for
#' \code{FALSE} they are omitted.
#' @param parallel.backend Can be "none" (the default) for single-threaded
#' execution, "multicore"  (using \code{mclapply} in package \code{parallel})
#' for multicore processing (NB. not currently implemented on MS Windows and
#' automatically falls back on "snow" on this platform),  or "snow",
#' "snow(SOCK)" (equivalent to "snow"), "snow(MPI)" or "snow(NWS)" specifying
#' various types of SNOW clusters (can be on the local machine or more complex
#' setups -- please see the documentation of package \code{snow} for details;
#' the last two require packages \code{Rmpi} and \code{nws}, respectively, not
#' automatically installed with \code{AdhereR}).
#' @param parallel.threads Can be "auto" (for \code{parallel.backend} ==
#' "multicore", defaults to the number of cores in the system as given by
#' \code{options("cores")}, while for \code{parallel.backend} == "snow",
#' defaults to 2), a strictly positive integer specifying the number of parallel
#' threads, or a more complex specification of the SNOW cluster nodes for
#' \code{parallel.backend} == "snow" (see the documentation of package
#' \code{snow} for details).
#' @param suppress.warnings \emph{Logical}, if \code{TRUE} don't show any
#' warnings.
#' @param suppress.special.argument.checks \emph{Logical} parameter for internal
#' use; if \code{FALSE} (default) check if the important columns in the \code{data}
#' have some of the reserved names, if \code{TRUE} this check is not performed.
#' @param arguments.that.should.not.be.defined a \emph{list} of argument names
#' and pre-defined valuesfor which a warning should be thrown if passed to the
#' function.
#' @param ... other possible parameters
#' @return An \code{S3} object of class \code{CMA1} (derived from \code{CMA0})
#' with the following fields:
#' \itemize{
#'  \item \code{data} The actual event data, as given by the \code{data}
#'  parameter.
#'  \item \code{ID.colname} the name of the column in \code{data} containing the
#'  unique patient ID, as given by the \code{ID.colname} parameter.
#'  \item \code{event.date.colname} the name of the column in \code{data}
#'  containing the start date of the event (in the format given in the
#'  \code{date.format} parameter), as given by the \code{event.date.colname}
#'  parameter.
#'  \item \code{event.duration.colname} the name of the column in \code{data}
#'  containing the event duration (in days), as given by the
#'  \code{event.duration.colname} parameter.
#'  \item \code{event.daily.dose.colname} the name of the column in \code{data}
#'  containing the prescribed daily dose, as given by the
#'  \code{event.daily.dose.colname} parameter.
#'  \item \code{medication.class.colname} the name of the column in \code{data}
#'  containing the classes/types/groups of medication, as given by the
#'  \code{medication.class.colname} parameter.
#'  \item \code{followup.window.start} the beginning of the follow-up window, as
#'  given by the \code{followup.window.start} parameter.
#'  \item \code{followup.window.start.unit} the time unit of the
#'  \code{followup.window.start}, as given by the
#'  \code{followup.window.start.unit} parameter.
#'  \item \code{followup.window.duration} the duration of the follow-up window,
#'  as given by the \code{followup.window.duration} parameter.
#'  \item \code{followup.window.duration.unit} the time unit of the
#'  \code{followup.window.duration}, as given by the
#'  \code{followup.window.duration.unit} parameter.
#'  \item \code{observation.window.start} the beginning of the observation
#'  window, as given by the \code{observation.window.start} parameter.
#'  \item \code{observation.window.start.unit} the time unit of the
#'  \code{observation.window.start}, as given by the
#'  \code{observation.window.start.unit} parameter.
#'  \item \code{observation.window.duration} the duration of the observation
#'  window, as given by the \code{observation.window.duration} parameter.
#'  \item \code{observation.window.duration.unit} the time unit of the
#'  \code{observation.window.duration}, as given by the
#'  \code{observation.window.duration.unit} parameter.
#'  \item \code{date.format} the format of the dates, as given by the
#'  \code{date.format} parameter.
#'  \item \code{summary} the metadata, as given by the \code{summary}
#'  parameter.
#'  \item \code{event.info} the \code{data.frame} containing the event info
#'  (irrelevant for most users; see \code{\link{compute.event.int.gaps}} for
#'  details).
#'  \item \code{CMA} the \code{data.frame} containing the actual \code{CMA}
#'  estimates for each participant (the \code{ID.colname} column).
#' }
#' Please note that if \code{medication.groups} are defined and
#' \code{flatten.medication.groups} is \code{FALSE}, then the \code{CMA}
#' and \code{event.info} are named lists, each element containing the CMA and
#' event.info corresponding to a single medication group (the element's name),
#' but if \code{flatten.medication.groups} is \code{FALSE} then they are
#' \code{data.frame}s with an extra column giving the medication group (the
#' column's name is given by \code{medication.groups.colname}).
#' @seealso CMAs 1 to 8 are described in:
#'
#' Vollmer, W. M., Xu, M., Feldstein, A., Smith, D., Waterbury, A., & Rand, C.
#' (2012). Comparison of pharmacy-based measures of medication adherence.
#' \emph{BMC Health Services Research}, \strong{12}, 155.
#' \doi{10.1186/1472-6963-12-155}.
#'
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
#' cma3 <- CMA3(data=med.events,
#'              ID.colname="PATIENT_ID",
#'              event.date.colname="DATE",
#'              event.duration.colname="DURATION",
#'              followup.window.start=30,
#'              observation.window.start=30,
#'              observation.window.duration=365,
#'              date.format="%m/%d/%Y"
#'             );
#' @export
CMA1 <- function( data=NULL, # the data used to compute the CMA on
                  # Important columns in the data
                  ID.colname=NA, # the name of the column containing the unique patient ID (NA = undefined)
                  event.date.colname=NA, # the start date of the event in the date.format format (NA = undefined)
                  event.duration.colname=NA, # the event duration in days (NA = undefined)
                  # Groups of medication classes:
                  medication.groups=NULL, # a named vector of medication group definitions, the name of a column in the data that defines the groups, or NULL
                  flatten.medication.groups=FALSE, medication.groups.colname=".MED_GROUP_ID", # if medication.groups were defined, return CMAs and event.info as single data.frame?
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
                  summary=NA,
                  # The description of the output (added) columns:
                  event.interval.colname="event.interval", # contains number of days between the start of current event and the start of the next
                  gap.days.colname="gap.days", # contains the number of days when medication was not available
                  # Dealing with failed estimates:
                  force.NA.CMA.for.failed.patients=TRUE, # force the failed patients to have NA CMA estimates?
                  # Parallel processing:
                  parallel.backend=c("none","multicore","snow","snow(SOCK)","snow(MPI)","snow(NWS)")[1], # parallel backend to use
                  parallel.threads="auto", # specification (or number) of parallel threads
                  # Misc:
                  suppress.warnings=FALSE,
                  suppress.special.argument.checks=TRUE, # used internally to suppress the check that we don't use special argument names
                  arguments.that.should.not.be.defined=c("carryover.within.obs.window"=FALSE,
                                                         "carryover.into.obs.window"=FALSE,
                                                         "carry.only.for.same.medication"=FALSE,
                                                         "consider.dosage.change"=FALSE), # the list of argument names and values for which a warning should be thrown if passed to the function
                  ...
                )
{
  # The summary:
  if( is.na(summary) ) summary <- "The ratio of days with medication available in the observation window excluding the last event; durations of all events added up and divided by number of days from first to last event, possibly resulting in a value >1.0";

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
        .report.ewms(paste0("Please note that '",args.list[[1]],"' overrides argument '",names(arguments.that.should.not.be.defined)[i],"' with value '",arguments.that.should.not.be.defined[i],"'!\n"), "warning", "CMA1", "AdhereR");
      }
    }
  }

  # Create the CMA0 object:
  ret.val <- CMA0(data=data,
                  ID.colname=ID.colname,
                  event.date.colname=event.date.colname,
                  event.duration.colname=event.duration.colname,
                  medication.groups=medication.groups,
                  flatten.medication.groups=flatten.medication.groups,
                  medication.groups.colname=medication.groups.colname,
                  followup.window.start=followup.window.start,
                  followup.window.start.unit=followup.window.start.unit,
                  followup.window.start.per.medication.group=followup.window.start.per.medication.group,
                  followup.window.duration=followup.window.duration,
                  followup.window.duration.unit=followup.window.duration.unit,
                  observation.window.start=observation.window.start,
                  observation.window.start.unit=observation.window.start.unit,
                  observation.window.duration=observation.window.duration,
                  observation.window.duration.unit=observation.window.duration.unit,
                  date.format=date.format,
                  summary=summary,
                  suppress.warnings=suppress.warnings,
                  suppress.special.argument.checks=suppress.special.argument.checks);
  if( is.null(ret.val) ) return (NULL); # some serious error upstream
  # The followup.window.start and observation.window.start might have been converted to Date:
  followup.window.start <- ret.val$followup.window.start; observation.window.start <- ret.val$observation.window.start;

  # The workhorse auxiliary function:
  .workhorse.function <- function(data=NULL,
                                  ID.colname=NULL,
                                  event.date.colname=NULL,
                                  event.duration.colname=NULL,
                                  event.daily.dose.colname=NULL,
                                  medication.class.colname=NULL,
                                  event.interval.colname=NULL,
                                  gap.days.colname=NULL,
                                  carryover.within.obs.window=NULL,
                                  carryover.into.obs.window=NULL,
                                  carry.only.for.same.medication=NULL,
                                  consider.dosage.change=NULL,
                                  followup.window.start=NULL,
                                  followup.window.start.unit=NULL,
                                  followup.window.duration=NULL,
                                  followup.window.duration.unit=NULL,
                                  observation.window.start=NULL,
                                  observation.window.start.unit=NULL,
                                  observation.window.duration=NULL,
                                  observation.window.duration.unit=NULL,
                                  date.format=NULL,
                                  suppress.warnings=NULL,
                                  suppress.special.argument.checks=NULL
  )
  {
    # Auxiliary internal function: Compute the CMA for a given patient:
    .process.patient <- function(data4ID)
    {
      # Force the selection, evaluation of promises and caching of the needed columns:
      # ... which columns to select (with their indices):
      columns.to.cache <- c(".EVENT.STARTS.BEFORE.OBS.WINDOW", ".EVENT.STARTS.AFTER.OBS.WINDOW", ".DATE.as.Date", event.duration.colname);
      # ... select these columns:
      data4ID.selected.columns <- data4ID[, columns.to.cache, with=FALSE]; # alternative to: data4ID[,..columns.to.cache];
      # ... cache the columns based on their indices:
      .EVENT.STARTS.BEFORE.OBS.WINDOW <- data4ID.selected.columns[[1]];
      .EVENT.STARTS.AFTER.OBS.WINDOW <- data4ID.selected.columns[[2]];
      .DATE.as.Date <- data4ID.selected.columns[[3]];
      event.duration.column <- data4ID.selected.columns[[4]];

      # which data to consider:
      s <- which(!(.EVENT.STARTS.BEFORE.OBS.WINDOW | .EVENT.STARTS.AFTER.OBS.WINDOW));
      s.len <- length(s); s1 <- s[1]; ss.len <- s[s.len];

      if( s.len < 2 || (.date.diff <- .difftime.Dates.as.days(.DATE.as.Date[ss.len], .DATE.as.Date[s1])) == 0 )
      {
        # For less than two events or when the first and the last events are on the same day, CMA1 does not make sense
        return (list("CMA"=NA_real_,
                     "included.events.original.row.order"=NA_integer_));
      } else
      {
        # Otherwise, the sum of durations of the events excluding the last divided by the number of days between the first and the last event
        return (list("CMA"=as.numeric(sum(event.duration.column[s[-s.len]],na.rm=TRUE) / .date.diff),
                     "included.events.original.row.order"=data4ID$..ORIGINAL.ROW.ORDER..[s[-s.len]]));
      }
    }

    # Call the compute.event.int.gaps() function and use the results:
    event.info <- compute.event.int.gaps(data=as.data.frame(data),
                                         ID.colname=ID.colname,
                                         event.date.colname=event.date.colname,
                                         event.duration.colname=event.duration.colname,
                                         event.daily.dose.colname=event.daily.dose.colname,
                                         medication.class.colname=medication.class.colname,
                                         event.interval.colname=event.interval.colname,
                                         gap.days.colname=gap.days.colname,
                                         carryover.within.obs.window=carryover.within.obs.window,
                                         carryover.into.obs.window=carryover.into.obs.window,
                                         carry.only.for.same.medication=carry.only.for.same.medication,
                                         consider.dosage.change=consider.dosage.change,
                                         followup.window.start=followup.window.start,
                                         followup.window.start.unit=followup.window.start.unit,
                                         followup.window.duration=followup.window.duration,
                                         followup.window.duration.unit=followup.window.duration.unit,
                                         observation.window.start=observation.window.start,
                                         observation.window.start.unit=observation.window.start.unit,
                                         observation.window.duration=observation.window.duration,
                                         observation.window.duration.unit=observation.window.duration.unit,
                                         date.format=date.format,
                                         keep.window.start.end.dates=TRUE,
                                         parallel.backend="none", # make sure this runs sequentially!
                                         parallel.threads=1,
                                         suppress.warnings=suppress.warnings,
                                         suppress.special.argument.checks=suppress.special.argument.checks,
                                         return.data.table=TRUE);
    if( is.null(event.info) ) return (list("CMA"=NA, "event.info"=NULL));

    CMA <- event.info[, .process.patient(.SD), by=ID.colname];
    event.info$.EVENT.USED.IN.CMA <- (event.info$..ORIGINAL.ROW.ORDER.. %in% CMA$included.events.original.row.order); # store if the event was used to compute the CMA or not
    return (list("CMA"=unique(CMA[,1:2]), "event.info"=event.info));
  }

  ret.val <- .cma.skeleton(data=data,
                           ret.val=ret.val,
                           cma.class.name="CMA1",

                           ID.colname=ID.colname,
                           event.date.colname=event.date.colname,
                           event.duration.colname=event.duration.colname,
                           event.daily.dose.colname=NA, # not relevant
                           medication.class.colname=NA, # not relevant
                           event.interval.colname=event.interval.colname,
                           gap.days.colname=gap.days.colname,
                           carryover.within.obs.window=FALSE, # if TRUE consider the carry-over within the observation window
                           carryover.into.obs.window=FALSE, # if TRUE consider the carry-over from before the starting date of the observation window
                           carry.only.for.same.medication=FALSE, # if TRUE the carry-over applies only across medication of same type
                           consider.dosage.change=FALSE, # if TRUE carry-over is adjusted to reflect changes in dosage
                           followup.window.start=followup.window.start,
                           followup.window.start.unit=followup.window.start.unit,
                           followup.window.duration=followup.window.duration,
                           followup.window.duration.unit=followup.window.duration.unit,
                           observation.window.start=observation.window.start,
                           observation.window.start.unit=observation.window.start.unit,
                           observation.window.duration=observation.window.duration,
                           observation.window.duration.unit=observation.window.duration.unit,
                           date.format=date.format,

                           flatten.medication.groups=flatten.medication.groups,
                           followup.window.start.per.medication.group=followup.window.start.per.medication.group,

                           suppress.warnings=suppress.warnings,
                           suppress.special.argument.checks=suppress.special.argument.checks,
                           force.NA.CMA.for.failed.patients=force.NA.CMA.for.failed.patients,
                           parallel.backend=parallel.backend,
                           parallel.threads=parallel.threads,
                           .workhorse.function=.workhorse.function);

  return (ret.val);
}


#' @rdname print.CMA0
#' @export
print.CMA1 <- function(...) print.CMA0(...)

#' Plot CMA0-derived objects.
#'
#' Plots the event data and estimated CMA encapsulated in objects derived from
#' \code{CMA0}.
#'
#' Please note that this function plots objects inheriting from \code{CMA0} but
#' not objects of type \code{CMA0} itself (these are plotted by
#' \code{\link{plot.CMA0}}).
#'
#' The x-axis represents time (either in days since the earliest date or as
#' actual dates), with consecutive events represented as ascending on the y-axis.
#'
#' Each event is represented as a segment with style \code{lty.event} and line
#' width \code{lwd.event} starting with a \code{pch.start.event} and ending with
#' a \code{pch.end.event} character, coloured with a unique color as given by
#' \code{col.cats}, extending from its start date until its end date.
#' Superimposed on these are shown the event intervals and gap days as estimated
#' by the particular CMA method, more precisely plotting the start and end of
#' the available events as solid filled-in rectangles, and the event gaps as
#' shaded rectangles.
#'
#' The follow-up and the observation windows are plotted as an empty rectangle
#' and as shaded rectangle, respectively (for some CMAs the observation window
#' might be adjusted in which case the adjustment may also be plotted using a
#'  different shading).
#'
#' The CMA estimates can be visually represented as well in the left side of the
#' figure using bars (sometimes the estimates can go above 100\%, in which case
#' the maximum possible bar filling is adjusted to reflect this).
#'
#' When several patients are displayed on the same plot, they are organized
#' vertically, and alternating bands (white and gray) help distinguish
#' consecutive patients.
#' Implicitely, all patients contained in the \code{cma} object will be plotted,
#' but the \code{patients.to.plot} parameter allows the selection of a subset of
#' patients.
#'
#' Finally, the y-axis shows the patient ID and possibly the CMA estimate as
#' well.
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
#' @param cex,cex.axis,cex.lab,legend.cex,legend.cex.title,CMA.cex \emph{numeric}
#' values specifying the \code{cex} of the various types of text.
#' @param show.cma \emph{Logical}, should the CMA type be shown in the title?
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
#' @param show.event.intervals \emph{Logical}, should the actual event intervals
#' be shown?
#' @param col.na The colour used for missing event data.
#' @param bw.plot \emph{Logical}, should the plot use grayscale only (i.e., the
#' \code{\link[grDevices]{gray.colors}} function)?
#' @param rotate.text \emph{Numeric}, the angle by which certain text elements
#' (e.g., axis labels) should be rotated.
#' @param force.draw.text \emph{Logical}, if \code{TRUE}, always draw text even
#' if too big or too small
#' @param print.CMA \emph{Logical}, should the CMA values be printed?
#' @param plot.CMA \emph{Logical}, should the CMA values be represented
#' graphically?
#' @param CMA.plot.ratio A \emph{number}, the proportion of the total horizontal
#' plot space to be allocated to the CMA plot.
#' @param CMA.plot.col,CMA.plot.border,CMA.plot.bkg,CMA.plot.text \emph{Strings}
#' giving the colours of the various components of the CMA plot.
#' @param highlight.followup.window \emph{Logical}, should the follow-up window
#' be plotted?
#' @param followup.window.col The follow-up window's colour.
#' @param highlight.observation.window \emph{Logical}, should the observation
#' window be plotted?
#' @param observation.window.col,observation.window.density,observation.window.angle,observation.window.opacity
#' Attributes of the observation window (colour, shading density, angle and
#' opacity).
#' @param show.real.obs.window.start,real.obs.window.density,real.obs.window.angle For some CMAs, the observation window might
#' be adjusted, in which case should it be plotted and with that attributes?
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
#' @param alternating.bands.cols The colors of the alternating vertical bands
#' distinguishing the patients; can be \code{NULL} = don't draw the bandes;
#' or a vector of colors.
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
#' cma1 <- CMA1(data=med.events,
#'              ID.colname="PATIENT_ID",
#'              event.date.colname="DATE",
#'              event.duration.colname="DURATION",
#'              followup.window.start=30,
#'              observation.window.start=30,
#'              observation.window.duration=365,
#'              date.format="%m/%d/%Y"
#'             );
#' plot(cma1, patients.to.plot=c("1","2"));
#' @export
plot.CMA1 <- function(x,                                     # the CMA1 (or derived) object
                      ...,                                   # required for S3 consistency
                      patients.to.plot=NULL,                 # list of patient IDs to plot or NULL for all
                      duration=NA,                           # duration to plot in days (if missing, determined from the data)
                      align.all.patients=FALSE, align.first.event.at.zero=FALSE, # should all patients be aligned? and, if so, place the first event as the horizintal 0?
                      show.period=c("dates","days")[2],      # draw vertical bars at regular interval as dates or days?
                      period.in.days=90,                     # the interval (in days) at which to draw veritcal lines
                      show.legend=TRUE, legend.x="right", legend.y="bottom", legend.bkg.opacity=0.5, legend.cex=0.75, legend.cex.title=1.0, # legend params and position
                      cex=1.0, cex.axis=0.75, cex.lab=1.0,   # various graphical params
                      show.cma=TRUE,                         # show the CMA type
                      col.cats=rainbow,                      # single color or a function mapping the categories to colors
                      unspecified.category.label="drug",     # the label of the unspecified category of medication
                      medication.groups.to.plot=NULL,        # the names of the medication groups to plot (by default, all)
                      medication.groups.separator.show=TRUE, medication.groups.separator.lty="solid", medication.groups.separator.lwd=2, medication.groups.separator.color="blue", # group medication events by patient?
                      medication.groups.allother.label="*",  # the label to use for the __ALL_OTHERS__ medication class (defaults to *)
                      lty.event="solid", lwd.event=2, pch.start.event=15, pch.end.event=16, # event style
                      show.event.intervals=TRUE,             # show the actual rpescription intervals
                      col.na="lightgray",                    # color for mising data
                      #col.continuation="black", lty.continuation="dotted", lwd.continuation=1, # style of the contuniation lines connecting consecutive events
                      print.CMA=TRUE, CMA.cex=0.50,           # print CMA next to the participant's ID?
                      plot.CMA=TRUE,                   # plot the CMA next to the participant ID?
                      CMA.plot.ratio=0.10,             # the proportion of the total horizontal plot to be taken by the CMA plot
                      CMA.plot.col="lightgreen", CMA.plot.border="darkgreen", CMA.plot.bkg="aquamarine", CMA.plot.text=CMA.plot.border, # attributes of the CMA plot
                      highlight.followup.window=TRUE, followup.window.col="green",
                      highlight.observation.window=TRUE, observation.window.col="yellow", observation.window.density=35, observation.window.angle=-30, observation.window.opacity=0.3,
                      show.real.obs.window.start=TRUE, real.obs.window.density=35, real.obs.window.angle=30, # for some CMAs, the real observation window starts at a different date
                      print.dose=FALSE, cex.dose=0.75, print.dose.outline.col="white", print.dose.centered=FALSE, # print daily dose
                      plot.dose=FALSE, lwd.event.max.dose=8, plot.dose.lwd.across.medication.classes=FALSE, # draw daily dose as line width
                      alternating.bands.cols=c("white", "gray95"), # the colors of the alternating vertical bands across patients (NULL=don't draw any; can be >= 1 color)
                      bw.plot=FALSE,                         # if TRUE, override all user-given colors and replace them with a scheme suitable for grayscale plotting
                      rotate.text=-60,                       # some text (e.g., axis labels) may be rotated by this much degrees
                      force.draw.text=FALSE,                 # if true, always draw text even if too big or too small
                      min.plot.size.in.characters.horiz=0, min.plot.size.in.characters.vert=0, # the minimum plot size (in characters: horizontally, for the whole duration, vertically, per event)
                      max.patients.to.plot=100,              # maximum number of patients to plot
                      export.formats=NULL,                   # the formats to export the figure to (by default, none); can be any subset of "svg" (just SVG file), "html" (SVG + HTML + CSS + JavaScript all embedded within the HTML document), "jpg", "png", "webp", "ps" and "pdf"
                      export.formats.fileprefix="AdhereR-plot", # the file name prefix for the exported formats
                      export.formats.height=NA, export.formats.width=NA, # desired dimensions (in pixels) for the exported figure (defaults to sane values)
                      export.formats.save.svg.placeholder=TRUE,
                      export.formats.svg.placeholder.type=c("jpg", "png", "webp")[2],
                      export.formats.svg.placeholder.embed=FALSE, # save a placeholder for the SVG image?
                      export.formats.directory=NA,           # if exporting, which directory to export to (if not give, creates files in the temporary directory)
                      export.formats.html.template=NULL, export.formats.html.javascript=NULL, export.formats.html.css=NULL, # HTML, JavaScript and CSS templates for exporting HTML+SVG
                      generate.R.plot=TRUE,                  # generate standard (base R) plot for plotting within R?
                      do.not.draw.plot=FALSE                 # if TRUE, don't draw the actual plot, but only the legend (if required)
)
{
  .plot.CMA1plus(cma=x,
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
                 show.cma=show.cma,
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
                 pch.start.event=pch.start.event,
                 pch.end.event=pch.end.event,
                 show.event.intervals=show.event.intervals,
                 col.na=col.na,
                 print.CMA=print.CMA,
                 CMA.cex=CMA.cex,
                 plot.CMA=plot.CMA,
                 CMA.plot.ratio=CMA.plot.ratio,
                 CMA.plot.col=CMA.plot.col,
                 CMA.plot.border=CMA.plot.border,
                 CMA.plot.bkg=CMA.plot.bkg,
                 CMA.plot.text=CMA.plot.text,
                 highlight.followup.window=highlight.followup.window,
                 followup.window.col=followup.window.col,
                 highlight.observation.window=highlight.observation.window,
                 observation.window.col=observation.window.col,
                 observation.window.density=observation.window.density,
                 observation.window.angle=observation.window.angle,
                 observation.window.opacity=observation.window.opacity,
                 show.real.obs.window.start=show.real.obs.window.start,
                 real.obs.window.density=real.obs.window.density,
                 real.obs.window.angle=real.obs.window.angle,
                 print.dose=print.dose,
                 cex.dose=cex.dose,
                 print.dose.outline.col=print.dose.outline.col,
                 print.dose.centered=print.dose.centered,
                 plot.dose=plot.dose,
                 lwd.event.max.dose=lwd.event.max.dose,
                 plot.dose.lwd.across.medication.classes=plot.dose.lwd.across.medication.classes,
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
                 ...)
}
