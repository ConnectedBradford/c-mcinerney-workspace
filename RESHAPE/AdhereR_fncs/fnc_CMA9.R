#' CMA9 constructor.
#'
#' Constructs a CMA (continuous multiple-interval measures of medication
#' availability/gaps) type 9 object.
#'
#'
#' \code{CMA9} is similar to \code{CMA7} and \code{CMA8} in that it accounts for
#' carry-over within and before the observation window assuming that new
#' medication is "banked" until needed (oversupply from previous events is used
#' first, followed new medication supply). Yet, unlike these previous CMAs, it
#' does not assume the medication is used as prescribed; in longitudinal studies
#' with multiple CMA measures, this assumption may introduce additional
#' variation in CMA estimates depending on when the observation window starts
#' in relation to the previous medication event. A shorter time distance from
#' the previous event (and longer to the first event in the observation window)
#' results in higher values even if the number of gap days is the same, and it
#' may also be that the patient has had a similar use pattern for that time
#' interval, rather than perfect adherence followed by no medication use.
#' \code{CMA9} applies a different adjustment: it computes a ratio of days
#' supply over each interval between two prescriptions and considers this
#' applies for each day of that interval, up to 100\% (moving oversupply to the
#' next event interval). All medication events in the follow up window before
#' observation window are considered for carry-over calculation. The last
#' interval ends at the end of the follow-up window.
#' Thus, it accounts for timing within  the observation window, as well as
#' before (but differently from \code{CMA7} and \code{CMA8}), and excludes the
#' remaining supply at the end of the observation window, if any.
#'
#' The formula is
#' \deqn{(number of days in the observation window, each weighted by the ratio
#' of days supply applicable to their event interval) / (start to end of
#' observation window)}
#'
#' Observations:
#' \itemize{
#'  \item the \code{carry.only.for.same.medication} parameter controls the
#'  transmission of carry-over across medication changes, producing a "standard"
#'  \code{CMA7} (default value is FALSE), and an "alternative" \code{CMA7b},
#'  respectively;
#'  \item the \code{consider.dosage.change} parameter controls if dosage changes
#'  are taken into account, i.e. if set as TRUE and a new medication event has a
#'  different daily dosage recommendation, carry-over is recomputed assuming
#'  medication use according to the new prescribed dosage (default value is
#'  FALSE).
#' }
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
#' @param event.daily.dose.colname A \emph{string}, the name of the column in
#' \code{data} containing the prescribed daily dose, or \code{NA} if not defined.
#' @param medication.class.colname A \emph{string}, the name of the column in
#' \code{data} containing the medication type, or \code{NA} if not defined.
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
#' @param carry.only.for.same.medication \emph{Logical}, if \code{TRUE}, the
#' carry-over applies only across medication of the same type.
#' @param consider.dosage.change \emph{Logical}, if \code{TRUE}, the carry-over
#' is adjusted to also reflect changes in dosage.
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
#' relevant treatment, loss-to-follow-up or change to a health care provider not
#' covered by the database).
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
#' @return An \code{S3} object of class \code{CMA9} (derived from \code{CMA0})
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
#'  \item \code{carry.only.for.same.medication} whether the carry-over applies
#'  only across medication of the same type, as given by the
#'  \code{carry.only.for.same.medication} parameter.
#'  \item \code{consider.dosage.change} whether the carry-over is adjusted to
#'  reflect changes in dosage, as given by the \code{consider.dosage.change}
#'  parameter.
#'  \item \code{followup.window.start} the beginning of the follow-up window,
#'  as given by the \code{followup.window.start} parameter.
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
#'  \item \code{summary} the metadata, as given by the \code{summary} parameter.
#'  \item \code{event.info} the \code{data.frame} containing the event info
#'  (irrelevant for most users; see \code{\link{compute.event.int.gaps}} for
#'  details).
#'  \item \code{CMA} the \code{data.frame} containing the actual \code{CMA}
#'  estimates for each participant (the \code{ID.colname} column).
#' }
#' Please note that if \code{medication.groups} are defined, then the \code{CMA}
#' and \code{event.info} are named lists, each element containing the CMA and
#' event.info corresponding to a single medication group (the element's name).
#' @examples
#' cma9 <- CMA9(data=med.events,
#'              ID.colname="PATIENT_ID",
#'              event.date.colname="DATE",
#'              event.duration.colname="DURATION",
#'              event.daily.dose.colname="PERDAY",
#'              medication.class.colname="CATEGORY",
#'              carry.only.for.same.medication=FALSE,
#'              consider.dosage.change=FALSE,
#'              followup.window.start=30,
#'              observation.window.start=30,
#'              observation.window.duration=365,
#'              date.format="%m/%d/%Y"
#'             );
#' @export
CMA9 <- function( data=NULL, # the data used to compute the CMA on
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
                  carry.only.for.same.medication=FALSE, # if TRUE the carry-over applies only across medication of same type (NA = undefined)
                  consider.dosage.change=FALSE, # if TRUE carry-over is adjusted to reflect changes in dosage (NA = undefined)
                  # The follow-up window:
                  followup.window.start=0, # if a number is the earliest event per participant date + number of units, or a Date object, or a column name in data (NA = undefined)
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
                  force.NA.CMA.for.failed.patients=TRUE, # force the failed patients to have NA CM estimate?
                  # Parallel processing:
                  parallel.backend=c("none","multicore","snow","snow(SOCK)","snow(MPI)","snow(NWS)")[1], # parallel backend to use
                  parallel.threads="auto", # specification (or number) of parallel threads
                  # Misc:
                  suppress.warnings=FALSE,
                  suppress.special.argument.checks=TRUE, # used internally to suppress the check that we don't use special argument names
                  arguments.that.should.not.be.defined=c("carryover.within.obs.window"=TRUE,
                                                         "carryover.into.obs.window"=TRUE), # the list of argument names and values for which a warning should be thrown if passed to the function
                  ...
                )
{
  # The summary:
  if( is.na(summary) ) summary <- "The ratio of days with medication available in the observation window; the supply of each event is evenly spread until the next event (ratio days supply up to 1), then oversupply carried over to the next event; the each day in the observation window is weighted by its ratio of days supply, and the sum divided by the duration of the observation window";

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
        .report.ewms(paste0("Please note that '",args.list[[1]],"' overrides argument '",names(arguments.that.should.not.be.defined)[i],"' with value '",arguments.that.should.not.be.defined[i],"'!\n"), "warning", "CMA9", "AdhereR");
      }
    }
  }

  # Create the CMA0 object:
  ret.val <- CMA0(data=data,
                  ID.colname=ID.colname,
                  event.date.colname=event.date.colname,
                  event.duration.colname=event.duration.colname,
                  event.daily.dose.colname=event.daily.dose.colname,
                  medication.class.colname=medication.class.colname,
                  carry.only.for.same.medication=carry.only.for.same.medication,
                  consider.dosage.change=consider.dosage.change,
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
  if( is.null(ret.val) ) return (NULL); # some error upstream
  # The followup.window.start and observation.window.start might have been converted to Date:
  followup.window.start <- ret.val$followup.window.start; observation.window.start <- ret.val$observation.window.start;

  # The workhorse auxiliary function: For a given (subset) of data, compute the event intervals and gaps:
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
      n.events <- nrow(data4ID); # cache number of events
      # Caching;
      .obs.start.actual <- data4ID$.OBS.START.DATE.ACTUAL[1]; .obs.end.actual <- data4ID$.OBS.END.DATE.ACTUAL[1];
      # Select the events within the observation window:
      s <- which(data4ID$.DATE.as.Date <= .obs.end.actual);
      if( length(s) == 0 )
      {
        return (list("CMA"=NA_real_,
                     "included.events.original.row.order"=NA_integer_));
      } else
      {
        # Compute the ratios of CMA within the FUW:
        .CMA.PER.INTERVAL <- (1 - data4ID[,get(gap.days.colname)] / data4ID[,get(event.interval.colname)]);
        event.start <- data4ID$.DATE.as.Date; event.end <- (event.start + data4ID[,get(event.interval.colname)]);
        int.start <- pmax(event.start, .obs.start.actual); int.end <- pmin(event.end, .obs.end.actual);
        .EVENT.INTERVAL.IN.OBS.WIN <- pmax(0, as.numeric(int.end - int.start));
        # Weight the days by the ratio:
        return (list("CMA"=sum(.CMA.PER.INTERVAL * .EVENT.INTERVAL.IN.OBS.WIN, na.rm=TRUE) /
                       as.numeric(.obs.end.actual - .obs.start.actual),
                     "included.events.original.row.order"=data4ID$..ORIGINAL.ROW.ORDER..[s]));
      }
    }

    # ATTENTION: this is done for an observation window that is identical to the follow-up window to accommodate events that overshoot the observation window!
    event.info <- compute.event.int.gaps(data=as.data.frame(data),
                                         ID.colname=ID.colname,
                                         event.date.colname=event.date.colname,
                                         event.duration.colname=event.duration.colname,
                                         event.daily.dose.colname=event.daily.dose.colname,
                                         medication.class.colname=medication.class.colname,
                                         event.interval.colname=event.interval.colname,
                                         gap.days.colname=gap.days.colname,
                                         carryover.within.obs.window=TRUE,
                                         carryover.into.obs.window=TRUE,
                                         carry.only.for.same.medication=carry.only.for.same.medication,
                                         consider.dosage.change=consider.dosage.change,
                                         followup.window.start=followup.window.start,
                                         followup.window.start.unit=followup.window.start.unit,
                                         followup.window.duration=followup.window.duration,
                                         followup.window.duration.unit=followup.window.duration.unit,
                                         observation.window.start=0, # force follow-up window start at 0!
                                         observation.window.start.unit="days",
                                         observation.window.duration=followup.window.duration,
                                         observation.window.duration.unit=followup.window.duration.unit,
                                         date.format=date.format,
                                         keep.window.start.end.dates=TRUE,
                                         keep.event.interval.for.all.events=TRUE,
                                         parallel.backend="none", # make sure this runs sequentially!
                                         parallel.threads=1,
                                         suppress.warnings=suppress.warnings,
                                         suppress.special.argument.checks=suppress.special.argument.checks,
                                         return.data.table=TRUE);
    if( is.null(event.info) ) return (list("CMA"=NA, "event.info"=NULL));

    # Also compute the actual observation window start and end:
    #patinfo.cols <- which(names(event.info) %in% c(ID.colname, event.date.colname, event.duration.colname, event.daily.dose.colname, medication.class.colname));
    #tmp <- as.data.frame(data); tmp <- tmp[!duplicated(tmp[,ID.colname]),patinfo.cols]; # the reduced dataset for comuting the actual OW:
    tmp <- as.data.frame(data); tmp <- tmp[!duplicated(tmp[,ID.colname]),names(data)]; # the reduced dataset for comuting the actual OW:
    actual.obs.win <- compute.event.int.gaps(data=tmp,
                                             ID.colname=ID.colname,
                                             event.date.colname=event.date.colname,
                                             event.duration.colname=event.duration.colname,
                                             event.daily.dose.colname=event.daily.dose.colname,
                                             medication.class.colname=medication.class.colname,
                                             event.interval.colname=event.interval.colname,
                                             gap.days.colname=gap.days.colname,
                                             carryover.within.obs.window=FALSE,
                                             carryover.into.obs.window=FALSE,
                                             carry.only.for.same.medication=FALSE,
                                             consider.dosage.change=FALSE,
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
                                             remove.events.outside.followup.window=FALSE,
                                             parallel.backend="none", # make sure this runs sequentially!
                                             parallel.threads=1,
                                             suppress.warnings=suppress.warnings,
                                             suppress.special.argument.checks=suppress.special.argument.checks,
                                             return.data.table=TRUE);
    if( is.null(actual.obs.win) ) return (list("CMA"=NA, "event.info"=NULL));

    # Add the actual OW dates to event.info:
    actual.obs.win <- actual.obs.win[,c(ID.colname,".OBS.START.DATE",".OBS.END.DATE"),with=FALSE]
    setkeyv(actual.obs.win, ID.colname);
    event.info <- merge(event.info, actual.obs.win, all.x=TRUE); setnames(event.info, ncol(event.info)-c(1,0), c(".OBS.START.DATE.ACTUAL", ".OBS.END.DATE.ACTUAL"));

    CMA <- event.info[, .process.patient(.SD), by=ID.colname];

    # Make sure the information reflects the real observation window:
    event.info <- compute.event.int.gaps(data=as.data.frame(data),
                                         ID.colname=ID.colname,
                                         event.date.colname=event.date.colname,
                                         event.duration.colname=event.duration.colname,
                                         event.daily.dose.colname=event.daily.dose.colname,
                                         medication.class.colname=medication.class.colname,
                                         event.interval.colname=event.interval.colname,
                                         gap.days.colname=gap.days.colname,
                                         carryover.within.obs.window=TRUE,
                                         carryover.into.obs.window=TRUE, # if TRUE consider the carry-over from before the starting date of the observation window
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
                                         keep.event.interval.for.all.events=TRUE,
                                         parallel.backend="none", # make sure this runs sequentially!
                                         parallel.threads=1,
                                         suppress.warnings=suppress.warnings,
                                         suppress.special.argument.checks=suppress.special.argument.checks,
                                         return.data.table=TRUE);
    event.info$.EVENT.USED.IN.CMA <- (event.info$..ORIGINAL.ROW.ORDER.. %in% CMA$included.events.original.row.order); # store if the event was used to compute the CMA or not

    return (list("CMA"=unique(CMA[,1:2]), "event.info"=event.info));
  }

  ret.val <- .cma.skeleton(data=data,
                           ret.val=ret.val,
                           cma.class.name=c("CMA9","CMA1"),

                           ID.colname=ID.colname,
                           event.date.colname=event.date.colname,
                           event.duration.colname=event.duration.colname,
                           event.daily.dose.colname=event.daily.dose.colname,
                           medication.class.colname=medication.class.colname,
                           event.interval.colname=event.interval.colname,
                           gap.days.colname=gap.days.colname,
                           carryover.within.obs.window=TRUE,
                           carryover.into.obs.window=TRUE,
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
print.CMA9 <- function(...) print.CMA0(...)

#' @rdname plot.CMA1
#' @export
plot.CMA9 <- function(...) .plot.CMA1plus(...)