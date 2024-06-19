#' CMA_sliding_window constructor.
#'
#' Applies a given CMA to each sliding window and constructs a
#' CMA_sliding_window object.
#'
#' \code{CMA_sliding_window} first computes a set of fixed-size (possibly partly
#' overlapping) sliding windows,
#' each sliding to the right by a fixed timelag,
#' and then, for each of them, it computes the given "simple" CMA.
#' Thus, as opposed to the "simple" CMAs 1 to 9, it returns a set of CMAs, with
#' possibly more than one element.
#'
#' It is highly similar to \code{\link{CMA_per_episode}} which computes a CMA
#' for a set of treatment episodes.
#'
#' @param CMA.to.apply A \emph{string} giving the name of the CMA function (1 to
#' 9) that will be computed for each treatment episode.
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
#' @param sliding.window.start,sliding.window.start.unit,sliding.window.duration,sliding.window.duration.unit the definition of the first sliding window
#' (see the follow-up window parameters above for details).
#' @param sliding.window.step.duration,sliding.window.step.unit if not missing
#' (\code{NA}), these give the step (or "jump") to the right of the sliding
#' window in time units.
#' @param sliding.window.no.steps a \emph{integer} specifying the desired number
#' of sliding windows to cover the observation window (if possible); trumps
#' \code{sliding.window.step.duration} and \code{sliding.window.step.unit}.
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
#' @param return.inner.event.info \emph{Logical} specifying if the function
#' should also return the event.info for all the individual events in each
#' sliding window; by default it is \code{FALSE} as this information is useful
#' only in very specific cases (e.g., plotting the event intervals) and adds a
#' small but non-negligible computational overhead.
#' @param force.NA.CMA.for.failed.patients \emph{Logical} describing how the
#' patients for which the CMA estimation fails are treated: if \code{TRUE}
#' they are returned with an \code{NA} CMA estimate, while for
#' \code{FALSE} they are omitted.
#' @param return.mapping.events.sliding.window A \emph{Logical}, if \code{TRUE} then
#' the mapping between events and sliding windows is returned as the component
#' \code{mapping.windows.to.events}, which is a \code{data.table} giving, for
#' each sliding window, the events that belong to it (an event is given by its row
#' number in the \code{data}).
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
#' @param ... other possible parameters
#' @return An \code{S3} object of class \code{CMA_sliding_window} with the
#' following fields:
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
#'  reflect changes in dosage, as given by the \code{consider.dosage.change} parameter.
#'  \item \code{followup.window.start} the beginning of the follow-up window,
#'   as given by the \code{followup.window.start} parameter.
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
#'  \item \code{computed.CMA} the class name of the computed CMA.
#'  \item \code{CMA} the \code{data.frame} containing the actual \code{CMA}
#'  estimates for each participant (the \code{ID.colname} column) and sliding
#'  window, with columns:
#'    \itemize{
#'      \item \code{ID.colname} the patient ID as given by the \code{ID.colname}
#'      parameter.
#'      \item \code{window.ID} the unique window ID (within patients).
#'      \item \code{window.start} the window's start date (as a \code{Date}
#'      object).
#'      \item \code{window.end} the window's end date (as a \code{Date} object).
#'      \item \code{CMA} the window's estimated CMA.
#'    }
#' }
#' Please note that if \code{medication.groups} are defined, then the \code{CMA}
#' and \code{event.info} are named lists, each element containing the CMA and
#' event.info corresponding to a single medication group (the element's name).
#' @seealso \code{\link{CMA_per_episode}} is very similar, computing a "simple"
#' CMA for each of the treatment episodes.
#' The "simple" CMAs that can be computed comprise \code{\link{CMA1}},
#' \code{\link{CMA2}}, \code{\link{CMA3}}, \code{\link{CMA4}},
#' \code{\link{CMA5}}, \code{\link{CMA6}}, \code{\link{CMA7}},
#' \code{\link{CMA8}}, \code{\link{CMA9}}, as well as user-defined classes
#' derived from \code{\link{CMA0}} that have a \code{CMA} component giving the
#' estimated CMA per patient as a \code{data.frame}.
#' If \code{return.mapping.events.sliding.window} is \code{TRUE}, then this also has an
#' attribute \code{mapping.windows.to.events} that gives the mapping between
#' episodes and events as a \code{data.table} with the following columns:
#' \itemize{
#'  \item \code{patid} the patient ID.
#'  \item \code{window.ID} the episode unique ID (increasing sequentially).
#'  \item \code{event.index.in.data} the event given by its row number in the \code{data}.
#' }
#' @examples
#' \dontrun{
#' cmaW <- CMA_sliding_window(CMA="CMA1",
#'                            data=med.events,
#'                            ID.colname="PATIENT_ID",
#'                            event.date.colname="DATE",
#'                            event.duration.colname="DURATION",
#'                            event.daily.dose.colname="PERDAY",
#'                            medication.class.colname="CATEGORY",
#'                            carry.only.for.same.medication=FALSE,
#'                            consider.dosage.change=FALSE,
#'                            followup.window.start=0,
#'                            observation.window.start=0,
#'                            observation.window.duration=365,
#'                            sliding.window.start=0,
#'                            sliding.window.start.unit="days",
#'                            sliding.window.duration=90,
#'                            sliding.window.duration.unit="days",
#'                            sliding.window.step.duration=7,
#'                            sliding.window.step.unit="days",
#'                            sliding.window.no.steps=NA,
#'                            date.format="%m/%d/%Y"
#'                           );}
#' @export
CMA_sliding_window <- function( CMA.to.apply,  # the name of the CMA function (e.g., "CMA1") to be used
                                data, # the data used to compute the CMA on
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
                                carry.only.for.same.medication=NA, # if TRUE the carry-over applies only across medication of same type (NA = undefined)
                                consider.dosage.change=NA, # if TRUE carry-over is adjusted to reflect changes in dosage (NA = undefined)
                                # The follow-up window:
                                followup.window.start=0, # if a number is the earliest event per participant date + number of units, or a Date object, or a column name in data (NA = undefined)
                                followup.window.start.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!) (NA = undefined)
                                followup.window.start.per.medication.group=FALSE, # if there are medication groups and this is TRUE, then the first event is relative to each medication group separately, otherwise is relative to the patient
                                followup.window.duration=365*2, # the duration of the follow-up window in the time units given below (NA = undefined)
                                followup.window.duration.unit=c("days", "weeks", "months", "years")[1],# the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!)  (NA = undefined)
                                # The observation window (embedded in the follow-up window):
                                observation.window.start=0, # the number of time units relative to followup.window.start (NA = undefined)
                                observation.window.start.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!) (NA = undefined)
                                observation.window.duration=365*2, # the duration of the observation window in time units (NA = undefined)
                                observation.window.duration.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!) (NA = undefined)
                                # Sliding window:
                                sliding.window.start=0, # if a number is the earliest event per participant date + number of units, or a Date object, or a column name in data (NA = undefined)
                                sliding.window.start.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!) (NA = undefined)
                                sliding.window.duration=90,  # the duration of the sliding window in time units (NA = undefined)
                                sliding.window.duration.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!) (NA = undefined)
                                sliding.window.step.duration=30, # the step ("jump") of the sliding window in time units (NA = undefined)
                                sliding.window.step.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!) (NA = undefined)
                                sliding.window.no.steps=NA, # the number of steps to jump; has priority over setp specification
                                return.inner.event.info=FALSE, # should we return the event.info for all the individual events in each sliding window?
                                # Date format:
                                date.format="%m/%d/%Y", # the format of the dates used in this function (NA = undefined)
                                # Comments and metadata:
                                summary="CMA per sliding window",
                                # The description of the output (added) columns:
                                event.interval.colname="event.interval", # contains number of days between the start of current event and the start of the next
                                gap.days.colname="gap.days", # contains the number of days when medication was not available
                                # Dealing with failed estimates:
                                force.NA.CMA.for.failed.patients=TRUE, # force the failed patients to have NA CM estimate?
                                return.mapping.events.sliding.window=FALSE, # return the mapping of events to sliding windows?
                                # Parallel processing:
                                parallel.backend=c("none","multicore","snow","snow(SOCK)","snow(MPI)","snow(NWS)")[1], # parallel backend to use
                                parallel.threads="auto", # specification (or number) of parallel threads
                                # Misc:
                                suppress.warnings=FALSE,
                                suppress.special.argument.checks=FALSE, # used internally to suppress the check that we don't use special argument names
                                # extra parameters to be sent to the CMA function:
                                ...
)
{
  # Preconditions:
  if( is.numeric(sliding.window.start) && sliding.window.start < 0 )
  {
    if( !suppress.warnings ) .report.ewms("The sliding window must start a positive number of time units after the start of the observation window!\n", "error", "CMA_sliding_window", "AdhereR")
    return (NULL);
  }
  if( !inherits(sliding.window.start,"Date") && !is.numeric(sliding.window.start) && !(sliding.window.start %in% names(data)) )
  {
    if( !suppress.warnings ) .report.ewms("The sliding window start must be a valid column name!\n", "error", "CMA_sliding_window", "AdhereR")
    return (NULL);
  }
  if( !(sliding.window.start.unit %in% c("days", "weeks", "months", "years") ) )
  {
    if( !suppress.warnings ) .report.ewms("The sliding window start unit is not recognized!\n", "error", "CMA_sliding_window", "AdhereR")
    return (NULL);
  }
  if( !is.numeric(sliding.window.duration) || sliding.window.duration <= 0 )
  {
    if( !suppress.warnings ) .report.ewms("The sliding window duration must be greater than 0!\n", "error", "CMA_sliding_window", "AdhereR")
    return (NULL);
  }
  if( !(sliding.window.duration.unit %in% c("days", "weeks", "months", "years") ) )
  {
    if( !suppress.warnings ) .report.ewms("The sliding window duration unit is not recognized!\n", "error", "CMA_sliding_window", "AdhereR")
    return (NULL);
  }
  if( is.numeric(sliding.window.step.duration) && sliding.window.step.duration < 1 )
  {
    if( !suppress.warnings ) .report.ewms("The sliding window's step duration must be at least 1 unit!\n", "error", "CMA_sliding_window", "AdhereR")
    return (NULL);
  }
  if( !(sliding.window.step.unit %in% c("days", "weeks", "months", "years") ) )
  {
    if( !suppress.warnings ) .report.ewms("The sliding window's step duration unit is not recognized!\n", "error", "CMA_sliding_window", "AdhereR")
    return (NULL);
  }
  if( is.numeric(sliding.window.no.steps) && sliding.window.no.steps < 1 )
  {
    if( !suppress.warnings ) .report.ewms("The sliding window must move at least once!\n", "error", "CMA_sliding_window", "AdhereR")
    return (NULL);
  }

  # Get the CMA function corresponding to the name:
  if( !(is.character(CMA.to.apply) || is.factor(CMA.to.apply)) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("'CMA.to.apply' must be a string contining the name of the simple CMA to apply!\n)"), "error", "CMA_sliding_window", "AdhereR");
    return (NULL);
  }
  CMA.FNC <- switch(as.character(CMA.to.apply), # if factor, force it to string
                    "CMA1" = CMA1,
                    "CMA2" = CMA2,
                    "CMA3" = CMA3,
                    "CMA4" = CMA4,
                    "CMA5" = CMA5,
                    "CMA6" = CMA6,
                    "CMA7" = CMA7,
                    "CMA8" = CMA8,
                    "CMA9" = CMA9,
                    {if( !suppress.warnings ) .report.ewms(paste0("Unknown 'CMA.to.apply' '",CMA.to.apply,"': defaulting to CMA0!\n)"), "warning", "CMA_sliding_window", "AdhereR"); CMA0;}); # by default, fall back to CMA0

  # Default argument values and overrides:
  def.vals <- formals(CMA.FNC);
  if( CMA.to.apply %in% c("CMA1", "CMA2", "CMA3", "CMA4") )
  {
    carryover.into.obs.window <- carryover.within.obs.window <- FALSE;
    if( !is.na(carry.only.for.same.medication) && carry.only.for.same.medication && !suppress.warnings ) .report.ewms("'carry.only.for.same.medication' cannot be defined for CMAs 1-4!\n", "warning", "CMA_sliding_window", "AdhereR");
    carry.only.for.same.medication <- FALSE;
    if( !is.na(consider.dosage.change) && consider.dosage.change && !suppress.warnings ) .report.ewms("'consider.dosage.change' cannot be defined for CMAs 1-4!\n", "warning", "CMA_sliding_window", "AdhereR");
    consider.dosage.change <- FALSE;
  } else if( CMA.to.apply %in% c("CMA5", "CMA6") )
  {
    carryover.into.obs.window <- FALSE;
    carryover.within.obs.window <- TRUE;
    if( is.na(carry.only.for.same.medication) ) carry.only.for.same.medication <- def.vals[["carry.only.for.same.medication"]]; # use the default value from CMA
    if( is.na(consider.dosage.change) ) consider.dosage.change <- def.vals[["consider.dosage.change"]]; # use the default value from CMA
  } else if( CMA.to.apply %in% c("CMA7", "CMA8", "CMA9") )
  {
    carryover.into.obs.window <- carryover.within.obs.window <- TRUE;
    if( is.na(carry.only.for.same.medication) ) carry.only.for.same.medication <- def.vals[["carry.only.for.same.medication"]]; # use the default value from CMA
    if( is.na(consider.dosage.change) ) consider.dosage.change <- def.vals[["consider.dosage.change"]]; # use the default value from CMA
  } else
  {
    if( !suppress.warnings ) .report.ewms("I know how to do CMA sliding windows only for CMAs 1 to 9!\n", "error", "CMA_sliding_window", "AdhereR");
    return (NULL);
  }

  # Create the return value skeleton and check consistency:
  ret.val <- CMA0(data,
                  ID.colname=ID.colname,
                  event.date.colname=event.date.colname,
                  event.duration.colname=event.duration.colname,
                  event.daily.dose.colname=event.daily.dose.colname,
                  medication.class.colname=medication.class.colname,
                  medication.groups=medication.groups,
                  flatten.medication.groups=flatten.medication.groups,
                  medication.groups.colname=medication.groups.colname,
                  carryover.within.obs.window=carryover.within.obs.window,
                  carryover.into.obs.window=carryover.into.obs.window,
                  carry.only.for.same.medication=carry.only.for.same.medication,
                  consider.dosage.change=consider.dosage.change,
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
                  suppress.warnings=suppress.warnings,
                  suppress.special.argument.checks=suppress.special.argument.checks,
                  summary=NA);
  if( is.null(ret.val) ) return (NULL);
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
      if( n.events < 1 ) return (NULL);

      # Compute the sliding windows for this patient:
      start.date <- .add.time.interval.to.date(data4ID$.OBS.START.DATE[1], sliding.window.start, sliding.window.start.unit, suppress.warnings); # when do the windows start?
      sliding.duration <- as.numeric(data4ID$.OBS.END.DATE[1] - start.date) - sliding.window.duration.in.days; # the effective duration to be covered with sliding windows
      if( sliding.duration < 0)  return (NULL); # the sliding window is longer than the available time in the observation window
      if( is.na(sliding.window.no.steps) )
      {
        # Compute the number of steps required from the step size:
        sliding.window.no.steps <- (sliding.duration / sliding.window.step.duration.in.days) + 1;
      } else if( sliding.window.no.steps > 1 )
      {
        # Compute the step size to optimally cover the duration (possibly adjust the number of steps, too):
        sliding.window.step.duration.in.days <- (sliding.duration / (sliding.window.no.steps - 1));
        sliding.window.step.duration.in.days <- max(1, min(sliding.window.duration.in.days, sliding.window.step.duration.in.days)); # make sure we don't overdue it
        sliding.window.no.steps <- min(((sliding.duration / sliding.window.step.duration.in.days) + 1), sliding.duration); # adjust the number of steps just in case
      } else
      {
        # Only one sliding window:
        sliding.window.step.duration.in.days <- 0;
      }
      if( sliding.window.no.steps < 1 || sliding.window.step.duration.in.days < 0 ) return (NULL);

      # Expand the participant data by replicating it for each consecutive sliding window:
      data4ID.wnds <- cbind(".WND.ID"        =rep(1:sliding.window.no.steps, each=n.events),
                            data4ID,
                            ".WND.START.DATE"=rep(as.Date(vapply(1:sliding.window.no.steps,
                                                                 function(i) .add.time.interval.to.date(start.date, (i-1)*sliding.window.step.duration.in.days, "days", suppress.warnings),
                                                                 numeric(1)),
                                                          origin=lubridate::origin),
                                                      each=n.events),
                            ".WND.DURATION"  =sliding.window.duration.in.days);
      setkeyv(data4ID.wnds, ".WND.ID");

      # Apply the desired CMA to all the windows:
      if(length(dot.args <- list(...)) > 0 && "arguments.that.should.not.be.defined" %in% names(dot.args)) # check if arguments.that.should.not.be.defined is passed in the ... argument
      {
        # Avoid passing again arguments.that.should.not.be.defined to the CMA function:
        cma <- CMA.FNC(data=as.data.frame(data4ID.wnds),
                       ID.colname=".WND.ID",
                       event.date.colname=event.date.colname,
                       event.duration.colname=event.duration.colname,
                       event.daily.dose.colname=event.daily.dose.colname,
                       medication.class.colname=medication.class.colname,
                       carry.only.for.same.medication=carry.only.for.same.medication,
                       consider.dosage.change=consider.dosage.change,
                       followup.window.start=followup.window.start,
                       followup.window.start.unit=followup.window.start.unit,
                       followup.window.duration=followup.window.duration,
                       followup.window.duration.unit=followup.window.duration.unit,
                       observation.window.start=".WND.START.DATE",
                       observation.window.duration=".WND.DURATION",
                       observation.window.duration.unit="days",
                       date.format=date.format,
                       parallel.backend="none", # make sure this runs sequentially!
                       parallel.threads=1,
                       suppress.warnings=suppress.warnings,
                       suppress.special.argument.checks=TRUE, # we know we're using special columns, so supress this check
                       ...);
      } else
      {
        # Temporarily avoid warnings linked to rewriting some CMA arguments by passing it arguments.that.should.not.be.defined=NULL:
        cma <- CMA.FNC(data=as.data.frame(data4ID.wnds),
                       ID.colname=".WND.ID",
                       event.date.colname=event.date.colname,
                       event.duration.colname=event.duration.colname,
                       event.daily.dose.colname=event.daily.dose.colname,
                       medication.class.colname=medication.class.colname,
                       carry.only.for.same.medication=carry.only.for.same.medication,
                       consider.dosage.change=consider.dosage.change,
                       followup.window.start=followup.window.start,
                       followup.window.start.unit=followup.window.start.unit,
                       followup.window.duration=followup.window.duration,
                       followup.window.duration.unit=followup.window.duration.unit,
                       observation.window.start=".WND.START.DATE",
                       observation.window.duration=".WND.DURATION",
                       observation.window.duration.unit="days",
                       date.format=date.format,
                       parallel.backend="none", # make sure this runs sequentially!
                       parallel.threads=1,
                       suppress.warnings=suppress.warnings,
                       suppress.special.argument.checks=TRUE, # we know we're using special columns, so supress this check
                       arguments.that.should.not.be.defined=NULL,
                       ...);
      }
      if( is.null(cma) ) return (NULL);

      # Unpack the data fo returning:
      wnd.info <- data4ID.wnds[,
                               list(".WND.START"=.WND.START.DATE[1], # the start
                                    ".WND.END"=.add.time.interval.to.date(.WND.START.DATE[1], sliding.window.duration.in.days, "days", suppress.warnings)), # and end
                               by=.WND.ID]; # for each window
      wnd.info <- cbind(merge(wnd.info, cma$CMA, by=".WND.ID", all=TRUE),
                        "CMA.to.apply"=class(cma)[1]);
      setnames(wnd.info, c("window.ID", "window.start", "window.end", "CMA", "CMA.to.apply"));

      cma$event.info$..ORIGINAL.ROW.ORDER.IN.DATA.. <- data4ID.wnds$..ORIGINAL.ROW.ORDER..[ cma$event.info$..ORIGINAL.ROW.ORDER.. ]; # restore the row numbers in the original data

      if( return.inner.event.info || return.mapping.events.sliding.window )
      {
        # Add the event.info to also return it:
        wnd.info <- merge(wnd.info, cma$event.info, by.x="window.ID", by.y=".WND.ID", all=TRUE);
      }

      return (as.data.frame(wnd.info));
    }

    # Compute the real observation windows (might differ per patient) only once per patient (speed things up & the observation window is the same for all events within a patient):
    #patinfo.cols <- which(names(data) %in% c(ID.colname, event.date.colname, event.duration.colname, event.daily.dose.colname, medication.class.colname));
    #tmp <- as.data.frame(data); tmp <- tmp[!duplicated(tmp[,ID.colname]),patinfo.cols]; # the reduced dataset for computing the actual OW:
    tmp <- as.data.frame(data); tmp <- tmp[!duplicated(tmp[,ID.colname]),]; # the reduced dataset for computing the actual OW:
    event.info2 <- compute.event.int.gaps(data=tmp,
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
    if( is.null(event.info2) ) return (NULL);
    # Merge the observation window start and end dates back into the data:
    data <- merge(data, event.info2[,c(ID.colname, ".OBS.START.DATE", ".OBS.END.DATE"),with=FALSE], all.x=TRUE);
    setnames(data, ncol(data)-c(1,0), c(".OBS.START.DATE.PRECOMPUTED", ".OBS.END.DATE.PRECOMPUTED"));

    if( return.inner.event.info || return.mapping.events.sliding.window )
    {
      CMA.and.event.info <- data[, .process.patient(.SD), by=ID.colname ]; # contains both the CMAs per sliding window and the event.info's for each event
      CMA <- unique(CMA.and.event.info[,1:6,with=FALSE]); # the per-window CMAs...
      inner.event.info <- CMA.and.event.info[,c(1,2,7:ncol(CMA.and.event.info)),with=FALSE]; # ...and the event info for each event within the sliding windows
      return (list("CMA"=CMA, "event.info"=event.info2[,c(ID.colname, ".FU.START.DATE", ".FU.END.DATE", ".OBS.START.DATE", ".OBS.END.DATE"), with=FALSE], "inner.event.info"=inner.event.info));
    } else
    {
      CMA <- data[, .process.patient(.SD), by=ID.colname ]; # contains just the CMAs per sliding window
      return (list("CMA"=CMA, "event.info"=event.info2[,c(ID.colname, ".FU.START.DATE", ".FU.END.DATE", ".OBS.START.DATE", ".OBS.END.DATE"), with=FALSE]));
    }
  }

  # Convert the sliding window duration and step in days:
  sliding.window.duration.in.days <- switch(sliding.window.duration.unit,
                                            "days"=sliding.window.duration,
                                            "weeks"=sliding.window.duration * 7,
                                            "months"=sliding.window.duration * 30,
                                            "years"=sliding.window.duration * 365,
                                            sliding.window.duration);
  sliding.window.step.duration.in.days <- switch(sliding.window.step.unit,
                                                 "days"=sliding.window.step.duration,
                                                 "weeks"=sliding.window.step.duration * 7,
                                                 "months"=sliding.window.step.duration * 30,
                                                 "years"=sliding.window.step.duration * 365,
                                                 sliding.window.step.duration);

  # Convert to data.table, cache event date as Date objects, and key by patient ID and event date
  data.copy <- data.table(data);
  data.copy[, .DATE.as.Date := as.Date(get(event.date.colname),format=date.format)]; # .DATE.as.Date: convert event.date.colname from formatted string to Date
  data.copy$..ORIGINAL.ROW.ORDER.. <- 1:nrow(data.copy); # preserve the original order of the rows (needed for medication groups)
  setkeyv(data.copy, c(ID.colname, ".DATE.as.Date")); # key (and sorting) by patient ID and event date

  ## Computing the inner.event.info makes sense only for non-overlapping sliding windows:
  #if( return.inner.event.info && (sliding.window.duration.in.days > sliding.window.step.duration.in.days) )
  #{
  #  if( !suppress.warnings ) .report.ewms("Computing the inner.event.info makes sense only for non-overlapping sliding windows.\n", "warning", "CMA_sliding_windows", "AdhereR");
  #  return.inner.event.info <- FALSE;
  #}


  # Are there medication groups?
  if( is.null(mg <- getMGs(ret.val)) )
  {
    # Nope: do a single estimation on the whole dataset:

    # Compute the workhorse function:
    tmp <- .compute.function(.workhorse.function, fnc.ret.vals=2,
                             parallel.backend=parallel.backend,
                             parallel.threads=parallel.threads,
                             data=data.copy,
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
                             suppress.warnings=suppress.warnings,
                             suppress.special.argument.checks=suppress.special.argument.checks);
    if( is.null(tmp) || is.null(tmp$CMA) || !inherits(tmp$CMA,"data.frame") || is.null(tmp$event.info) ) return (NULL);

    # Construct the return object:
    class(ret.val) <- "CMA_sliding_window";
    ret.val$event.info <- as.data.frame(tmp$event.info);
    ret.val$computed.CMA <- as.character(tmp$CMA$CMA.to.apply[1]);
    ret.val$sliding.window.start <- sliding.window.start;
    ret.val$sliding.window.start.unit <- sliding.window.start.unit;
    ret.val$sliding.window.duration <- sliding.window.duration;
    ret.val$sliding.window.duration.unit <- sliding.window.duration.unit;
    ret.val$sliding.window.step.duration <- sliding.window.step.duration;
    ret.val$sliding.window.step.unit <- sliding.window.step.unit;
    ret.val$sliding.window.no.steps <- sliding.window.no.steps;
    ret.val$summary <- summary;
    ret.val$CMA <- as.data.frame(tmp$CMA); setnames(ret.val$CMA, 1, ID.colname); ret.val$CMA <- ret.val$CMA[,-ncol(ret.val$CMA)];
    if( return.inner.event.info && !is.null(tmp$inner.event.info) ) ret.val$inner.event.info <- as.data.frame(tmp$inner.event.info);
    if( return.mapping.events.sliding.window && !is.null(tmp$inner.event.info) && ".EVENT.USED.IN.CMA" %in% names(tmp$inner.event.info) )
    {
      # Unpack the info in the needed format:
      ret.val$mapping.windows.to.events <- as.data.frame(tmp$inner.event.info)[ tmp$inner.event.info$.EVENT.USED.IN.CMA, c(ID.colname, "window.ID", "..ORIGINAL.ROW.ORDER.IN.DATA..") ];
      names(ret.val$mapping.windows.to.events)[3] <- "event.index.in.data";
    }
    return (ret.val);

  } else
  {
    # Yes

    # Make sure the group's observations reflect the potentially new order of the observations in the data:
    mb.obs <- mg$obs[data.copy$..ORIGINAL.ROW.ORDER.., ];

    # Focus only on the non-trivial ones:
    mg.to.eval <- (colSums(!is.na(mb.obs) & mb.obs) > 0);
    if( sum(mg.to.eval) == 0 )
    {
      # None selects not even one observation!
      .report.ewms(paste0("None of the medication classes (included __ALL_OTHERS__) selects any observation!\n"), "warning", "CMA1", "AdhereR");
      return (NULL);
    }
    mb.obs <- mb.obs[,mg.to.eval]; # keep only the non-trivial ones

    # Check if there are medication classes that refer to the same observations (they would result in the same estimates):
    mb.obs.dupl <- duplicated(mb.obs, MARGIN=2);

    # Estimate each separately:
    tmp <- lapply(1:nrow(mg$defs), function(i)
    {
      # Check if these are to be evaluated:
      if( !mg.to.eval[i] )
      {
        return (list("CMA"=NULL, "event.info"=NULL, "CMA.to.apply"=NA));
      }

      # Translate into the index of the classes to be evaluated:
      ii <- sum(mg.to.eval[1:i]);

      # Cache the selected observations:
      mg.sel.obs <- mb.obs[,ii];

      # Check if this is a duplicated medication class:
      if( mb.obs.dupl[ii] )
      {
        # Find which one is the original:
        for( j in 1:(ii-1) ) # ii=1 never should be TRUE
        {
          if( identical(mb.obs[,j], mg.sel.obs) )
          {
            # This is the original: return it and stop
            return (c("identical.to"=j));
          }
        }
      }

      # Compute the workhorse function:
      tmp <- .compute.function(.workhorse.function, fnc.ret.vals=2,
                               parallel.backend=parallel.backend,
                               parallel.threads=parallel.threads,
                               data=data.copy[mg.sel.obs,], # apply it on the subset of observations covered by this medication class
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
                               suppress.warnings=suppress.warnings,
                               suppress.special.argument.checks=suppress.special.argument.checks);
      if( is.null(tmp) || is.null(tmp$CMA) || !inherits(tmp$CMA,"data.frame") || is.null(tmp$event.info) ) return (NULL);

      # Convert to data.frame and return:
      tmp$CMA.to.apply <- tmp$CMA$CMA.to.apply[1];
      tmp$CMA <- as.data.frame(tmp$CMA); setnames(tmp$CMA, 1, ID.colname); tmp$CMA <- tmp$CMA[,-ncol(tmp$CMA)];
      tmp$event.info <- as.data.frame(tmp$event.info);
      if( return.inner.event.info && !is.null(tmp$inner.event.info) ) tmp$inner.event.info <- as.data.frame(tmp$inner.event.info);
      return (tmp);

    });

    # Set the names:
    names(tmp) <- mg$defs$name;

    # Solve the duplicates:
    for( i in seq_along(tmp) )
    {
      if( is.numeric(tmp[[i]]) && length(tmp[[i]]) == 1 && names(tmp[[i]]) == "identical.to" ) tmp[[i]] <- tmp[[ tmp[[i]] ]];
    }

    # Rearrange these and return:
    ret.val[["CMA"]]        <- lapply(tmp, function(x) x$CMA);
    ret.val[["event.info"]] <- lapply(tmp, function(x) x$event.info);
    if( return.inner.event.info && !is.null(tmp$inner.event.info) ) ret.val[["inner.event.info"]] <- lapply(tmp, function(x) x$inner.event.info);
    ret.val$computed.CMA <- unique(vapply(tmp, function(x) if(is.null(x) || is.na(x$CMA.to.apply)){return (NA_character_)}else{return(x$CMA.to.apply)}, character(1))); ret.val$computed.CMA <- ret.val$computed.CMA[ !is.na(ret.val$computed.CMA) ];
    if( flatten.medication.groups && !is.na(medication.groups.colname) )
    {
      # Flatten the CMA:
      tmp <- do.call(rbind, ret.val[["CMA"]]);
      if( is.null(tmp) || nrow(tmp) == 0 )
      {
        ret.val[["CMA"]] <- NULL;
      } else
      {
        tmp <- cbind(tmp, unlist(lapply(1:length(ret.val[["CMA"]]), function(i) if(!is.null(ret.val[["CMA"]][[i]])){rep(names(ret.val[["CMA"]])[i], nrow(ret.val[["CMA"]][[i]]))}else{NULL})));
        names(tmp)[ncol(tmp)] <- medication.groups.colname; rownames(tmp) <- NULL;
        ret.val[["CMA"]] <- tmp;
      }

      # ... and the event.info:
      tmp <- do.call(rbind, ret.val[["event.info"]]);
      if( is.null(tmp) || nrow(tmp) == 0 )
      {
        ret.val[["event.info"]] <- NULL;
      } else
      {
        tmp <- cbind(tmp, unlist(lapply(1:length(ret.val[["event.info"]]), function(i) if(!is.null(ret.val[["event.info"]][[i]])){rep(names(ret.val[["event.info"]])[i], nrow(ret.val[["event.info"]][[i]]))}else{NULL})));
        names(tmp)[ncol(tmp)] <- medication.groups.colname; rownames(tmp) <- NULL;
        ret.val[["event.info"]] <- tmp;
      }

      # ... and the inner.event.info:
      if( return.inner.event.info && !is.null(ret.val[["inner.event.info"]]) )
      {
        tmp <- do.call(rbind, ret.val[["inner.event.info"]]);
        if( is.null(tmp) || nrow(tmp) == 0 )
        {
          ret.val[["inner.event.info"]] <- NULL;
        } else
        {
          tmp <- cbind(tmp, unlist(lapply(1:length(ret.val[["inner.event.info"]]), function(i) if(!is.null(ret.val[["inner.event.info"]][[i]])){rep(names(ret.val[["inner.event.info"]])[i], nrow(ret.val[["inner.event.info"]][[i]]))}else{NULL})));
          names(tmp)[ncol(tmp)] <- medication.groups.colname; rownames(tmp) <- NULL;
          ret.val[["inner.event.info"]] <- tmp;
        }
      }
    }
    class(ret.val) <- "CMA_sliding_window";
    ret.val$sliding.window.start.unit <- sliding.window.start.unit;
    ret.val$sliding.window.duration <- sliding.window.duration;
    ret.val$sliding.window.duration.unit <- sliding.window.duration.unit;
    ret.val$sliding.window.step.duration <- sliding.window.step.duration;
    ret.val$sliding.window.step.unit <- sliding.window.step.unit;
    ret.val$sliding.window.no.steps <- sliding.window.no.steps;
    ret.val$summary <- summary;
    return (ret.val);

  }
}

#' @export
getMGs.CMA_sliding_window <- function(x)
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA_sliding_window") || is.null(cma$medication.groups) ) return (NULL);
  return (cma$medication.groups);
}

#' getEventsToSlidingWindowsMapping
#'
#' This function returns the event-to-sliding-window mapping, if this information exists.
#'
#' There are cases where it is interesting to know which events belong to which
#' sliding windows and which events have been actually used in computing the simple CMA
#' for each sliding window
#' This information can be returned by \code{CMA_sliding_window()} if the
#' parameter \code{return.mapping.events.episodes} is set to \code{TRUE}.
#'
#' @param x is an \code{CMA_sliding_window object}.
#' @return The mapping between events and episodes, if it exists as the
#' \code{mapping.windows.to.events} component of the \code{CMA_sliding_window object}
#' object, or \code{NULL} otherwise.
#' @export
getEventsToSlidingWindowsMapping <- function(x)
{
  if( is.null(x) )
  {
    return (NULL);
  } else if( inherits(x, "CMA_sliding_window") && ("mapping.windows.to.events" %in% names(x)) && !is.null(x$mapping.windows.to.events) && inherits(x$mapping.windows.to.events, "data.frame") )
  {
    return (x$mapping.windows.to.events);
  } else
  {
    return (NULL);
  }
}

#' @export
getCMA.CMA_sliding_window <- function(x, flatten.medication.groups=FALSE, medication.groups.colname=".MED_GROUP_ID")
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA_sliding_window") || !("CMA" %in% names(cma)) || is.null(cma$CMA) ) return (NULL);
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

#' @export
getEventInfo.CMA_sliding_window <- function(x, flatten.medication.groups=FALSE, medication.groups.colname=".MED_GROUP_ID")
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA_sliding_window") || !("event.info" %in% names(cma)) || is.null(cma$event.info) ) return (NULL);
  if( inherits(cma$event.info, "data.frame") || !flatten.medication.groups )
  {
    return (cma$event.info);
  } else
  {
    # Flatten the medication groups into a single data.frame:
    ret.val <- do.call(rbind, cma$event.info);
    if( is.null(ret.val) || nrow(ret.val) == 0 ) return (NULL);
    ret.val <- cbind(ret.val, unlist(lapply(1:length(cma$event.info), function(i) if(!is.null(cma$event.info[[i]])){rep(names(cma$event.info)[i], nrow(cma$event.info[[i]]))}else{NULL})));
    names(ret.val)[ncol(ret.val)] <- medication.groups.colname; rownames(ret.val) <- NULL;
    return (ret.val);
  }
}

#' @export
getInnerEventInfo.CMA_sliding_window <- function(x, flatten.medication.groups=FALSE, medication.groups.colname=".MED_GROUP_ID")
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA_sliding_window") || !("event.info" %in% names(cma)) || is.null(cma$inner.event.info) ) return (NULL);
  if( inherits(cma$inner.event.info, "data.frame") || !flatten.medication.groups )
  {
    return (cma$inner.event.info);
  } else
  {
    # Flatten the medication groups into a single data.frame:
    ret.val <- do.call(rbind, cma$inner.event.info);
    if( is.null(ret.val) || nrow(ret.val) == 0 ) return (NULL);
    ret.val <- cbind(ret.val, unlist(lapply(1:length(cma$inner.event.info), function(i) if(!is.null(cma$inner.event.info[[i]])){rep(names(cma$inner.event.info)[i], nrow(cma$inner.event.info[[i]]))}else{NULL})));
    names(ret.val)[ncol(ret.val)] <- medication.groups.colname; rownames(ret.val) <- NULL;
    return (ret.val);
  }
}

#' @export
subsetCMA.CMA_sliding_window <- function(cma, patients, suppress.warnings=FALSE)
{
  if( inherits(patients, "factor") ) patients <- as.character(patients);
  patients.to.keep <- intersect(patients, unique(cma$data[,cma$ID.colname]));
  if( length(patients.to.keep) == 0 )
  {
    if( !suppress.warnings ) .report.ewms("No patients to subset on!\n", "error", "subsetCMA.CMA_sliding_window", "AdhereR");
    return (NULL);
  }
  if( length(patients.to.keep) < length(patients) && !suppress.warnings ) .report.ewms("Some patients in the subsetting set are not in the CMA itsefl and are ignored!\n", "warning", "subsetCMA.CMA_sliding_window", "AdhereR");

  ret.val <- cma;
  ret.val$data <- ret.val$data[ ret.val$data[,ret.val$ID.colname] %in% patients.to.keep, ];
  if( !is.null(ret.val$event.info) )
  {
    if( inherits(ret.val$event.info, "data.frame") )
    {
      ret.val$event.info <- ret.val$event.info[ ret.val$event.info[,ret.val$ID.colname] %in% patients.to.keep, ]; if( nrow(ret.val$event.info) == 0 ) ret.val$event.info <- NULL;
    } else if( is.list(ret.val$event.info) && length(ret.val$event.info) > 0 )
    {
      ret.val$event.info <- lapply(ret.val$event.info, function(x){tmp <- x[ x[,ret.val$ID.colname] %in% patients.to.keep, ]; if(!is.null(tmp) && nrow(tmp) > 0){tmp}else{NULL}});
    }
  }
  if( !is.null(ret.val$inner.event.info) )
  {
    if( inherits(ret.val$inner.event.info, "data.frame") )
    {
      ret.val$inner.event.info <- ret.val$inner.event.info[ ret.val$inner.event.info[,ret.val$ID.colname] %in% patients.to.keep, ]; if( nrow(ret.val$inner.event.info) == 0 ) ret.val$inner.event.info <- NULL;
    } else if( is.list(ret.val$inner.event.info) && length(ret.val$inner.event.info) > 0 )
    {
      ret.val$inner.event.info <- lapply(ret.val$inner.event.info, function(x){tmp <- x[ x[,ret.val$ID.colname] %in% patients.to.keep, ]; if(!is.null(tmp) && nrow(tmp) > 0){tmp}else{NULL}});
    }
  }
  if( ("CMA" %in% names(ret.val)) && !is.null(ret.val$CMA) )
  {
    if( inherits(ret.val$CMA, "data.frame") )
    {
      ret.val$CMA <- ret.val$CMA[ ret.val$CMA[,ret.val$ID.colname] %in% patients.to.keep, ];
    } else if( is.list(ret.val$CMA) && length(ret.val$CMA) > 0 )
    {
      ret.val$CMA <- lapply(ret.val$CMA, function(x){tmp <- x[ x[,ret.val$ID.colname] %in% patients.to.keep, ]; if(!is.null(tmp) && nrow(tmp) > 0){tmp}else{NULL}});
    }
  }
  return (ret.val);
}


#' @rdname print.CMA0
#' @export
print.CMA_sliding_window <- function(...) print.CMA_per_episode(...)

#' @rdname plot.CMA_per_episode
#' @export
plot.CMA_sliding_window <- plot.CMA_per_episode;