#' CMA_per_episode constructor.
#'
#' Applies a given CMA to each treatment episode and constructs a
#' CMA_per_episode object.
#'
#' \code{CMA_per_episode} first identifies the treatment episodes for the whole
#' follo-up window (using the \code{\link{compute.treatment.episodes}} function),
#' and then computes the given "simple" CMA for each treatment episode that
#' intersects with the observation window. NB: the CMA is computed for the
#' period of the episode that is part of the observations window; thus, if an
#' episode starts earlier or ends later than the observation window, CMA will
#' be computed for a section of that episode.
#' Thus, as opposed to the "simple" CMAs 1 to 9, it returns a set of CMAs, with
#' possibly more than one element.
#'
#' It is highly similar to \code{\link{CMA_sliding_window}} which computes a CMA
#' for a set of sliding windows.
#'
#' @param CMA.to.apply A \emph{string} giving the name of the CMA function (1 to
#' 9) that will be computed for each treatment episode.
#' @param data A \emph{\code{data.frame}} containing the events (prescribing or
#' dispensing) used to compute the CMA. Must contain, at a minimum, the patient
#' unique ID, the event date and duration, and might also contain the daily
#' dosage and medication type (the actual column names are defined in the
#' following four parameters).
#' @param treat.epi A \emph{\code{data.frame}} containing the treatment episodes.
#' Must contain the patient ID (as given in \code{ID.colname}), the episode unique ID
#' (increasing sequentially, \code{episode.ID}), the episode start date
#' (\code{episode.start}), the episode duration in days (\code{episode.duration}),
#' and the episode end date (\code{episode.end}).
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
#' carry-over applies only across medication of the same type; valid only for
#' CMAs 5 to 9, in which case it is coupled (i.e., the same value is used for
#' computing the treatment episodes and the CMA on each treatment episode).
#' @param consider.dosage.change \emph{Logical}, if \code{TRUE}, the carry-over
#' is adjusted to also reflect changes in dosage; valid only for CMAs 5 to 9, in
#' which case it is coupled (i.e., the same value is used for computing the
#' treatment episodes and the CMA on each treatment episode).
#' @param medication.change.means.new.treatment.episode \emph{Logical}, should a
#' change in medication automatically start a new treatment episode?
#' @param dosage.change.means.new.treatment.episode \emph{Logical}, should a
#' change in dosage automatically start a new treatment episode?
#' @param maximum.permissible.gap The \emph{number} of units given by
#' \code{maximum.permissible.gap.unit} representing the maximum duration of
#' permissible gaps between treatment episodes (can also be a percent, see
#' \code{maximum.permissible.gap.unit} for details).
#' @param maximum.permissible.gap.unit can be either \emph{"days"},
#' \emph{"weeks"}, \emph{"months"}, \emph{"years"} or \emph{"percent"}, and
#' represents the time units that \code{maximum.permissible.gap} refers to;
#' if \emph{percent}, then  \code{maximum.permissible.gap} is interpreted as a
#' percent (can be greater than 100\%) of the duration of the current
#' prescription.
#' @param maximum.permissible.gap.append.to.episode a \emph{logical} value
#' specifying of the \code{maximum.permissible.gap} should be append at the
#' end of an episode with a gap larger than the \code{maximum.permissible.gap};
#' \code{FALSE} (the default) mean no addition, while \code{TRUE} mean that the
#' full \code{maximum.permissible.gap} is added.
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
#' @param return.inner.event.info \emph{Logical} specifying if the function
#' should also return the event.info for all the individual events in each
#' sliding window; by default it is \code{FALSE} as this information is useful
#' only in very specific cases (e.g., plotting the event intervals) and adds a
#' small but non-negligible computational overhead.
#' @param force.NA.CMA.for.failed.patients \emph{Logical} describing how the
#' patients for which the CMA estimation fails are treated: if \code{TRUE}
#' they are returned with an \code{NA} CMA estimate, while for
#' \code{FALSE} they are omitted.
#' @param return.mapping.events.episodes A \emph{Logical}, if \code{TRUE} then
#' the mapping between events and episodes is returned as an extra component
#' \code{mapping.episodes.to.events}, which is a \code{data.table} giving, for
#' each episode, the events that belong to it (an event is given by its row
#' number in the \code{data}). Please note that the episodes returned are
#' specific to the particular simple CMA used, and should preferentially used
#' over those returned by \code{compute.treatment.episodes()}. This component
#' can also be accessed using the \code{getEventsToEpisodesMapping()} function.
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
#' @return An \code{S3} object of class \code{CMA_per_episode} with the
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
#'  reflect changes in dosage, as given by the \code{consider.dosage.change}
#'  parameter.
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
#'  \item \code{summary} the metadata, as given by the \code{summary} parameter.
#'  \item \code{event.info} the \code{data.frame} containing the event info
#'  (irrelevant for most users; see \code{\link{compute.event.int.gaps}} for
#'  details).
#'  \item \code{computed.CMA} the class name of the computed CMA.
#'  \item \code{CMA} the \code{data.frame} containing the actual \code{CMA}
#'  estimates for each participant (the \code{ID.colname} column) and treatment
#'  episode, with columns:
#'    \itemize{
#'      \item \code{ID.colname} the patient ID as given by the \code{ID.colname}
#'      parameter.
#'      \item \code{episode.ID} the unique treatment episode ID (within
#'      patients).
#'      \item \code{episode.start} the treatment episode's start date (as a
#'      \code{Date} object).
#'      \item \code{end.episode.gap.days} the corresponding gap days of the last
#'      event in this episode.
#'      \item \code{episode.duration} the treatment episode's duration in days.
#'      \item \code{episode.end} the treatment episode's end date (as a
#'      \code{Date} object).
#'      \item \code{CMA} the treatment episode's estimated CMA.
#'    }
#' }
#' Please note that if \code{medication.groups} are defined, then the \code{CMA}
#' and \code{event.info} are named lists, each element containing the CMA and
#' event.info corresponding to a single medication group (the element's name).
#' @seealso \code{\link{CMA_sliding_window}} is very similar, computing a
#' "simple" CMA for each of a set of same-size sliding windows.
#' The "simple" CMAs that can be computed comprise \code{\link{CMA1}},
#' \code{\link{CMA2}}, \code{\link{CMA3}}, \code{\link{CMA4}},
#' \code{\link{CMA5}}, \code{\link{CMA6}}, \code{\link{CMA7}},
#' \code{\link{CMA8}}, \code{\link{CMA9}}, as well as user-defined classes
#' derived from \code{\link{CMA0}} that have a \code{CMA} component giving the
#' estimated CMA per patient as a \code{data.frame}.
#' If \code{return.mapping.events.episodes} is \code{TRUE}, then this also has a
#' component \code{mapping.episodes.to.events} that gives the mapping between
#' episodes and events as a \code{data.table} with the following columns:
#' \itemize{
#'  \item \code{patid} the patient ID.
#'  \item \code{episode.ID} the episode unique ID (increasing sequentially).
#'  \item \code{event.index.in.data} the event given by its row number in the \code{data}.
#' }
#' @examples
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
#'                        );}
#' @export
CMA_per_episode <- function( CMA.to.apply,  # the name of the CMA function (e.g., "CMA1") to be used
                             data, # the data used to compute the CMA on
                             treat.epi=NULL, # the treatment episodes, if available
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
                             carry.only.for.same.medication=NA, # if TRUE the carry-over applies only across medication of same type (NA = use the CMA's values)
                             consider.dosage.change=NA, # if TRUE carry-over is adjusted to reflect changes in dosage (NA = use the CMA's values)
                             # Treatment episodes:
                             medication.change.means.new.treatment.episode=TRUE, # does a change in medication automatically start a new treatment episode?
                             dosage.change.means.new.treatment.episode=FALSE, # does a change in dosage automatically start a new treatment episode?
                             maximum.permissible.gap=180, # if a number, is the duration in units of max. permissible gaps between treatment episodes
                             maximum.permissible.gap.unit=c("days", "weeks", "months", "years", "percent")[1], # time units; can be "days", "weeks" (fixed at 7 days), "months" (fixed at 30 days), "years" (fixed at 365 days), or "percent", in which case maximum.permissible.gap is interpreted as a percent (can be > 100%) of the duration of the current prescription
                             maximum.permissible.gap.append.to.episode=FALSE, # should the maximum permissible gap be appended at the end of an episode with a gap larger than the maximum permissible gap? FALSE = no addition (the default), TRUE = the full maximum permissible gap is added
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
                             return.inner.event.info=FALSE, # should we return the event.info for all the individual events in each sliding window?
                             # Date format:
                             date.format="%m/%d/%Y", # the format of the dates used in this function (NA = undefined)
                             # Comments and metadata:
                             summary="CMA per treatment episode",
                             # The description of the output (added) columns:
                             event.interval.colname="event.interval", # contains number of days between the start of current event and the start of the next
                             gap.days.colname="gap.days", # contains the number of days when medication was not available
                             # Dealing with failed estimates:
                             force.NA.CMA.for.failed.patients=TRUE, # force the failed patients to have NA CM estimate?
                             # Return also the mapping between episodes and events?
                             return.mapping.events.episodes=FALSE, # if TRUE, return the mapping
                             # Parallel processing:
                             parallel.backend=c("none","multicore","snow","snow(SOCK)","snow(MPI)","snow(NWS)")[1], # parallel backend to use
                             parallel.threads="auto", # specification (or number) of parallel threads
                             # Misc:
                             suppress.warnings=FALSE,
                             suppress.special.argument.checks=TRUE, # used internally to suppress the check that we don't use special argument names
                             # extra parameters to be sent to the CMA function:
                             ...
)
{
  # Get the CMA function corresponding to the name:
  if( !(is.character(CMA.to.apply) || is.factor(CMA.to.apply)) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("'CMA.to.apply' must be a string contining the name of the simple CMA to apply!\n)"), "error", "CMA_per_episode", "AdhereR");
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
                    {if( !suppress.warnings ) .report.ewms(paste0("Unknown 'CMA.to.apply' '",CMA.to.apply,"': defaulting to CMA0!\n)"), "warning", "CMA_per_episode", "AdhereR"); CMA0;}); # by default, fall back to CMA0

  # Default argument values and overrides:
  def.vals <- formals(CMA.FNC);
  if( CMA.to.apply %in% c("CMA1", "CMA2", "CMA3", "CMA4") )
  {
    carryover.into.obs.window <- carryover.within.obs.window <- FALSE;
    if( !is.na(carry.only.for.same.medication) && carry.only.for.same.medication && !suppress.warnings ) .report.ewms("'carry.only.for.same.medication' cannot be defined for CMAs 1-4!\n", "warning", "CMA_per_episode", "AdhereR");
    carry.only.for.same.medication <- FALSE;
    if( !is.na(consider.dosage.change) && consider.dosage.change && !suppress.warnings ) .report.ewms("'consider.dosage.change' cannot be defined for CMAs 1-4!\n", "warning", "CMA_per_episode", "AdhereR");
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
    if( !suppress.warnings ) .report.ewms("I know how to do CMA per episodes only for CMAs 1 to 9!\n", "error", "CMA_per_episode", "AdhereR");
    return (NULL);
  }

  ## Force data to data.table
  #if( !inherits(data,"data.table") ) data <- as.data.table(data);

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

  ## retain only necessary columns of data
  #data <- data[,c(ID.colname,
  #                event.date.colname,
  #                event.duration.colname,
  #                event.daily.dose.colname,
  #                medication.class.colname),
  #             with=FALSE];

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
    mapping.episodes.to.events <- NULL; # the mapping events <-> episodes, if any
    if(is.null(treat.epi))
    {

      # Compute the treatment episodes:
      treat.epi <- compute.treatment.episodes(data=data,
                                              ID.colname=ID.colname,
                                              event.date.colname=event.date.colname,
                                              event.duration.colname=event.duration.colname,
                                              event.daily.dose.colname=event.daily.dose.colname,
                                              medication.class.colname=medication.class.colname,
                                              carryover.within.obs.window=carryover.within.obs.window,
                                              carry.only.for.same.medication=carry.only.for.same.medication,
                                              consider.dosage.change=consider.dosage.change,
                                              medication.change.means.new.treatment.episode=medication.change.means.new.treatment.episode,
                                              dosage.change.means.new.treatment.episode=dosage.change.means.new.treatment.episode,
                                              maximum.permissible.gap=maximum.permissible.gap,
                                              maximum.permissible.gap.unit=maximum.permissible.gap.unit,
                                              maximum.permissible.gap.append.to.episode=maximum.permissible.gap.append.to.episode,
                                              followup.window.start=followup.window.start,
                                              followup.window.start.unit=followup.window.start.unit,
                                              followup.window.duration=followup.window.duration,
                                              followup.window.duration.unit=followup.window.duration.unit,
                                              return.mapping.events.episodes=TRUE, # map the events to episodes anyways
                                              date.format=date.format,
                                              parallel.backend="none", # make sure this runs sequentially!
                                              parallel.threads=1,
                                              suppress.warnings=suppress.warnings,
                                              suppress.special.argument.checks=suppress.special.argument.checks,
                                              return.data.table=TRUE);
      if("mapping.episodes.to.events" %in% names(attributes(treat.epi))) mapping.episodes.to.events <- attr(treat.epi, "mapping.episodes.to.events"); # store the mapping events <-> episodes
    } else
    {

      # various checks

      # Convert treat.epi to data.table, cache event date as Date objects, and key by patient ID and event date
      treat.epi <- as.data.table(treat.epi);
      treat.epi[, `:=` (episode.start = as.Date(episode.start,format=date.format),
                        episode.end = as.Date(episode.end,format=date.format)
                        )]; # .DATE.as.Date: convert event.date.colname from formatted string to Date
      setkeyv(treat.epi, c(ID.colname, "episode.ID")); # key (and sorting) by patient and episode ID
    }

    if( is.null(treat.epi) || nrow(treat.epi) == 0 ) return (NULL);

    # Compute the real observation windows (might differ per patient) only once per patient (speed things up & the observation window is the same for all events within a patient):
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

    # Merge the observation window start and end dates back into the treatment episodes:
    treat.epi <- merge(treat.epi, event.info2[,c(ID.colname, ".OBS.START.DATE", ".OBS.END.DATE"),with=FALSE],
                       all.x=TRUE,
                       by = c(ID.colname));
    setnames(treat.epi, ncol(treat.epi)-c(1,0), c(".OBS.START.DATE.PRECOMPUTED", ".OBS.END.DATE.PRECOMPUTED"));
    # Get the intersection between the episode and the observation window:
    treat.epi[, c(".INTERSECT.EPISODE.OBS.WIN.START",
                  ".INTERSECT.EPISODE.OBS.WIN.END")
                := list(max(episode.start, .OBS.START.DATE.PRECOMPUTED),
                        min(episode.end,   .OBS.END.DATE.PRECOMPUTED)),
              by=c(ID.colname,"episode.ID")];
    treat.epi <- treat.epi[ .INTERSECT.EPISODE.OBS.WIN.START < .INTERSECT.EPISODE.OBS.WIN.END, ]; # keep only the episodes which fall within the OW
    treat.epi[, c("episode.duration",
                  ".INTERSECT.EPISODE.OBS.WIN.DURATION",
                  ".PATIENT.EPISODE.ID")
                := list(as.numeric(episode.end - episode.start),
                        as.numeric(.INTERSECT.EPISODE.OBS.WIN.END - .INTERSECT.EPISODE.OBS.WIN.START),
                        paste(get(ID.colname),episode.ID,sep="*"))];

    # Merge the data and the treatment episodes info:
    data.epi <- merge(treat.epi, data, allow.cartesian=TRUE);
    setkeyv(data.epi, c(".PATIENT.EPISODE.ID", ".DATE.as.Date"));

    # compute end.episode.gap.days, if treat.epi are supplied
    if(!"end.episode.gap.days" %in% colnames(treat.epi))
    {
      data.epi2 <- compute.event.int.gaps(data=as.data.frame(data.epi),
                                          ID.colname=".PATIENT.EPISODE.ID",
                                          event.date.colname=event.date.colname,
                                          event.duration.colname=event.duration.colname,
                                          event.daily.dose.colname=event.daily.dose.colname,
                                          medication.class.colname=medication.class.colname,
                                          carryover.within.obs.window=carryover.within.obs.window,
                                          carryover.into.obs.window=carryover.into.obs.window,
                                          carry.only.for.same.medication=carry.only.for.same.medication,
                                          consider.dosage.change=consider.dosage.change,
                                          followup.window.start="episode.start",
                                          followup.window.start.unit=followup.window.start.unit,
                                          followup.window.duration="episode.duration",
                                          followup.window.duration.unit=followup.window.duration.unit,
                                          observation.window.start=".INTERSECT.EPISODE.OBS.WIN.START",
                                          observation.window.duration=".INTERSECT.EPISODE.OBS.WIN.DURATION",
                                          observation.window.duration.unit="days",
                                          date.format=date.format,
                                          keep.window.start.end.dates=TRUE,
                                          remove.events.outside.followup.window=FALSE,
                                          parallel.backend="none", # make sure this runs sequentially!
                                          parallel.threads=1,
                                          suppress.warnings=suppress.warnings,
                                          suppress.special.argument.checks=suppress.special.argument.checks,
                                          return.data.table=TRUE);

      episode.gap.days <- data.epi2[which(.EVENT.WITHIN.FU.WINDOW), c(ID.colname, "episode.ID", gap.days.colname), by = c(ID.colname, "episode.ID"), with = FALSE]; # gap days during the follow-up window
      end.episode.gap.days <- episode.gap.days[,last(get(gap.days.colname)), by = c(ID.colname, "episode.ID")]; # gap days during the last event
      setnames(end.episode.gap.days, old = "V1", new = "end.episode.gap.days")

      treat.epi <- merge(treat.epi, end.episode.gap.days, all.x = TRUE, by = c(ID.colname, "episode.ID")); # merge end.episode.gap.days back to data.epi
      treat.epi[, episode.duration := as.numeric(.INTERSECT.EPISODE.OBS.WIN.END-.INTERSECT.EPISODE.OBS.WIN.START)];
    }

    # Compute the required CMA on this new combined database:
    if(length(dot.args <- list(...)) > 0 && "arguments.that.should.not.be.defined" %in% names(dot.args)) # check if arguments.that.should.not.be.defined is passed in the ... argument
    {
      # Avoid passing again arguments.that.should.not.be.defined to the CMA function:
      cma <- CMA.FNC(data=as.data.frame(data.epi),
                     ID.colname=".PATIENT.EPISODE.ID",
                     event.date.colname=event.date.colname,
                     event.duration.colname=event.duration.colname,
                     event.daily.dose.colname=event.daily.dose.colname,
                     medication.class.colname=medication.class.colname,
                     carryover.within.obs.window=carryover.within.obs.window,
                     carryover.into.obs.window=carryover.into.obs.window,
                     carry.only.for.same.medication=carry.only.for.same.medication,
                     consider.dosage.change=consider.dosage.change,
                     followup.window.start="episode.start",
                     followup.window.start.unit=followup.window.start.unit,
                     followup.window.duration="episode.duration",
                     followup.window.duration.unit=followup.window.duration.unit,
                     observation.window.start=".INTERSECT.EPISODE.OBS.WIN.START",
                     observation.window.duration=".INTERSECT.EPISODE.OBS.WIN.DURATION",
                     observation.window.duration.unit="days",
                     date.format=date.format,
                     parallel.backend="none", # make sure this runs sequentially!
                     parallel.threads=1,
                     suppress.warnings=suppress.warnings,
                     suppress.special.argument.checks=TRUE, # we know we're using special arguments, so don't do an extra check
                     ...);
    } else
    {
      # Temporarily avoid warnings linked to rewriting some CMA arguments by passing it arguments.that.should.not.be.defined=NULL:
      cma <- CMA.FNC(data=as.data.frame(data.epi),
                     ID.colname=".PATIENT.EPISODE.ID",
                     event.date.colname=event.date.colname,
                     event.duration.colname=event.duration.colname,
                     event.daily.dose.colname=event.daily.dose.colname,
                     medication.class.colname=medication.class.colname,
                     carryover.within.obs.window=carryover.within.obs.window,
                     carryover.into.obs.window=carryover.into.obs.window,
                     carry.only.for.same.medication=carry.only.for.same.medication,
                     consider.dosage.change=consider.dosage.change,
                     followup.window.start="episode.start",
                     followup.window.start.unit=followup.window.start.unit,
                     followup.window.duration="episode.duration",
                     followup.window.duration.unit=followup.window.duration.unit,
                     observation.window.start=".INTERSECT.EPISODE.OBS.WIN.START",
                     observation.window.duration=".INTERSECT.EPISODE.OBS.WIN.DURATION",
                     observation.window.duration.unit="days",
                     date.format=date.format,
                     parallel.backend="none", # make sure this runs sequentially!
                     parallel.threads=1,
                     suppress.warnings=suppress.warnings,
                     suppress.special.argument.checks=TRUE, # we know we're using special arguments, so don't do an extra check
                     arguments.that.should.not.be.defined=NULL,
                     ...);
    }

    # adjust episode start- and end dates
    treat.epi[, `:=` (episode.start = .INTERSECT.EPISODE.OBS.WIN.START,
                      episode.end = .INTERSECT.EPISODE.OBS.WIN.END)]

    # Add back the patient and episode IDs:
    tmp <- as.data.table(merge(cma$CMA, treat.epi)[,c(ID.colname, "episode.ID", "episode.start", "end.episode.gap.days", "episode.duration", "episode.end", "CMA")]);
    setkeyv(tmp, c(ID.colname,"episode.ID"));

    # The return value:
    ret.val <- list("CMA"=as.data.frame(tmp),
                    "event.info"=as.data.frame(event.info2)[,c(ID.colname, ".FU.START.DATE", ".FU.END.DATE", ".OBS.START.DATE", ".OBS.END.DATE")]);
    if( return.inner.event.info )
    {
      ret.val[["inner.event.info"]] <- as.data.frame(event.info2);
    }
    if( return.mapping.events.episodes )
    {
      events.actually.used <- cma$data$..ORIGINAL.ROW.ORDER..[ cma$event.info$..ORIGINAL.ROW.ORDER..[ cma$event.info$.EVENT.USED.IN.CMA ] ]; # translate back to the original data row numbers of the used events
      mapping.episodes.to.events <- mapping.episodes.to.events[ mapping.episodes.to.events$event.index.in.data %in% events.actually.used, ]; # keep only those events that are used in the episodes

      ret.val[["mapping.episodes.to.events"]] <- merge(mapping.episodes.to.events, tmp[,c(ID.colname, "episode.ID"), with=FALSE], by=c(ID.colname, "episode.ID"), all.x=FALSE, all.y=TRUE);
    }
    return (ret.val);
  }

  # Convert to data.table, cache event date as Date objects, and key by patient ID and event date
  # columns to keep:
  columns.to.keep <- c();
  if( !is.null(ID.colname) && !is.na(ID.colname) && length(ID.colname) == 1 && ID.colname %in% names(data) ) columns.to.keep <- c(columns.to.keep, ID.colname);
  if( !is.null(event.date.colname) && !is.na(event.date.colname) && length(event.date.colname) == 1 && event.date.colname %in% names(data) ) columns.to.keep <- c(columns.to.keep, event.date.colname);
  if( !is.null(event.duration.colname) && !is.na(event.duration.colname) && length(event.duration.colname) == 1  && event.duration.colname %in% names(data)) columns.to.keep <- c(columns.to.keep, event.duration.colname);
  if( !is.null(event.daily.dose.colname) && !is.na(event.daily.dose.colname) && length(event.daily.dose.colname) == 1 && event.daily.dose.colname %in% names(data) ) columns.to.keep <- c(columns.to.keep, event.daily.dose.colname);
  if( !is.null(medication.class.colname) && !is.na(medication.class.colname) && length(medication.class.colname) == 1 && medication.class.colname %in% names(data) ) columns.to.keep <- c(columns.to.keep, medication.class.colname);
  if( !is.null(mg <- getMGs(ret.val)) && !is.null(medication.groups.colname) && !is.na(medication.groups.colname) && length(medication.groups.colname) == 1  && medication.groups.colname %in% names(data)) columns.to.keep <- c(columns.to.keep, medication.groups.colname);
  if( !is.null(followup.window.start) && !is.na(followup.window.start) && length(followup.window.start) == 1 && (is.character(followup.window.start) || (is.factor(followup.window.start) && is.character(followup.window.start <- as.character(followup.window.start)))) && followup.window.start %in% names(data) ) columns.to.keep <- c(columns.to.keep, followup.window.start);
  if( !is.null(followup.window.duration) && !is.na(followup.window.duration) && length(followup.window.duration) == 1 && (is.character(followup.window.duration) || (is.factor(followup.window.duration) && is.character(followup.window.duration <- as.character(followup.window.duration)))) && followup.window.duration %in% names(data) ) columns.to.keep <- c(columns.to.keep, followup.window.duration);
  if( !is.null(observation.window.start) && !is.na(observation.window.start) && length(observation.window.start) == 1 && (is.character(observation.window.start) || (is.factor(observation.window.start) && is.character(observation.window.start <- as.character(observation.window.start)))) && observation.window.start %in% names(data) ) columns.to.keep <- c(columns.to.keep, observation.window.start);
  if( !is.null(observation.window.duration) && !is.na(observation.window.duration) && length(observation.window.duration) == 1 && (is.character(observation.window.duration) || (is.factor(observation.window.duration) && is.character(observation.window.duration <- as.character(observation.window.duration)))) && observation.window.duration %in% names(data) ) columns.to.keep <- c(columns.to.keep, observation.window.duration);
  # special column names that should not be used:
  if( !suppress.special.argument.checks )
  {
    ..special.colnames <- c(.special.colnames, event.interval.colname, gap.days.colname); # don't forget event.interval.colname and gap.days.colname!
    if( any(s <- ..special.colnames %in% columns.to.keep) )
    {
      .report.ewms(paste0("Column name(s) ",paste0("'",..special.colnames[s],"'",collapse=", "),"' are reserved: please don't use them in your input data!\n"), "error", cma.class.name, "AdhereR");
      return (NULL);
    }
  }
  # copy only the relevant bits of the data:
  data.copy <- data.table(data)[, columns.to.keep, with=FALSE];
  data.copy[, .DATE.as.Date := as.Date(get(event.date.colname),format=date.format)]; # .DATE.as.Date: convert event.date.colname from formatted string to Date
  data.copy$..ORIGINAL.ROW.ORDER.. <- 1:nrow(data.copy); # preserve the original order of the rows (needed for medication groups)
  setkeyv(data.copy, c(ID.colname, ".DATE.as.Date")); # key (and sorting) by patient ID and event date


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
    class(ret.val) <- "CMA_per_episode";
    ret.val$event.info <- as.data.frame(tmp$event.info);
    ret.val$computed.CMA <- CMA.to.apply;
    ret.val$summary <- summary;
    ret.val$CMA <- as.data.frame(tmp$CMA);
    setnames(ret.val$CMA, 1, ID.colname);
    if( return.inner.event.info && !is.null(tmp$inner.event.info) ) ret.val$inner.event.info <- as.data.frame(tmp$inner.event.info);
    if( return.mapping.events.episodes && !is.null(tmp$mapping.episodes.to.events) ) ret.val$mapping.episodes.to.events <- as.data.frame(tmp$mapping.episodes.to.events);

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
        return (list("CMA"=NULL, "event.info"=NULL));
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
      tmp$CMA <- as.data.frame(tmp$CMA); setnames(tmp$CMA, 1, ID.colname);
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
    ret.val$computed.CMA <- CMA.to.apply;
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
    class(ret.val) <- "CMA_per_episode";
    ret.val$summary <- summary;
    return (ret.val);

  }
}

#' @export
getMGs.CMA_per_episode <- function(x)
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA_per_episode") || is.null(cma$medication.groups) ) return (NULL);
  return (cma$medication.groups);
}

#' getEventsToEpisodesMapping
#'
#' This function returns the event-to-episode mapping, if this information exists.
#'
#' There are cases where it is interesting to know which events belong to which
#' episodes and which events have been actually used in computing the simple CMA
#' for each episode.
#' This information can be returned by \code{compute.treatment.episodes()} and
#' \code{CMA_per_episode()} if the parameter \code{return.mapping.events.episodes}
#' is set to \code{TRUE}.
#'
#' @param x Either a \code{data.frame} as returned by \code{compute.treatment.episodes()}
#' or an \code{CMA_per_episode object}.
#' @return The mapping between events and episodes, if it exists either as the
#' attribute \code{mapping.episodes.to.events} of the \code{data.frame} or as the
#' \code{mapping.episodes.to.events} component of the \code{CMA_per_episode object}
#' object, or \code{NULL} otherwise.
#' @export
getEventsToEpisodesMapping <- function(x)
{
  if( is.null(x) )
  {
    return (NULL);
  } else if( inherits(x, "CMA_per_episode") && ("mapping.episodes.to.events" %in% names(x)) && !is.null(x$mapping.episodes.to.events) && inherits(x$mapping.episodes.to.events, "data.frame") )
  {
    return (x$mapping.episodes.to.events);
  } else if( inherits(x, "data.frame") && ("mapping.episodes.to.events" %in% names(attributes(x))) && !is.null(attr(x,"mapping.episodes.to.events")) && inherits(attr(x,"mapping.episodes.to.events"), "data.frame") )
  {
    return (attr(x,"mapping.episodes.to.events"));
  } else
  {
    return (NULL);
  }
}

#' @export
getCMA.CMA_per_episode <- function(x, flatten.medication.groups=FALSE, medication.groups.colname=".MED_GROUP_ID")
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA_per_episode") || !("CMA" %in% names(cma)) || is.null(cma$CMA) ) return (NULL);
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
getEventInfo.CMA_per_episode <- function(x, flatten.medication.groups=FALSE, medication.groups.colname=".MED_GROUP_ID")
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA_per_episode") || !("event.info" %in% names(cma)) || is.null(cma$event.info) ) return (NULL);
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
getInnerEventInfo.CMA_per_episode <- function(x, flatten.medication.groups=FALSE, medication.groups.colname=".MED_GROUP_ID")
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) || !inherits(cma, "CMA_per_episode") || !("inner.event.info" %in% names(cma)) || is.null(cma$inner.event.info) ) return (NULL);
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
subsetCMA.CMA_per_episode <- function(cma, patients, suppress.warnings=FALSE)
{
  if( inherits(patients, "factor") ) patients <- as.character(patients);
  all.patients <- unique(cma$data[,cma$ID.colname]);
  patients.to.keep <- intersect(patients, all.patients);
  if( length(patients.to.keep) == length(all.patients) )
  {
    # Keep all patients:
    return (cma);
  }
  if( length(patients.to.keep) == 0 )
  {
    if( !suppress.warnings ) .report.ewms("No patients to subset on!\n", "error", "subsetCMA.CMA_per_episode", "AdhereR");
    return (NULL);
  }
  if( length(patients.to.keep) < length(patients) && !suppress.warnings ) .report.ewms("Some patients in the subsetting set are not in the CMA itself and are ignored!\n", "warning", "subsetCMA.CMA_per_episode", "AdhereR");

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
print.CMA_per_episode <- function(x,                                     # the CMA_per_episode (or derived) object
                                  ...,                                   # required for S3 consistency
                                  inline=FALSE,                          # print inside a line of text or as a separate, extended object?
                                  format=c("text", "latex", "markdown"), # the format to print to
                                  print.params=TRUE,                     # show the parameters?
                                  print.data=TRUE,                       # show the summary of the data?
                                  exclude.params=c("event.info", "inner.event.info", "mapping.episodes.to.events"), # if so, should I not print some?
                                  skip.header=FALSE,                     # should I print the generic header?
                                  cma.type=class(x)[1]
)
{
  cma <- x; # parameter x is required for S3 consistency, but I like cma more
  if( is.null(cma) ) return (invisible(NULL));

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
    .report.ewms("Unknown format for printing!\n", "error", "print.CMA_per_episode", "AdhereR");
    return (invisible(NULL));
  }
}


#' Plot CMA_per_episode and CMA_sliding_window objects.
#'
#' Plots the event data and the estimated CMA per treatment episode and sliding
#' window, respectively.
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
#' sequence of events, and carry no information about the avilability of
#' medicine in this interval.
#'
#' Above these, the treatment episodes or the sliding windows are represented in
#' a stacked manner from the earlieast (left, bottom of the stack) to the latest
#' (right, top of the stack), each showing the CMA as percent fill (capped at
#' 100\% even if CMA values may be higher) and also as text.
#'
#' The follow-up and the observation windows are plotted as empty an rectangle
#' and as shaded rectangle, respectively (for some CMAs the observation window
#' might be adjusted in which case the adjustment may also be plotted using a
#' different shading).
#'
#' The kernel density ("smoothed histogram") of the CMA estimates across
#' treatment episodes/sliding windows (if more than 2) can be visually
#' represented as well in the left side of the figure (NB, their horizontal
#' scales may be different across patients).
#'
#' When several patients are displayed on the same plot, they are organized
#' vertically, and alternating bands (white and gray) help distinguish
#' consecutive patients.
#' Implicitely, all patients contained in the \code{cma} object will be plotted,
#' but the \code{patients.to.plot} parameter allows the selection of a subset
#' of patients.
#'
#' Finally, the y-axis shows the patient ID and possibly the CMA estimate as
#' well.
#'
#' Any not explicitely defined arguments are passed to the simple CMA estimation
#' and plotting function; therefore, for more info about possible estimation
#' parameters plese see the help for the appropriate simple CMA, and for possible
#' aesthetic tweaks, please see the help for their plotting.
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
#' @param legend.cex,legend.cex.title The legend and legend title font sizes.
#' @param cex,cex.axis,cex.lab \emph{numeric} values specifying the cex of the
#' various types of text.
#' @param show.cma \emph{Logical}, should the CMA type be shown in the title?
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
#' @param show.event.intervals \emph{Logical}, should the actual event intervals
#' be shown? As per-episode and sliding windows might have overlapping intervals,
#' it is better not to show them by default (\code{FALSE}).
#' @param show.overlapping.event.intervals specifies how to plot the event
#' intervals that appear in multiple sliding windows or episodes. We can plot
#' how thye look in the \emph{first} sliding window or episode (the default),
#' how they appear in the \emph{last}, pick the one that minimizes the gap
#' (\emph{min gap}) or maximizes it (\emph{max gap}), or compute their
#' \emph{average} across all sliding windows or episodes containing them.
#' @param plot.events.vertically.displaced Should consecutive events be plotted
#' on separate rows (i.e., separated vertically, the default) or on the same row?
#' @param print.dose,cex.dose,print.dose.outline.col,print.dose.centered Print daily
#' dose as a number and, if so, how (color, size, position...).
#' @param plot.dose,lwd.event.max.dose,plot.dose.lwd.across.medication.classes
#' Show dose through the width of the event lines and, if so, what the maximum
#' width should be, and should this maximum be by medication class or overall.
#' @param col.na The colour used for missing event data.
#' @param col.continuation,lty.continuation,lwd.continuation The color, style
#' and width of the contuniation lines connecting consecutive events.
#' @param alternating.bands.cols The colors of the alternating vertical bands
#' distinguishing the patients; can be \code{NULL} = don't draw the bandes;
#' or a vector of colors.
#' @param print.episode.or.sliding.window \emph{Logical}, should we show which
#' events belong to which episode or sliding window? To work, the CMA must have
#' been constructed with \code{return.mapping.events.episodes} or
#' \code{return.mapping.events.sliding.window} set to \code{TRUE}, respectively.
#' @param bw.plot \emph{Logical}, should the plot use grayscale only (i.e., the
#' \code{\link[grDevices]{gray.colors}} function)?
#' @param rotate.text \emph{Numeric}, the angle by which certain text elements
#' (e.g., axis labels) should be rotated.
#' @param force.draw.text \emph{Logical}, if \code{TRUE}, always draw text even
#' if too big or too small
#' @param print.CMA \emph{Logical}, should the CMA values be printed?
#' @param CMA.cex ... and, if printed, what cex (\emph{numeric}) to use?
#' @param plot.CMA \emph{Logical}, should the distribution of the CMA values
#' across episodes/sliding windows be plotted? If \code{TRUE} (the default), the
#' distribution is shown on the left-hand side of the plot, otherwise it is not.
#' @param plot.CMA.as.histogram \emph{Logical}, should the CMA plot be a
#' histogram or a (truncated) density plot? Please note that it is TRUE by
#' deafult for CMA_per_episode and FALSE for CMA_sliding_window, because
#' usually there are more sliding windows than episodes. Also, the density
#' estimate cannot be estimated for less than three different values.
#' @param plot.partial.CMAs.as Should the partial CMAs be plotted? Possible values
#' are "stacked", "overlapping" or "timeseries", or \code{NULL} for no partial
#' CMA plots. Please note that \code{plot.CMA} and \code{plot.partial.CMAs.as}
#' are independent of each other.
#' @param plot.partial.CMAs.as.stacked.col.bars,plot.partial.CMAs.as.stacked.col.border,plot.partial.CMAs.as.stacked.col.text
#' If plotting the partial CMAs as stacked bars, define their graphical attributes.
#' @param plot.partial.CMAs.as.timeseries.vspace,plot.partial.CMAs.as.timeseries.start.from.zero,plot.partial.CMAs.as.timeseries.col.dot,plot.partial.CMAs.as.timeseries.col.interval,plot.partial.CMAs.as.timeseries.col.text,plot.partial.CMAs.as.timeseries.interval.type,plot.partial.CMAs.as.timeseries.lwd.interval,plot.partial.CMAs.as.timeseries.alpha.interval,plot.partial.CMAs.as.timeseries.show.0perc,plot.partial.CMAs.as.timeseries.show.100perc
#' If plotting the partial CMAs as imeseries, these are their graphical attributes.
#' @param plot.partial.CMAs.as.overlapping.alternate,plot.partial.CMAs.as.overlapping.col.interval,plot.partial.CMAs.as.overlapping.col.text
#' If plotting the partial CMAs as overlapping segments, these are their
#' graphical attributes.
#' @param CMA.plot.ratio A \emph{number}, the proportion of the total horizontal
#' plot space to be allocated to the CMA plot.
#' @param CMA.plot.col,CMA.plot.border,CMA.plot.bkg,CMA.plot.text \emph{Strings}
#' giving the colours of the various components of the CMA plot.
#' @param highlight.followup.window \emph{Logical}, should the follow-up window
#' be plotted?
#' @param followup.window.col The follow-up window colour.
#' @param highlight.observation.window \emph{Logical}, should the observation
#' window be plotted?
#' @param observation.window.col,observation.window.opacity
#' Attributes of the observation window (colour, transparency).
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
#' @param suppress.warnings \emph{Logical}, if \code{TRUE} don't show any
#' warnings.
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
#' @param ... other parameters (to be passed to the estimation and plotting of
#' the simple CMA)
#'
#' @seealso See the simple CMA estimation \code{\link[AdhereR]{CMA1}}
#' to \code{\link[AdhereR]{CMA9}} and plotting \code{\link[AdhereR]{plot.CMA1}}
#' functions for extra parameters.
#'
#' @examples
#' \dontrun{
#' cmaW <- CMA_sliding_window(CMA=CMA1,
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
#'                         sliding.window.start=0,
#'                         sliding.window.start.unit="days",
#'                         sliding.window.duration=90,
#'                         sliding.window.duration.unit="days",
#'                         sliding.window.step.duration=7,
#'                         sliding.window.step.unit="days",
#'                         sliding.window.no.steps=NA,
#'                         date.format="%m/%d/%Y"
#'                        );
#' plot(cmaW, patients.to.plot=c("1","2"));
#' cmaE <- CMA_per_episode(CMA=CMA1,
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
#' plot(cmaE, patients.to.plot=c("1","2"));}
#' @export
plot.CMA_per_episode <- function(x,                                     # the CMA_per_episode or CMA_sliding_window (or derived) object
                                 patients.to.plot=NULL,                 # list of patient IDs to plot or NULL for all
                                 duration=NA,                           # duration and end period to plot in days (if missing, determined from the data)
                                 align.all.patients=FALSE, align.first.event.at.zero=FALSE, # should all patients be aligned? and, if so, place the first event as the horizintal 0?
                                 show.period=c("dates","days")[2],      # draw vertical bars at regular interval as dates or days?
                                 period.in.days=90,                     # the interval (in days) at which to draw veritcal lines
                                 show.legend=TRUE, legend.x="right", legend.y="bottom", legend.bkg.opacity=0.5, legend.cex=0.75, legend.cex.title=1.0, # legend params and position
                                 cex=1.0, cex.axis=0.75, cex.lab=1.0,   # various graphical params
                                 show.cma=TRUE,                         # show the CMA type
                                 xlab=c("dates"="Date", "days"="Days"), # Vector of x labels to show for the two types of periods, or a single value for both, or NULL for nothing
                                 ylab=c("withoutCMA"="patient", "withCMA"="patient (& CMA)"), # Vector of y labels to show without and with CMA estimates, or a single value for both, or NULL ofr nonthing
                                 title=c("aligned"="Event patterns (all patients aligned)", "notaligned"="Event patterns"), # Vector of titles to show for and without alignment, or a single value for both, or NULL for nonthing
                                 col.cats=rainbow,                      # single color or a function mapping the categories to colors
                                 unspecified.category.label="drug",     # the label of the unspecified category of medication
                                 medication.groups.to.plot=NULL,        # the names of the medication groups to plot (by default, all)
                                 medication.groups.separator.show=TRUE, medication.groups.separator.lty="solid", medication.groups.separator.lwd=2, medication.groups.separator.color="blue", # group medication events by patient?
                                 medication.groups.allother.label="*",  # the label to use for the __ALL_OTHERS__ medication class (defaults to *)
                                 lty.event="solid", lwd.event=2, pch.start.event=15, pch.end.event=16, # event style
                                 show.event.intervals=FALSE,            # per-episode and sliding windows might have overlapping intervals, so better not to show them by default
                                 show.overlapping.event.intervals=c("first", "last", "min gap", "max gap", "average")[1], # how to plot overlapping event intervals (relevant for sliding windows and per episode)
                                 plot.events.vertically.displaced=TRUE, # display the events on different lines (vertical displacement) or not (defaults to TRUE)?
                                 print.dose=FALSE, cex.dose=0.75, print.dose.outline.col="white", print.dose.centered=FALSE, # print daily dose
                                 plot.dose=FALSE, lwd.event.max.dose=8, plot.dose.lwd.across.medication.classes=FALSE, # draw daily dose as line width
                                 col.na="lightgray",                    # color for mising data
                                 col.continuation="black", lty.continuation="dotted", lwd.continuation=1, # style of the contuniation lines connecting consecutive events
                                 print.CMA=TRUE, CMA.cex=0.50,    # print CMA next to the participant's ID?
                                 plot.CMA=TRUE,                   # plot the CMA next to the participant ID?
                                 plot.CMA.as.histogram=TRUE,      # plot CMA as a histogram or as a density plot?
                                 plot.partial.CMAs.as=c("stacked", "overlapping", "timeseries")[1], # how to plot the "partial" (i.e., intervals/episodes) CMAs (NULL for none)?
                                 plot.partial.CMAs.as.stacked.col.bars="gray90", plot.partial.CMAs.as.stacked.col.border="gray30", plot.partial.CMAs.as.stacked.col.text="black",
                                 plot.partial.CMAs.as.timeseries.vspace=7, # how much vertical space to reserve for the timeseries plot (in character lines)
                                 plot.partial.CMAs.as.timeseries.start.from.zero=TRUE, #show the vertical axis start at 0 or at the minimum actual value (if positive)?
                                 plot.partial.CMAs.as.timeseries.col.dot="darkblue", plot.partial.CMAs.as.timeseries.col.interval="gray70", plot.partial.CMAs.as.timeseries.col.text="firebrick", # setting any of these to NA results in them not being plotted
                                 plot.partial.CMAs.as.timeseries.interval.type=c("none", "segments", "arrows", "lines", "rectangles")[2], # how to show the covered intervals
                                 plot.partial.CMAs.as.timeseries.lwd.interval=1, # line width for some types of intervals
                                 plot.partial.CMAs.as.timeseries.alpha.interval=0.25, # the transparency of the intervales (when drawn as rectangles)
                                 plot.partial.CMAs.as.timeseries.show.0perc=TRUE, plot.partial.CMAs.as.timeseries.show.100perc=FALSE, #show the 0% and 100% lines?
                                 plot.partial.CMAs.as.overlapping.alternate=TRUE, # should successive intervals be plotted low/high?
                                 plot.partial.CMAs.as.overlapping.col.interval="gray70", plot.partial.CMAs.as.overlapping.col.text="firebrick", # setting any of these to NA results in them not being plotted
                                 CMA.plot.ratio=0.10,             # the proportion of the total horizontal plot to be taken by the CMA plot
                                 CMA.plot.col="lightgreen", CMA.plot.border="darkgreen", CMA.plot.bkg="aquamarine", CMA.plot.text=CMA.plot.border, # attributes of the CMA plot
                                 highlight.followup.window=TRUE, followup.window.col="green",
                                 highlight.observation.window=TRUE, observation.window.col="yellow", observation.window.opacity=0.3,
                                 print.episode.or.sliding.window=FALSE, # should we print the episode or sliding window to which an event belongs?
                                 alternating.bands.cols=c("white", "gray95"), # the colors of the alternating vertical bands across patients (NULL=don't draw any; can be >= 1 color)
                                 bw.plot=FALSE,                   # if TRUE, override all user-given colors and replace them with a scheme suitable for grayscale plotting
                                 rotate.text=-60,                 # some text (e.g., axis labels) may be rotated by this much degrees
                                 force.draw.text=FALSE,           # if true, always draw text even if too big or too small
                                 min.plot.size.in.characters.horiz=0, min.plot.size.in.characters.vert=0, # the minimum plot size (in characters: horizontally, for the whole duration, vertically, per event (and, if shown, per episode/sliding window))
                                 max.patients.to.plot=100,        # maximum number of patients to plot
                                 export.formats=NULL,                   # the formats to export the figure to (by default, none); can be any subset of "svg" (just SVG file), "html" (SVG + HTML + CSS + JavaScript all embedded within the HTML document), "jpg", "png", "webp", "ps" and "pdf"
                                 export.formats.fileprefix="AdhereR-plot", # the file name prefix for the exported formats
                                 export.formats.height=NA, export.formats.width=NA, # desired dimensions (in pixels) for the exported figure (defaults to sane values)
                                 export.formats.save.svg.placeholder=TRUE,
                                 export.formats.svg.placeholder.type=c("jpg", "png", "webp")[2],
                                 export.formats.svg.placeholder.embed=FALSE, # save a placeholder for the SVG image?
                                 export.formats.html.template=NULL, export.formats.html.javascript=NULL, export.formats.html.css=NULL, # HTML, JavaScript and CSS templates for exporting HTML+SVG
                                 export.formats.directory=NA,           # if exporting, which directory to export to (if not give, creates files in the temporary directory)
                                 generate.R.plot=TRUE,                  # generate standard (base R) plot for plotting within R?
                                 do.not.draw.plot=FALSE,                # if TRUE, don't draw the actual plot, but only the legend (if required)
                                 suppress.warnings=FALSE,         # suppress warnings?
                                 ...
)
{
  #if( show.event.intervals )
  #{
  #  if( !suppress.warnings ) .report.ewms("show.event.intervals=TRUE is not yet implemented in plotting sliding windows and per episodes!\n", "message", "plot", "AdhereR");
  #  show.event.intervals <- FALSE;
  #}

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
             show.cma=show.cma,
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
             show.event.intervals=show.event.intervals,
             show.overlapping.event.intervals=show.overlapping.event.intervals,
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
             print.CMA=print.CMA,
             CMA.cex=CMA.cex,
             plot.CMA=plot.CMA,
             plot.CMA.as.histogram=plot.CMA.as.histogram,
             CMA.plot.ratio=CMA.plot.ratio,
             CMA.plot.col=CMA.plot.col,
             CMA.plot.border=CMA.plot.border,
             CMA.plot.bkg=CMA.plot.bkg,
             CMA.plot.text=CMA.plot.text,
             plot.partial.CMAs.as=plot.partial.CMAs.as,
             plot.partial.CMAs.as.stacked.col.bars=plot.partial.CMAs.as.stacked.col.bars,
             plot.partial.CMAs.as.stacked.col.border=plot.partial.CMAs.as.stacked.col.border,
             plot.partial.CMAs.as.stacked.col.text=plot.partial.CMAs.as.stacked.col.text,
             plot.partial.CMAs.as.timeseries.vspace=plot.partial.CMAs.as.timeseries.vspace,
             plot.partial.CMAs.as.timeseries.start.from.zero=plot.partial.CMAs.as.timeseries.start.from.zero,
             plot.partial.CMAs.as.timeseries.col.dot=plot.partial.CMAs.as.timeseries.col.dot,
             plot.partial.CMAs.as.timeseries.col.interval=plot.partial.CMAs.as.timeseries.col.interval,
             plot.partial.CMAs.as.timeseries.col.text=plot.partial.CMAs.as.timeseries.col.text,
             plot.partial.CMAs.as.timeseries.interval.type=plot.partial.CMAs.as.timeseries.interval.type,
             plot.partial.CMAs.as.timeseries.lwd.interval=plot.partial.CMAs.as.timeseries.lwd.interval,
             plot.partial.CMAs.as.timeseries.alpha.interval=plot.partial.CMAs.as.timeseries.alpha.interval,
             plot.partial.CMAs.as.timeseries.show.0perc=plot.partial.CMAs.as.timeseries.show.0perc,
             plot.partial.CMAs.as.timeseries.show.100perc=plot.partial.CMAs.as.timeseries.show.100perc,
             plot.partial.CMAs.as.overlapping.alternate=plot.partial.CMAs.as.overlapping.alternate,
             plot.partial.CMAs.as.overlapping.col.interval=plot.partial.CMAs.as.overlapping.col.interval,
             plot.partial.CMAs.as.overlapping.col.text=plot.partial.CMAs.as.overlapping.col.text,
             highlight.followup.window=highlight.followup.window,
             followup.window.col=followup.window.col,
             highlight.observation.window=highlight.observation.window,
             observation.window.col=observation.window.col,
             observation.window.opacity=observation.window.opacity,
             print.episode.or.sliding.window=print.episode.or.sliding.window,
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