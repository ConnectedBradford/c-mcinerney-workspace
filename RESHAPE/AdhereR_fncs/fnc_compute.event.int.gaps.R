#' Gap Days and Event (prescribing or dispensing) Intervals.
#'
#' For a given event (prescribing or dispensing) database, compute the gap days
#' and event intervals in various scenarious.
#'
#' This should in general not be called directly by the user, but is provided as
#' a basis for the extension to new CMAs.
#'
#' @param data A \emph{\code{data.frame}} containing the events used to
#' compute the CMA. Must contain, at a minimum, the patient unique ID, the event
#' date and duration, and might also contain the daily dosage and medication
#' type (the actual column names are defined in the following four parameters);
#' the \code{CMA} constructors call this parameter \code{data}.
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
#' \code{data} containing the classes/types/groups of medication, or \code{NA}
#' if not defined.
#' @param event.interval.colname A \emph{string}, the name of a newly-created
#' column storing the number of days between the start of the current event and
#' the start of the next one; the default value "event.interval" should be
#' changed only if there is a naming conflict with a pre-existing
#' "event.interval" column in \code{event.info}.
#' @param gap.days.colname A \emph{string}, the name of a newly-created column
#' storing the number of days when medication was not available (i.e., the
#' "gap days"); the default value "gap.days" should be changed only if there is
#' a naming conflict with a pre-existing "gap.days" column in \code{event.info}.
#' @param carryover.within.obs.window \emph{Logical}, if \code{TRUE} consider
#' the carry-over within the observation window, or \code{NA} if not defined.
#' @param carryover.into.obs.window \emph{Logical}, if \code{TRUE} consider the
#' carry-over from before the starting date of the observation window, or
#' \code{NA} if not defined.
#' @param carry.only.for.same.medication \emph{Logical}, if \code{TRUE} the
#' carry-over applies only across medication of the same type, or \code{NA}
#' if not defined.
#' @param consider.dosage.change \emph{Logical}, if \code{TRUE} the carry-over
#' is adjusted to reflect changes in dosage, or \code{NA} if not defined.
#' @param followup.window.start If a \emph{\code{Date}} object, it represents
#' the actual start date of the follow-up window; if a \emph{string} it is the
#' name of the column in \code{data} containing the start date of the follow-up
#' window either as the numbers of \code{followup.window.start.unit} units after
#' the first event (the column must be of type \code{numeric}) or as actual
#' dates (in which case the column must be of type \code{Date}); if a
#' \emph{number} it is the number of time units defined in the
#' \code{followup.window.start.unit} parameter after the begin of the
#' participant's first event.
#' @param followup.window.start.unit can be either \emph{"days"},
#' \emph{"weeks"}, \emph{"months"} or \emph{"years"}, and represents the time
#' units that \code{followup.window.start} refers to (when a number), or
#' \code{NA} if not defined.
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
#' @param keep.window.start.end.dates \emph{Logical}, should the computed start
#' and end dates of the windows be kept?
#' @param keep.event.interval.for.all.events \emph{Logical}, should the computed
#' event intervals be kept for all events, or \code{NA}'ed for those outside the
#' OW?
#' @param remove.events.outside.followup.window \emph{Logical}, should the
#' events that fall outside the follo-wup window be removed from the results?
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
#' @param return.data.table \emph{Logical}, if \code{TRUE} return a
#' \code{data.table} object, otherwise a \code{data.frame}.
#' @param ... extra arguments.
#' @return A \code{data.frame} or \code{data.table} extending the
#' \code{event.info} parameter with:
#' \itemize{
#'  \item \code{event.interval} Or any other name given in
#'  \code{event.interval.colname}, containing the number of days between the
#'  start of the current event and the start of the next one.
#'  \item \code{gap.days} Or any other name given in \code{gap.days.colname},
#'  containing the number of days when medication was not available for the
#'  current event (i.e., the "gap days").
#'  \item \code{.FU.START.DATE,.FU.END.DATE} if kept, the actual start and end
#'  dates of the follow-up window (after adjustments due to the various
#'  parameters).
#'  \item \code{.OBS.START.DATE,.OBS.END.DATE} if kept, the actual start and end
#'  dates of the observation window (after adjustments due to the various
#'  parameters).
#'  \item \code{.EVENT.STARTS.BEFORE.OBS.WINDOW} if kept, \code{TRUE} if the
#'  current event starts before the start of the observation window.
#'  \item \code{.TDIFF1,.TDIFF2} if kept, various auxiliary time differences
#'  (in days).
#'  \item \code{.EVENT.STARTS.AFTER.OBS.WINDOW} if kept, \code{TRUE} if the
#'  current event starts after the end of the observation window.
#'  \item \code{.CARRY.OVER.FROM.BEFORE} if kept, the carry-over (if any) from
#'  the previous events.
#'  \item \code{.EVENT.WITHIN.FU.WINDOW} if kept, \code{TRUE} if the current
#'  event is within the follow-up window.
#' }
#' @export
compute.event.int.gaps <- function(data, # this is a per-event data.frame with columns:
                                   ID.colname=NA, # the name of the column containing the unique patient ID
                                   event.date.colname=NA, # the start date of the event in the date.format format
                                   event.duration.colname=NA, # the event duration in days
                                   event.daily.dose.colname=NA, # the prescribed daily dose
                                   medication.class.colname=NA, # the classes/types/groups of medication
                                   # The description of the output (added) columns:
                                   event.interval.colname="event.interval", # contains number of days between the start of current event and the start of the next
                                   gap.days.colname="gap.days", # contains the number of days when medication was not available
                                   # Various types methods of computing gaps:
                                   carryover.within.obs.window=FALSE, # if TRUE consider the carry-over within the observation window
                                   carryover.into.obs.window=FALSE, # if TRUE consider the carry-over from before the starting date of the observation window
                                   carry.only.for.same.medication=FALSE, # if TRUE the carry-over applies only across medication of same type
                                   consider.dosage.change=FALSE, # if TRUE carry-over is adjusted to reflect changes in dosage
                                   # The follow-up window:
                                   followup.window.start=0, # if a number, date earliest event per participant + number of units, otherwise a date.format date or variable date
                                   followup.window.start.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!)
                                   followup.window.duration=365*2, # the duration of the follow-up window in the time units given below
                                   followup.window.duration.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!)
                                   # The observation window (embedded in the follow-up window):
                                   observation.window.start=0, # the number of time units relative to followup.window.start, otherwise a date.format date or variable date
                                   observation.window.start.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!)
                                   observation.window.duration=365*2, # the duration of the observation window in time units
                                   observation.window.duration.unit=c("days", "weeks", "months", "years")[1], # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!)
                                   # Date format:
                                   date.format="%m/%d/%Y", # the format of the dates used in this function
                                   # Keep the follow-up and observation window start and end dates?
                                   keep.window.start.end.dates=FALSE,
                                   remove.events.outside.followup.window=TRUE, # remove events outside the follow-up window?
                                   keep.event.interval.for.all.events=FALSE, # keep the event.interval estimates for all events
                                   # Parallel processing:
                                   parallel.backend=c("none","multicore","snow","snow(SOCK)","snow(MPI)","snow(NWS)")[1], # parallel backend to use
                                   parallel.threads="auto", # specification (or number) of parallel threads
                                   # Misc:
                                   suppress.warnings=FALSE,
                                   suppress.special.argument.checks=FALSE, # used internally to suppress the check that we don't use special argument names
                                   return.data.table=FALSE,  # should the result be converted to data.frame (default) or returned as a data.table (keyed by patient ID and event date)?
                                   ... # other stuff
)
{
  # preconditions concerning column names:
  if( is.null(data) || !inherits(data,"data.frame") || nrow(data) < 1 )
  {
    if( !suppress.warnings ) .report.ewms("Event data must be a non-empty data frame!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  data.names <- names(data); # cache names(data) as it is used a lot
  if( is.null(ID.colname) || is.na(ID.colname) ||                                           # avoid empty stuff
      !(is.character(ID.colname) ||                                                         # it must be a character...
        (is.factor(ID.colname) && is.character(ID.colname <- as.character(ID.colname)))) || # ...or a factor (forced to character)
      length(ID.colname) != 1 ||                                                            # make sure it's a single value
      !(ID.colname %in% data.names)                                                         # make sure it's a valid column name
      )
  {
    if( !suppress.warnings ) .report.ewms(paste0("The patient ID column \"",ID.colname,"\" cannot be empty, must be a single value, and must be present in the event data!\n"), "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( is.null(event.date.colname) || is.na(event.date.colname) ||                                                   # avoid empty stuff
      !(is.character(event.date.colname) ||                                                                         # it must be a character...
        (is.factor(event.date.colname) && is.character(event.date.colname <- as.character(event.date.colname)))) || # ...or a factor (forced to character)
      length(event.date.colname) != 1 ||                                                                            # make sure it's a single value
      !(event.date.colname %in% data.names)                                                                         # make sure it's a valid column name
      )
  {
    if( !suppress.warnings ) .report.ewms(paste0("The event date column \"",event.date.colname,"\" cannot be empty, must be a single value, and must be present in the event data!\n"), "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( is.null(event.duration.colname) || is.na(event.duration.colname) ||                                                       # avoid empty stuff
      !(is.character(event.duration.colname) ||                                                                                 # it must be a character...
        (is.factor(event.duration.colname) && is.character(event.duration.colname <- as.character(event.duration.colname)))) || # ...or a factor (forced to character)
      length(event.duration.colname) != 1 ||                                                                                    # make sure it's a single value
      !(event.duration.colname %in% data.names)                                                                                 # make sure it's a valid column name
      )
  {
    if( !suppress.warnings ) .report.ewms(paste0("The event duration column \"",event.duration.colname,"\" cannot be empty, must be a single value, and must be present in the event data!\n"), "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( (!is.null(event.daily.dose.colname) && !is.na(event.daily.dose.colname)) &&                                                      # if actually given:
      (!(is.character(event.daily.dose.colname) ||                                                                                     # it must be a character...
         (is.factor(event.daily.dose.colname) && is.character(event.daily.dose.colname <- as.character(event.daily.dose.colname)))) || # ...or a factor (forced to character)
       length(event.daily.dose.colname) != 1 ||                                                                                        # make sure it's a single value
       !(event.daily.dose.colname %in% data.names))                                                                                    # make sure it's a valid column name
      )
  {
    if( !suppress.warnings ) .report.ewms(paste0("If given, the event daily dose column \"",event.daily.dose.colname,"\" must be a single value and must be present in the event data!\n"), "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( (!is.null(medication.class.colname) && !is.na(medication.class.colname)) &&                                                      # if actually given:
      (!(is.character(medication.class.colname) ||                                                                                     # it must be a character...
         (is.factor(medication.class.colname) && is.character(medication.class.colname <- as.character(medication.class.colname)))) || # ...or a factor (forced to character)
       length(medication.class.colname) != 1 ||                                                                                        # make sure it's a single value
       !(medication.class.colname %in% data.names))                                                                                    # make sure it's a valid column name
      )
  {
    if( !suppress.warnings ) .report.ewms(paste0("If given, the event type column \"",medication.class.colname,"\" must be a single value and must be present in the event data!\n"), "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }

  # preconditions concerning carry-over:
  if( !is.logical(carryover.within.obs.window)    || is.na(carryover.within.obs.window)    || length(carryover.within.obs.window) != 1    ||
      !is.logical(carryover.into.obs.window)      || is.na(carryover.into.obs.window)      || length(carryover.into.obs.window) != 1      ||
      !is.logical(carry.only.for.same.medication) || is.na(carry.only.for.same.medication) || length(carry.only.for.same.medication) != 1 )
  {
    if( !suppress.warnings ) .report.ewms("Carry over arguments must be single value logicals!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( !carryover.within.obs.window && !carryover.into.obs.window && carry.only.for.same.medication )
  {
    if( !suppress.warnings ) .report.ewms("Cannot carry over only for same medication when no carry over at all is considered!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }

  # preconditions concerning dosage change:
  if( !is.logical(consider.dosage.change) || is.na(consider.dosage.change) || length(consider.dosage.change) != 1 )
  {
    if( !suppress.warnings ) .report.ewms("Consider dosage change must be single value logical!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }

  # preconditions concerning follow-up window (as all violations result in the same error, aggregate them in a single if):
  if( (is.null(followup.window.start) || is.na(followup.window.start) || length(followup.window.start) != 1) ||                   # cannot be missing or have more than one values
      (!inherits(followup.window.start,"Date") && !is.numeric(followup.window.start) &&                                           # not a Date or number:
          (!(is.character(followup.window.start) ||                                                                               # it must be a character...
             (is.factor(followup.window.start) && is.character(followup.window.start <- as.character(followup.window.start)))) || # ...or a factor (forced to character)
           !(followup.window.start %in% data.names))) )                                                                           # make sure it's a valid column name
  {
    if( !suppress.warnings ) .report.ewms("The follow-up window start must be a single value, either a number, a Date object, or a string giving a column name in the data!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( is.null(followup.window.start.unit) ||
      (is.na(followup.window.start.unit) && !(is.factor(followup.window.start) || is.character(followup.window.start))) ||
      length(followup.window.start.unit) != 1 ||
      ((is.factor(followup.window.start.unit) || is.character(followup.window.start.unit)) && !(followup.window.start.unit %in% c("days", "weeks", "months", "years"))) )
  {
    if( !suppress.warnings ) .report.ewms("The follow-up window start unit must be a single value, one of \"days\", \"weeks\", \"months\" or \"years\"!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( is.numeric(followup.window.duration) && (followup.window.duration <= 0 || length(followup.window.duration) != 1) ||               # cannot be missing or have more than one values
      (!is.numeric(followup.window.duration) &&
       (!(is.character(followup.window.duration) ||                                                                                     # it must be a character...
          (is.factor(followup.window.duration) && is.character(followup.window.duration <- as.character(followup.window.duration)))) || # ...or a factor (forced to character)
        !(followup.window.duration %in% data.names))))                                                                                  # make sure it's a valid column name
  {
    if( !suppress.warnings ) .report.ewms("The follow-up window duration must be a single value, either a positive number, or a string giving a column name in the data!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( is.null(followup.window.duration.unit) || is.na(followup.window.duration.unit) ||
      length(followup.window.duration.unit) != 1 ||
      !(followup.window.duration.unit %in% c("days", "weeks", "months", "years") ) )
  {
    if( !suppress.warnings ) .report.ewms("The follow-up window duration unit must be a single value, one of \"days\", \"weeks\", \"months\" or \"years\"!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }

  # preconditions concerning observation window (as all violations result in the same error, aggregate them in a single if):
  if( (is.null(observation.window.start) || is.na(observation.window.start) || length(observation.window.start) != 1) ||                   # cannot be missing or have more than one values
      (is.numeric(observation.window.start) && (observation.window.start < 0)) ||                                                          # if a number, must be a single positive one
      (!inherits(observation.window.start,"Date") && !is.numeric(observation.window.start) &&                                              # not a Date or number:
          (!(is.character(observation.window.start) ||                                                                                     # it must be a character...
             (is.factor(observation.window.start) && is.character(observation.window.start <- as.character(observation.window.start)))) || # ...or a factor (forced to character)
           !(observation.window.start %in% data.names))) )                                                                                 # make sure it's a valid column name
  {
    if( !suppress.warnings ) .report.ewms("The observation window start must be a single value, either a positive number, a Date object, or a string giving a column name in the data!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( is.null(observation.window.start.unit) ||
      (is.na(observation.window.start.unit) && !(is.factor(observation.window.start) || is.character(observation.window.start))) ||
      length(observation.window.start.unit) != 1 ||
      ((is.factor(observation.window.start.unit) || is.character(observation.window.start.unit)) && !(observation.window.start.unit %in% c("days", "weeks", "months", "years"))) )
  {
    if( !suppress.warnings ) .report.ewms("The observation window start unit must be a single value, one of \"days\", \"weeks\", \"months\" or \"years\"!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( is.numeric(observation.window.duration) && (observation.window.duration <= 0 || length(observation.window.duration) != 1) ||               # cannot be missing or have more than one values
      (!is.numeric(observation.window.duration) &&
       (!(is.character(observation.window.duration) ||                                                                                           # it must be a character...
          (is.factor(observation.window.duration) && is.character(observation.window.duration <- as.character(observation.window.duration)))) || # ...or a factor (forced to character)
        !(observation.window.duration %in% data.names))))                                                                                        # make sure it's a valid column name
  {
    if( !suppress.warnings ) .report.ewms("The observation window duration must be a single value, either a positive number, or a string giving a column name in the data!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }
  if( is.null(observation.window.duration.unit) || is.na(observation.window.duration.unit) ||
      length(observation.window.duration.unit) != 1 ||
      !(observation.window.duration.unit %in% c("days", "weeks", "months", "years") ) )
  {
    if( !suppress.warnings ) .report.ewms("The observation window duration unit must be a single value, one of \"days\", \"weeks\", \"months\" or \"years\"!\n", "error", "compute.event.int.gaps", "AdhereR")
    return (NULL);
  }

  # Check the patient IDs:
  if( anyNA(data[,ID.colname]) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("The patient unique identifiers in the \"",ID.colname,"\" column must not contain NAs; the first occurs on row ",min(which(is.na(data[,ID.colname]))),"!\n"), "error", "compute.event.int.gaps", "AdhereR");
    return (NULL);
  }

  # Check the date format (and save the conversion to Date() for later use):
  if( is.na(date.format) || is.null(date.format) || length(date.format) != 1 || !is.character(date.format) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("The date format must be a single string!\n"), "error", "compute.event.int.gaps", "AdhereR");
    return (NULL);
  }
  if( anyNA(Date.converted.to.DATE <- as.Date(data[,event.date.colname],format=date.format)) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("Not all entries in the event date \"",event.date.colname,"\" column are valid dates or conform to the date format \"",date.format,"\"; first issue occurs on row ",min(which(is.na(Date.converted.to.DATE))),"!\n"), "error", "compute.event.int.gaps", "AdhereR");
    return (NULL);
  }

  # Check the duration:
  tmp <- data[,event.duration.colname]; # caching for speed
  if( !is.numeric(tmp) || any(is.na(tmp) | tmp <= 0) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("The event durations in the \"",event.duration.colname,"\" column must be non-missing strictly positive numbers!\n"), "error", "compute.event.int.gaps", "AdhereR");
    return (NULL);
  }

  # Check the event daily dose:
  if( !is.na(event.daily.dose.colname) && !is.null(event.daily.dose.colname) &&             # if actually given:
      (!is.numeric(tmp <- data[,event.daily.dose.colname]) || any(is.na(tmp) | tmp <= 0)) ) # must be a non-missing strictly positive number (and cache it for speed)
  {
    if( !suppress.warnings ) .report.ewms(paste0("If given, the event daily dose in the \"",event.daily.dose.colname,"\" column must be a non-missing strictly positive numbers!\n"), "error", "compute.event.int.gaps", "AdhereR");
    return (NULL);
  }

  # Check the newly created columns:
  if( is.na(event.interval.colname) || is.null(event.interval.colname) || !is.character(event.interval.colname) || (event.interval.colname %in% data.names) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("The column name where the event interval will be stored \"",event.interval.colname,"\" cannot be missing nor already present in the event data!\n"), "error", "compute.event.int.gaps", "AdhereR");
    return (NULL);
  }
  if( is.na(gap.days.colname) || is.null(gap.days.colname) || !is.character(gap.days.colname) || (gap.days.colname %in% data.names) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("The column name where the gap days will be stored \"",gap.days.colname,"\" cannot be mising nor already present in the event data.\n"), "error", "compute.event.int.gaps", "AdhereR");
    return (NULL);
  }

  # Make a data.table copy of data so we can alter it without altering the original input data:
  ret.val <- data.table(data);

  event.date2.colname <- ".DATE.as.Date"; # name of column caching the event dates
  ret.val[,c(event.interval.colname,            # output column event.interval.colname
             gap.days.colname,                  # output column gap.days.colname
             ".DATE.as.Date",                   # cache column .DATE.as.Date
             ".FU.START.DATE",                  # cache FUW start date
             ".FU.END.DATE",                    # cache FUW end date
             ".OBS.START.DATE",                 # cache OW start date
             ".OBS.END.DATE",                   # cache OW end date
             ".OBS.WITHIN.FU",                  # cache if the OW falls within the FUW
             ".EVENT.STARTS.BEFORE.OBS.WINDOW", # cache if the event starts before the OW begins
             ".EVENT.STARTS.AFTER.OBS.WINDOW",  # cache if the event starts after the OW ends
             ".EVENT.WITHIN.FU.WINDOW",         # cache if the event is within the FUW
             ".TDIFF1",                         # cache time difference betwen the event duration and (OW start - event start)
             ".TDIFF2",                         # cache time difference betwen the next event start date and the current event start date
             ".CARRY.OVER.FROM.BEFORE"          # cache if there is carry over from before the current event
            ) :=
            list(NA_real_,    # event.interval.colname: initially NA (numeric)
                 NA_real_,    # gap.days.colname: initially NA (numeric)
                 Date.converted.to.DATE, # .DATE.as.Date: convert event.date.colname from formatted string to Date (speeding calendar computations)
                 as.Date(NA), # .FU.START.DATE: initially NA (of type Date)
                 as.Date(NA), # .FU.END.DATE: initially NA (of type Date)
                 as.Date(NA), # .OBS.START.DATE: initially NA (of type Date)
                 as.Date(NA), # .OBS.END.DATE: initially NA (of type Date)
                 NA,          # .OBS.WITHIN.FU: initially NA (logical)
                 NA,          # .EVENT.STARTS.BEFORE.OBS.WINDOW: initially NA (logical)
                 NA,          # .EVENT.STARTS.AFTER.OBS.WINDOW: initially NA (logical)
                 NA,          # .EVENT.WITHIN.FU.WINDOW: initially NA (logical)
                 NA_real_,    # .TDIFF1: initially NA (numeric)
                 NA_real_,    # .TDIFF2: initially NA (numeric)
                 NA_real_     # .CARRY.OVER.FROM.BEFORE: initially NA (numeric)
            )];

  # Cache types of followup.window.start and observation.window.start (1 = numeric, 2 = column with Dates, 3 = column with units, 4 = a Date), as well as durations:
  followup.window.start.type <- ifelse(is.numeric(followup.window.start),
                                       1, # a number in the appropriate units
                                       ifelse(followup.window.start %in% data.names,
                                              ifelse(inherits(data[,followup.window.start],"Date"),
                                                     2,  # column of Date objects
                                                     3), # column of numbers in the appropriate units
                                              4)); # a Date object
  followup.window.duration.is.number <- is.numeric(followup.window.duration)
  observation.window.start.type <- ifelse(is.numeric(observation.window.start),
                                          1, # a number in the appropriate units
                                          ifelse(observation.window.start %in% data.names,
                                                 ifelse(inherits(data[,observation.window.start],"Date"),
                                                        2,  # column of Date objects
                                                        3), # column of numbers in the appropriate units
                                                 4)); # a Date object
  observation.window.duration.is.number <- is.numeric(observation.window.duration)

  setkeyv(ret.val, c(ID.colname, ".DATE.as.Date")); # key (and sorting) by patient ID and event date

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
    # Auxiliary internal function: For a given patient, compute the gaps and return the required columns:
    .process.patient <- function(data4ID)
    {
      # Number of events:
      n.events <- nrow(data4ID);

      # Force the selection, evaluation of promises and caching of the needed columns:
      # ... which columns to select (with their indices):
      columns.to.cache <- c(event.date2.colname, event.duration.colname); event.date2.colname.index <- 1; event.duration.colname.index <- 2;
      curr.index <- 3;
      if( !is.na(medication.class.colname) ){ columns.to.cache <- c(columns.to.cache, medication.class.colname); medication.class.colname.index <- curr.index; curr.index <- curr.index + 1;}
      if( !is.na(event.daily.dose.colname) ){ columns.to.cache <- c(columns.to.cache, event.daily.dose.colname); event.daily.dose.colname.index <- curr.index; curr.index <- curr.index + 1;}
      if( followup.window.start.type %in% c(2,3) ){ columns.to.cache <- c(columns.to.cache, followup.window.start); followup.window.start.index <- curr.index; curr.index <- curr.index + 1;}
      if( !followup.window.duration.is.number ){ columns.to.cache <- c(columns.to.cache, followup.window.duration); followup.window.duration.index <- curr.index; curr.index <- curr.index + 1;}
      if( observation.window.start.type %in% c(2,3) ){ columns.to.cache <- c(columns.to.cache, observation.window.start); observation.window.start.index <- curr.index; curr.index <- curr.index + 1;}
      if( !observation.window.duration.is.number ){ columns.to.cache <- c(columns.to.cache, observation.window.duration); observation.window.duration.index <- curr.index; curr.index <- curr.index + 1;}
      # ... select these columns:
      data4ID.selected.columns <- data4ID[, columns.to.cache, with=FALSE]; # alternative to: data4ID[,..columns.to.cache]
      # ... cache the columns based on their indices:
      event.date2.column <- data4ID.selected.columns[[event.date2.colname.index]]; event.duration.column <- data4ID.selected.columns[[event.duration.colname.index]];
      if( !is.na(medication.class.colname) ) medication.class.column <- data4ID.selected.columns[[medication.class.colname.index]];
      if( !is.na(event.daily.dose.colname) ) event.daily.dose.column <- data4ID.selected.columns[[event.daily.dose.colname.index]];
      if( followup.window.start.type %in% c(2,3) ) followup.window.start.column <- data4ID.selected.columns[[followup.window.start.index]];
      if( !followup.window.duration.is.number ) followup.window.duration.column <- data4ID.selected.columns[[followup.window.duration.index]];
      if( observation.window.start.type %in% c(2,3) ) observation.window.start.column <- data4ID.selected.columns[[observation.window.start.index]];
      if( !observation.window.duration.is.number ) observation.window.duration.column <- data4ID.selected.columns[[observation.window.duration.index]];

      # Cache also follow-up window start and end dates:
      # start dates
      .FU.START.DATE <- switch(followup.window.start.type,
                               .add.time.interval.to.date(event.date2.column[1], followup.window.start, followup.window.start.unit, suppress.warnings),                 # 1
                               followup.window.start.column[1],                                                                                                         # 2
                               .add.time.interval.to.date(event.date2.column[1], followup.window.start.column[1], followup.window.start.unit, suppress.warnings),       # 3
                               followup.window.start);                                                                                                                  # 4
      if( is.na(.FU.START.DATE) )
      {
        # Return a valid but empty object:
        return (list(".FU.START.DATE"=as.Date(NA),
                     ".FU.END.DATE"=as.Date(NA),
                     ".OBS.START.DATE"=as.Date(NA),
                     ".OBS.END.DATE"=as.Date(NA),
                     ".EVENT.STARTS.BEFORE.OBS.WINDOW"=NA,
                     ".EVENT.STARTS.AFTER.OBS.WINDOW"=NA,
                     ".EVENT.WITHIN.FU.WINDOW"=NA,
                     ".TDIFF1"=NA_real_,
                     ".TDIFF2"=NA_real_,
                     ".OBS.WITHIN.FU"=FALSE,
                     ".EVENT.INTERVAL"=NA_real_,
                     ".GAP.DAYS"=NA_real_,
                     ".CARRY.OVER.FROM.BEFORE"=NA_real_));
      }
      # end dates
      .FU.END.DATE <- .add.time.interval.to.date(.FU.START.DATE,
                                                 ifelse(followup.window.duration.is.number, followup.window.duration, followup.window.duration.column[1]),
                                                 followup.window.duration.unit,
                                                 suppress.warnings);
      if( is.na(.FU.END.DATE) )
      {
        # Return a valid but empty object:
        return (list(".FU.START.DATE"=.FU.START.DATE,
                     ".FU.END.DATE"=as.Date(NA),
                     ".OBS.START.DATE"=as.Date(NA),
                     ".OBS.END.DATE"=as.Date(NA),
                     ".EVENT.STARTS.BEFORE.OBS.WINDOW"=NA,
                     ".EVENT.STARTS.AFTER.OBS.WINDOW"=NA,
                     ".EVENT.WITHIN.FU.WINDOW"=NA,
                     ".TDIFF1"=NA_real_,
                     ".TDIFF2"=NA_real_,
                     ".OBS.WITHIN.FU"=FALSE,
                     ".EVENT.INTERVAL"=NA_real_,
                     ".GAP.DAYS"=NA_real_,
                     ".CARRY.OVER.FROM.BEFORE"=NA_real_));
      }
      # Is the event within the FU window?
      .EVENT.WITHIN.FU.WINDOW <- (event.date2.column >= .FU.START.DATE) & (event.date2.column <= .FU.END.DATE);

      # Cache also observation window start and end dates:
      # start dates
      .OBS.START.DATE <- switch(observation.window.start.type,
                                .add.time.interval.to.date(.FU.START.DATE, observation.window.start, observation.window.start.unit, suppress.warnings),                 # 1
                                observation.window.start.column[1],                                                                                                     # 2
                                .add.time.interval.to.date(.FU.START.DATE, observation.window.start.column[1], observation.window.start.unit, suppress.warnings),       # 3
                                observation.window.start);                                                                                                              # 4
      if( is.na(.OBS.START.DATE) )
      {
        # Return a valid but empty object:
        return (list(".FU.START.DATE"=.FU.START.DATE,
                     ".FU.END.DATE"=.FU.END.DATE,
                     ".OBS.START.DATE"=as.Date(NA),
                     ".OBS.END.DATE"=as.Date(NA),
                     ".EVENT.STARTS.BEFORE.OBS.WINDOW"=NA,
                     ".EVENT.STARTS.AFTER.OBS.WINDOW"=NA,
                     ".EVENT.WITHIN.FU.WINDOW"=.EVENT.WITHIN.FU.WINDOW,
                     ".TDIFF1"=NA_real_,
                     ".TDIFF2"=NA_real_,
                     ".OBS.WITHIN.FU"=FALSE,
                     ".EVENT.INTERVAL"=NA_real_,
                     ".GAP.DAYS"=NA_real_,
                     ".CARRY.OVER.FROM.BEFORE"=NA_real_));
      }
      # end dates
      .OBS.END.DATE <- .add.time.interval.to.date(.OBS.START.DATE,
                                                  ifelse(observation.window.duration.is.number, observation.window.duration, observation.window.duration.column[1]),
                                                  observation.window.duration.unit,
                                                  suppress.warnings);
      if( is.na(.OBS.END.DATE) )
      {
        # Return a valid but empty object:
        return (list(".FU.START.DATE"=.FU.START.DATE,
                     ".FU.END.DATE"=.FU.END.DATE,
                     ".OBS.START.DATE"=.OBS.START.DATE,
                     ".OBS.END.DATE"=as.Date(NA),
                     ".EVENT.STARTS.BEFORE.OBS.WINDOW"=NA,
                     ".EVENT.STARTS.AFTER.OBS.WINDOW"=NA,
                     ".EVENT.WITHIN.FU.WINDOW"=.EVENT.WITHIN.FU.WINDOW,
                     ".TDIFF1"=NA_real_,
                     ".TDIFF2"=NA_real_,
                     ".OBS.WITHIN.FU"=FALSE,
                     ".EVENT.INTERVAL"=NA_real_,
                     ".GAP.DAYS"=NA_real_,
                     ".CARRY.OVER.FROM.BEFORE"=NA_real_));
      }

      # For each event, starts before/after the observation window start date?
      .EVENT.STARTS.BEFORE.OBS.WINDOW <- (event.date2.column < .OBS.START.DATE);
      .EVENT.STARTS.AFTER.OBS.WINDOW  <- (event.date2.column > .OBS.END.DATE);

      # Cache some time differences:
      # event.duration.colname - (.OBS.START.DATE - event.date2.colname):
      .TDIFF1 <- (event.duration.column - .difftime.Dates.as.days(.OBS.START.DATE, event.date2.column));
      # event.date2.colname[current+1] - event.date2.colname[current]:
      if( n.events > 1 )
      {
        .TDIFF2 <- c(.difftime.Dates.as.days(event.date2.column[-1], event.date2.column[-n.events]),NA_real_);
      } else
      {
        .TDIFF2 <- NA_real_;
      }

      # Make sure the observation window is included in the follow-up window:
      if( (.FU.START.DATE > .OBS.START.DATE) || (.FU.END.DATE < .OBS.END.DATE) )
      {
        # Return a valid but empty object:
        return (list(".FU.START.DATE"=.FU.START.DATE,
                     ".FU.END.DATE"=.FU.END.DATE,
                     ".OBS.START.DATE"=.OBS.START.DATE,
                     ".OBS.END.DATE"=.OBS.END.DATE,
                     ".EVENT.STARTS.BEFORE.OBS.WINDOW"=.EVENT.STARTS.BEFORE.OBS.WINDOW,
                     ".EVENT.STARTS.AFTER.OBS.WINDOW"=.EVENT.STARTS.AFTER.OBS.WINDOW,
                     ".EVENT.WITHIN.FU.WINDOW"=.EVENT.WITHIN.FU.WINDOW,
                     ".TDIFF1"=.TDIFF1,
                     ".TDIFF2"=.TDIFF2,
                     ".OBS.WITHIN.FU"=FALSE,
                     ".EVENT.INTERVAL"=NA_real_,
                     ".GAP.DAYS"=NA_real_,
                     ".CARRY.OVER.FROM.BEFORE"=NA_real_));
      } else
      {
        .OBS.WITHIN.FU <- TRUE;
      }

      .CARRY.OVER.FROM.BEFORE <- .EVENT.INTERVAL <- .GAP.DAYS <- rep(NA_real_,n.events); # initialize to NA
      # Select only those events within the FUW and OW (depending on carryover into OW):
      if( !carryover.into.obs.window )
      {
        s <- which(.EVENT.WITHIN.FU.WINDOW & !.EVENT.STARTS.BEFORE.OBS.WINDOW & !.EVENT.STARTS.AFTER.OBS.WINDOW); # select all events for this patient within the observation window
      } else
      {
        s <- which(.EVENT.WITHIN.FU.WINDOW & !.EVENT.STARTS.AFTER.OBS.WINDOW); # select all events for patient ID within the follow-up window
      }
      slen <- length(s); # cache it up
      if( slen == 1 ) # only one event in the observation window
      {
        # Computations happen within the observation window
        .EVENT.INTERVAL[s] <- .difftime.Dates.as.days(.OBS.END.DATE, event.date2.column[s]); # for last event, the interval ends at the end of OW
        .CARRY.OVER.FROM.BEFORE[s] <- 0.0; # no carry-over into this unique event
        .GAP.DAYS[s] <- max(0.0, (.EVENT.INTERVAL[s] - event.duration.column[s])); # the actual gap cannot be negative
      } else if( slen > 1 ) # at least one event in the observation window
      {
        # Computations happen within the observation window
        if( carryover.within.obs.window )
        {
          # do carry over within the observation window:

          # was there a change in medication?
          if( !is.na(medication.class.colname) )
          {
            medication.changed <- c((medication.class.column[s[-slen]] != medication.class.column[s[-1]]), FALSE);
          } else
          {
            medication.changed <- rep(FALSE,slen);
          }
          # was there a change in dosage?
          if( !is.na(event.daily.dose.colname) )
          {
            dosage.change.ratio <- c((event.daily.dose.column[s[-slen]] / event.daily.dose.column[s[-1]]), 1.0);
          } else
          {
            dosage.change.ratio <- rep(1.0,slen);
          }
          # event intervals:
          .EVENT.INTERVAL[s]       <- .TDIFF2[s]; # save the time difference between next and current event start dates, but
          .EVENT.INTERVAL[s[slen]] <- .difftime.Dates.as.days(.OBS.END.DATE, event.date2.column[s[slen]]); # for last event, the interval ends at the end of OW
          # event.interval - event.duration:
          .event.interval.minus.duration <- (.EVENT.INTERVAL[s] - event.duration.column[s]);
          # cache various medication and dosage change conditions:
          cond1 <- (carry.only.for.same.medication & medication.changed);
          cond2 <- ((carry.only.for.same.medication & consider.dosage.change & !medication.changed) | (!carry.only.for.same.medication & consider.dosage.change));

          carry.over <- 0; # Initially, no carry-over into the first event
          # for each event:
          for( i in seq_along(s) ) # this for loop is not a performance bottleneck!
          {
            si <- s[i]; # caching s[i] as it's used a lot

            # Save the carry-over into the event:
            .CARRY.OVER.FROM.BEFORE[si] <- carry.over;

            # Computing the gap between events:
            gap <- (.event.interval.minus.duration[i] - carry.over); # subtract the event duration and carry-over from event interval
            if( gap < 0.0 ) # the actual gap cannot be negative
            {
              .GAP.DAYS[si] <- 0.0; carry.over <- (-gap);
            } else
            {
              .GAP.DAYS[si] <- gap; carry.over <- 0.0;
            }

            if( cond1[i] )
            {
              # Do not carry over across medication changes:
              carry.over <- 0;
            } else if( cond2[i] )
            {
              # Apply the dosage change ratio:
              carry.over <- carry.over * dosage.change.ratio[i]; # adjust the carry-over relative to dosage change
            }
          }
        } else
        {
          # do not consider carry over within the observation window:

          # event intervals:
          .EVENT.INTERVAL[s]       <- .TDIFF2[s]; # save the time difference between next and current event start dates, but
          .EVENT.INTERVAL[s[slen]] <- .difftime.Dates.as.days(.OBS.END.DATE, event.date2.column[s[slen]]); # for last event, the interval ends at the end of OW
          # event.interval - event.duration:
          .event.interval.minus.duration <- (.EVENT.INTERVAL[s] - event.duration.column[s]);
          # no carry over in this case:
          .CARRY.OVER.FROM.BEFORE[s] <- 0.0;

          # for each event:
          for( i in seq_along(s) ) # this for loop is not a performance bottleneck!
          {
            si <- s[i]; # caching s[i] as it's used a lot

            # Computing the gap between events:
            gap <- .event.interval.minus.duration[i]; # the event duration
            if( gap < 0.0 ) # the actual gap cannot be negative
            {
              .GAP.DAYS[si] <- 0.0;
            } else
            {
              .GAP.DAYS[si] <- gap;
            }
          }
        }
      }

      # Return the computed columns:
      return (list(".FU.START.DATE"=.FU.START.DATE,
                   ".FU.END.DATE"=.FU.END.DATE,
                   ".OBS.START.DATE"=.OBS.START.DATE,
                   ".OBS.END.DATE"=.OBS.END.DATE,
                   ".EVENT.STARTS.BEFORE.OBS.WINDOW"=.EVENT.STARTS.BEFORE.OBS.WINDOW,
                   ".EVENT.STARTS.AFTER.OBS.WINDOW"=.EVENT.STARTS.AFTER.OBS.WINDOW,
                   ".EVENT.WITHIN.FU.WINDOW"=.EVENT.WITHIN.FU.WINDOW,
                   ".TDIFF1"=.TDIFF1,
                   ".TDIFF2"=.TDIFF2,
                   ".OBS.WITHIN.FU"=.OBS.WITHIN.FU,
                   ".EVENT.INTERVAL"=.EVENT.INTERVAL,
                   ".GAP.DAYS"=.GAP.DAYS,
                   ".CARRY.OVER.FROM.BEFORE"=.CARRY.OVER.FROM.BEFORE));
    }

    # Don't attempt to proces an empty dataset:
    if( is.null(data) || nrow(data) == 0 ) return (NULL);

    data[, c(".FU.START.DATE", ".FU.END.DATE",
             ".OBS.START.DATE", ".OBS.END.DATE",
             ".EVENT.STARTS.BEFORE.OBS.WINDOW", ".EVENT.STARTS.AFTER.OBS.WINDOW", ".EVENT.WITHIN.FU.WINDOW",
             ".TDIFF1", ".TDIFF2",
             ".OBS.WITHIN.FU",
             event.interval.colname, gap.days.colname,
             ".CARRY.OVER.FROM.BEFORE") :=
           .process.patient(.SD), # for each patient, compute the various columns and assign them back into the data.table
           by=ID.colname # group by patients
        ];
    return (data);
  }

  # Compute the workhorse function:
  ret.val <- .compute.function(.workhorse.function,
                               parallel.backend=parallel.backend,
                               parallel.threads=parallel.threads,
                               data=ret.val,
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

  if( is.null(ret.val) || nrow(ret.val) < 1 )
  {
    if( !suppress.warnings ) .report.ewms("Computing event intervals and gap days failed!\n", "error", "compute.event.int.gaps", "AdhereR");
    return (NULL);
  }
  if( any(!ret.val$.OBS.WITHIN.FU) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("The observation window is not within the follow-up window for participant(s) ",paste0(unique(ret.val[!ret.val$.OBS.WITHIN.FU,get(ID.colname)]),collapse=", ")," !\n"), "error", "compute.event.int.gaps", "AdhereR");
    return (NULL);
  }

  # Make sure events outside the observation window are NA'ed:
  if( !keep.event.interval.for.all.events ) ret.val[ .OBS.WITHIN.FU & (.EVENT.STARTS.BEFORE.OBS.WINDOW | .EVENT.STARTS.AFTER.OBS.WINDOW), c(event.interval.colname) := NA_real_ ];

  # Remove the events that fall outside the follow-up window:
  if( remove.events.outside.followup.window ) ret.val <- ret.val[ which(.OBS.WITHIN.FU & .EVENT.WITHIN.FU.WINDOW), ];
  # If the results are empty return NULL:
  if( is.null(ret.val) || nrow(ret.val) < 1 )
  {
    if( !suppress.warnings ) .report.ewms("No observations fall within the follow-up and observation windows!\n", "error", "compute.event.int.gaps", "AdhereR");
    return (NULL);
  }

  # Final clean-up: delete temporary columns:
  if( !keep.window.start.end.dates )
  {
    ret.val[,c(".DATE.as.Date",
               ".FU.START.DATE",
               ".FU.END.DATE",
               ".OBS.START.DATE",
               ".OBS.END.DATE",
               ".OBS.WITHIN.FU",
               ".EVENT.STARTS.BEFORE.OBS.WINDOW",
               ".TDIFF1",
               ".TDIFF2",
               ".EVENT.STARTS.AFTER.OBS.WINDOW",
               ".CARRY.OVER.FROM.BEFORE",
               ".EVENT.WITHIN.FU.WINDOW") := NULL];
  }

  # If the results are empty return NULL:
  if( is.null(ret.val) || nrow(ret.val) < 1 ) return (NULL);

  if( !return.data.table )
  {
    return (as.data.frame(ret.val));
  } else
  {
    if( ".DATE.as.Date" %in% names(ret.val) )
    {
      setkeyv(ret.val, c(ID.colname, ".DATE.as.Date")); # make sure it is keyed by patient ID and event date
    } else
    {
      setkeyv(ret.val, c(ID.colname)); # make sure it is keyed by patient ID (as event date was removed)
    }
    return (ret.val);
  }
}