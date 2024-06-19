compute_event_durations <- function(disp.data = NULL,
                                    presc.data = NULL,
                                    special.periods.data = NULL,
                                    ID.colname,
                                    medication.class.colnames,
                                    disp.date.colname,
                                    total.dose.colname,
                                    presc.date.colname,
                                    presc.daily.dose.colname,
                                    presc.duration.colname,
                                    visit.colname,
                                    split.on.dosage.change = TRUE,
                                    force.init.presc = FALSE,
                                    force.presc.renew = FALSE,
                                    trt.interruption = c("continue", "discard", "carryover")[1],
                                    special.periods.method = trt.interruption,
                                    carryover = FALSE,
                                    date.format = "%d.%m.%Y",
                                    suppress.warnings = FALSE,
                                    return.data.table = FALSE,
                                    progress.bar = TRUE,
                                    ...)
{

  # set carryover to false
  # carryover <- FALSE # remove when carryover argument is properly implemented

  # Preconditions:
  {
    # dispensing data class and dimensions:
    if( inherits(disp.data, "matrix") ) disp.data <- as.data.table(disp.data); # convert matrix to data.table
    if( !inherits(disp.data, "data.frame") )
    {
      if( !suppress.warnings ) warning("The dispensing data must be of type 'data.frame'!\n");
      return (NULL);
    }
    if( nrow(disp.data) < 1 )
    {
      if( !suppress.warnings ) warning("The dispensing data must have at least one row!\n");
      return (NULL);
    }
    # prescribing data class and dimensions:
    if( inherits(presc.data, "matrix") ) presc.data <- as.data.table(presc.data); # convert matrix to data.table
    if( !inherits(presc.data, "data.frame") )
    {
      if( !suppress.warnings ) warning("The prescribing data must be of type 'data.frame'!\n");
      return (NULL);
    }
    if( nrow(presc.data) < 1 )
    {
      if( !suppress.warnings ) warning("The prescribing data must have at least one row!\n");
      return (NULL);
    }
    # special period data class and dimensions:
    if(!is.null(special.periods.data))
    {
      special.periods.data <- data.table(special.periods.data)

      if( inherits(special.periods.data, "matrix") ) special.periods.data <- as.data.table(special.periods.data); # convert matrix to data.table
      if( !inherits(special.periods.data, "data.frame") )
      {
        if( !suppress.warnings ) warning("The special periods data must be of type 'data.frame'!\n");
        return (NULL);
      }
      if( nrow(special.periods.data) < 1 )
      {
        if( !suppress.warnings ) warning("The special periods data must have at least one row!\n");
        return (NULL);
      }
      if(!all(c(ID.colname, "DATE.IN", "DATE.OUT") %in% colnames(special.periods.data)))
      {
        if( !suppress.warnings ) warning(paste0("The special periods data must contain at least all
                                                columns with the names '", ID.colname, "', 'DATE.IN', and 'DATE.OUT'.\n
                                                Please refer to the documentation for more information.\n"));
        return (NULL);
      }

      # if(!all(colnames(special.periods.data) %in% c(ID.colname, "DATE.IN", "DATE.OUT", "TYPE", special.periods.method, medication.class.colnames)))
      # {
      #   if( !suppress.warnings ) warning(paste0("The special periods data can only contain columns
      #                                           with the names \"", ID.colname, "\", \"DATE.IN\", \"DATE.OUT\", \"TYPE\", ",
      #                                           paste(shQuote(medication.class.colnames), collapse = ", "), ", and a column with
      #                                           customized instructions how to handle a specific period.\n
      #                                           Please refer to the documentation for more information.\n"));
      #   return (NULL);
      # }

      if( !special.periods.method %in% c("continue", "discard", "carryover") && !special.periods.method %in% names(special.periods.data))
      {
        if( !suppress.warnings ) warning(paste0("special.periods.method must be either of 'continue', 'discard',
                                                'carryover', or a column name in the special periods data!\n"));
        return (NULL);
      }
      if(special.periods.method %in% names(special.periods.data) && any(!unique(special.periods.data[[special.periods.method]] %in% c("continue", "discard", "carryover"))))
      {
        unexpected.values <- unique(special.periods.data[[special.periods.method]][!special.periods.data[[special.periods.method]] %in% c("continue", "discard", "carryover")])

        if( !suppress.warnings ) warning(paste0("Column special.periods.method='",special.periods.method, "' in special periods data contains unexpected values: ",
                                                unexpected.values,"\n"));
        return (NULL);
      }
    }

    # the column names must exist in dispensing and prescription data:
    if( !is.na(ID.colname) && !(ID.colname %in% names(disp.data)) && !(ID.colname %in% names(presc.data)))
    {
      if( !suppress.warnings ) warning(paste0("Column ID.colname='",ID.colname,"' must appear in the dispensing and prescribing data!\n"));
      return (NULL);
    }
    if( !is.na(presc.date.colname) && !(presc.date.colname %in% names(presc.data)) )
    {
      if( !suppress.warnings ) warning(paste0("Column presc.date.colname='",presc.date.colname,"' must appear in the prescribing data!\n"));
      return (NULL);
    }
    if(anyNA(presc.data[[presc.date.colname]])){
      if( !suppress.warnings ) warning(paste0("Column presc.date.colname='",presc.date.colname,"' cannot contain missing values!\n"));
      return (NULL);
    }
    if( !is.na(disp.date.colname) && !(disp.date.colname %in% names(disp.data)) )
    {
      if( !suppress.warnings ) warning(paste0("Column disp.date.colname='",disp.date.colname,"' must appear in the dispensing data!\n"));
      return (NULL);
    }
    if(anyNA(disp.data[[disp.date.colname]])){
      if( !suppress.warnings ) warning(paste0("Column disp.date.colname='",disp.date.colname,"' cannot contain missing values!\n"));
      return (NULL);
    }
    if( any(!is.na(medication.class.colnames) & !(medication.class.colnames %in% names(disp.data)) & !(medication.class.colnames %in% names(presc.data))) ) # deal with the possibility of multiple column names
    {
      if( !suppress.warnings ) warning(paste0("Column(s) medication.class.colnames=",paste0("'",medication.class.colnames,"'",collapse=",")," must appear in the dispensing and prescribing data!\n"));
      return (NULL);
    }
    if( !is.na(total.dose.colname) && !(total.dose.colname %in% names(disp.data)) )
    {
      if( !suppress.warnings ) warning(paste0("Column total.dose.colname='",total.dose.colname,"' must appear in the dispensing data!\n"));
      return (NULL);
    }
    if(anyNA(disp.data[[total.dose.colname]])){
      if( !suppress.warnings ) warning(paste0("Column total.dose.colname='",total.dose.colname,"' cannot contain missing values!\n"));
      return (NULL);
    }
    if( !is.na(presc.daily.dose.colname) && !(presc.daily.dose.colname %in% names(presc.data)) )
    {
      if( !suppress.warnings ) warning(paste0("Column presc.daily.dose.colname='",presc.daily.dose.colname,"' must appear in the prescribing data!\n"));
      return (NULL);
    }
    if(anyNA(presc.data[[presc.daily.dose.colname]])){
      if( !suppress.warnings ) warning(paste0("Column presc.daily.dose.colname='",presc.daily.dose.colname,"' cannot contain missing values!\n"));
      return (NULL);
    }
    if( !is.na(presc.duration.colname) && !(presc.duration.colname %in% names(presc.data)) )
    {
      if( !suppress.warnings ) warning(paste0("Column presc.duration.colname='",presc.duration.colname,"' must appear in the prescribing data!\n"));
      return (NULL);
    }

    if( visit.colname %in% colnames(presc.data) ) {
      if(anyNA(presc.data[[visit.colname]])){
        if( !suppress.warnings ) warning(paste0("Column visit.colname='",visit.colname,"' cannot contain missing values!\n"));
        return (NULL);
      }
    }

    if( !is.logical(force.presc.renew) && !force.presc.renew %in% names(disp.data) )
    {
      if( !suppress.warnings ) warning(paste0("Column force.presc.renew='",force.presc.renew,"' must appear in the dispensing data!\n"));
      return (NULL);
    }

    if( !is.logical(split.on.dosage.change) && !split.on.dosage.change %in% names(disp.data) )
    {
      if( !suppress.warnings ) warning(paste0("Column split.on.dosage.change='",split.on.dosage.change,"' must appear in the dispensing data!\n"));
      return (NULL);
    }

    if( !trt.interruption %in% c("continue", "discard", "carryover") && !trt.interruption %in% names(disp.data))
    {
      if( !suppress.warnings ) warning(paste0("trt.interruption must be either of 'continue', 'discard',
                                              'carryover', or a column name in the dispensing data!\n"));
      return (NULL);
    }
    if(trt.interruption %in% names(disp.data) && any(!unique(disp.data[[trt.interruption]]) %in% c("continue", "discard", "carryover")))
    {
      unexpected.values <- unique(disp.data[[trt.interruption]][disp.data[[trt.interruption]] %in% c("continue", "discard", "carryover")])

      if( !suppress.warnings ) warning(paste0("Column trt.interruption='",trt.interruption, "' contains unexpected values: ",
                                              unexpected.values,"\n"));
      return (NULL);
    }

    if(".episode" %in% colnames(presc.data)){
      {
        if( !suppress.warnings ) warning("The column name \'.episode\' is used internally, please use another column name.");
        return (NULL);
      }
    }

    if( is.na(date.format) || is.null(date.format) || length(date.format) != 1 || !is.character(date.format) )
    {
      if( !suppress.warnings ) warning(paste0("The date format must be a single string!\n"));
      return (NULL);
    }

  }

  ## Force data to data.table
  if( !inherits(disp.data,"data.table") ) disp.data <- as.data.table(disp.data);
  if( !inherits(presc.data,"data.table") ) presc.data <- as.data.table(presc.data);

  # copy datasets
  disp.data.copy <- data.table(disp.data)
  presc.data.copy <- data.table(presc.data)

  # convert column names
  setnames(presc.data.copy,
           old = c(ID.colname,
                   presc.date.colname,
                   presc.daily.dose.colname,
                   presc.duration.colname),
           new = c("ID",
                   "PRESC.DATE",
                   "DAILY.DOSE",
                   "episode.duration"))

  setnames(disp.data.copy,
           old = c(ID.colname,
                   disp.date.colname,
                   total.dose.colname),
           new = c("ID",
                   "DISP.DATE",
                   "TOTAL.DOSE"))

  # convert dates
  disp.data.copy[,DISP.DATE := as.Date(DISP.DATE, format = date.format)];
  presc.data.copy[,PRESC.DATE := as.Date(PRESC.DATE, format = date.format)];
  if(!is.null(special.periods.data))
  {  ## Force data to data.table
    if( !inherits(special.periods.data,"data.table") ) special.periods.data <- as.data.table(special.periods.data);

    # copy datasets
    special.periods.data.copy <- data.table(special.periods.data)

    setnames(special.periods.data.copy,
             old = c(ID.colname),
             new = c("ID"))

    special.periods.data.copy[,`:=` (DATE.IN = as.Date(DATE.IN, format = date.format),
                                DATE.OUT = as.Date(DATE.OUT, format = date.format))];

    special.periods.data.copy[,SPECIAL.DURATION := as.numeric(DATE.OUT-DATE.IN)];
  } else {special.periods.data.copy <- NULL}

  # force medication class to character
  for(class.colname in medication.class.colnames)
  {
    if(inherits(disp.data.copy[[class.colname]], "factor"))
    {
      disp.data.copy[,(class.colname) := as.character(get(class.colname))];
    }

    if(inherits(presc.data.copy[[class.colname]], "factor"))
    {
      presc.data.copy[,(class.colname) := as.character(get(class.colname))];
    }
  }

  # add prescription duration column if NA is provided
  if( is.na(presc.duration.colname) )
  {
    presc.data.copy[,.PRESC.DURATION := NA]
    presc.duration.colname <- ".PRESC.DURATION"
  }

  # add event ID
  disp.data.copy[,EVENT.ID := 1]

  # helper function to process each patient
  process_patient <- function(pat)
  {
    # helper function to process each medication
    process_medication <- function(med)
    {

      # helper function to process each dispensing event
      process_dispensing_events <- function(event)
      {
        # helper function to compute special intervals
        compute.special.intervals <- function(data,
                                              DATE.IN.colname = "DATE.IN",
                                              DATE.OUT.colname = "DATE.OUT",
                                              TYPE.colname = "TYPE",
                                              CUSTOM.colname = special.periods.method)
          {

          if(CUSTOM.colname %in% colnames(data)){
            setnames(data, old = CUSTOM.colname, new = "CUSTOM")
          } else { data[,CUSTOM := special.periods.method]}


          # convert dates
          data[, (DATE.IN.colname) := as.Date(get(DATE.IN.colname), format = date.format)]
          data[, (DATE.OUT.colname) := as.Date(get(DATE.OUT.colname), format = date.format)]

          # add durations
          data[,DURATION := as.numeric(get(DATE.OUT.colname) - get(DATE.IN.colname))]

          # add episodes
          data[,.episode := seq_len(.N)]

          # melt special episodes
          data.melt <- melt(data,
                            measure.vars = c(DATE.IN.colname, DATE.OUT.colname),
                            variable.name = "EVENT",
                            value.name = "DATE")

          # sort by DATE.IN
          setkeyv(data.melt, cols = c("DATE", ".episode"))

          # add dispensing event
          data.melt <- rbind(data.melt,
                             data.table(ID = pat,
                                        DATE = disp.start.date.i,
                                        EVENT = "DISP.DATE",
                                        .episode = 0),
                             fill = TRUE)

          # find row with end of episode
          data.melt <- rbind(data.melt,
                             data.table(ID = pat,
                                        DATE = end.episode,
                                        EVENT = "episode.end",
                                        .episode = -1),
                             fill = TRUE)

          data.melt[, EVENT := factor(EVENT, levels = c("DATE.OUT", "DISP.DATE", "DATE.IN", "episode.end"))]

          setorderv(data.melt, cols = c("DATE", "EVENT", ".episode"), na.last = TRUE)

          # calculate durations of intersections
          data.melt[,`:=` (DISP.EVENT = 0,
                           CARRYOVER.DURATION = 0,
                           INT.DURATION = as.numeric(shift(DATE, n = 1, type = "lead")-DATE))]

          # find active period
          data.melt[,active.episode := sapply(seq(nrow(data.melt)), function(x) {

            dt <- data.melt[seq(x)]

            closed.episodes <- dt[duplicated(dt[,.episode]),.episode]

            active.episode <- dt[!.episode %in% closed.episodes, suppressWarnings(max(.episode))]

          })]

          # indicate intersections that should be counted
          data.melt[active.episode %in% unique(data.melt[CUSTOM == "continue", .episode]),
                         `:=` (SPECIAL.PERIOD = 1,
                               DISP.EVENT = 1)]
          data.melt[active.episode == 0, DISP.EVENT := 1]

          # calculat durations during carryover
          if( "carryover" %in% unique(data.melt$CUSTOM) ){


            data.melt[active.episode %in% unique(data.melt[CUSTOM == "carryover", .episode]),
                           CARRYOVER.DURATION := INT.DURATION]

            # remove duration during carryover
            data.melt[CARRYOVER.DURATION != 0,
                           INT.DURATION := 0]

         }

          # remove events before dispensing date and after end date
          first.row <- data.melt[EVENT == "DISP.DATE", which = TRUE]

          last.row <- data.melt[EVENT == "episode.end", which = TRUE]

          data.melt <- data.melt[first.row:last.row]

          # identify rows after discard
          data.melt[, .drop := 0]
          if("discard" %in% data$CUSTOM){
            data.melt[,DISP.EVENT := 0]
            data.melt[CUSTOM == "discard", DISP.EVENT := 1]
            data.melt[,.drop := cumsum(DISP.EVENT)]
            data.melt[CUSTOM == "discard" & EVENT == DATE.IN.colname, .drop := .drop-1]

            # remove durations after discard
            data.melt[CUSTOM == "discard", `:=` (DISP.EVENT = 0,
                                                      INT.DURATION = 0)]
          }

          # drop rows after discard
          data.melt.drop <- data.melt[.drop == 0]

          # create intervals of continuous use
          data.melt.drop[,.interval := rleidv(data.melt.drop, cols = "DISP.EVENT")]

          data.melt.drop[DISP.EVENT == 1,.interval := as.integer(.interval+1)]

          # calculate sum of all durations
          sum.duration <- sum(data.melt.drop$INT.DURATION, na.rm = TRUE);

          # if the supply duration is shorter than the sum of the duration
          if(duration.i <= sum.duration)
          {
            # calculate cumulative sum of durations
            data.melt.drop[,cum.duration := cumsum(INT.DURATION)];
            # subset to all rows until supply is exhaused and add 1
            .rows <- data.melt.drop[cum.duration <= duration.i,which=TRUE];
            if( length(.rows) == 0 ) {.rows <- 0};
            data.melt.drop <- data.melt.drop[c(.rows, tail(.rows,1)+1)];

            # calculate remaining duration for last row
            sum.duration <- sum(head(data.melt.drop,-1)$INT.DURATION);
            data.melt.drop[nrow(data.melt.drop), INT.DURATION := duration.i-sum.duration];
          }

          # calculate total duration
          data.melt.drop[,DURATION := sum(INT.DURATION, na.rm = TRUE), by = .interval]

          # calculate duration covered during special intervals
          data.melt.drop[,SPECIAL.DURATION := 0]

          data.melt.drop[SPECIAL.PERIOD == 1, SPECIAL.DURATION := sum(INT.DURATION, na.rm = TRUE), by = .interval]
          # data.melt.drop[(CUSTOM == "continue" & EVENT == DATE.IN.colname) | (shift(CUSTOM, n = 1, type = "lead") == "continue" & shift(EVENT, n = 1, type = "lead") == DATE.OUT.colname), # | (shift(CUSTOM, n = 1, type = "lag") == "continue" & shift(EVENT, n = 1, type = "lag") == DATE.IN.colname),
          #                SPECIAL.DURATION := sum(INT.DURATION, na.rm = TRUE), by = .interval]
          data.melt.drop[,SPECIAL.DURATION := max(SPECIAL.DURATION, na.rm = TRUE), by = .interval]

          # calculate duration NOT covered during special intervals
          data.melt.drop[,CARRYOVER.DURATION := sum(CARRYOVER.DURATION, na.rm = TRUE), by = .interval]

          # subset to first and last row of interval
          events <- data.melt.drop[ data.melt.drop[, .I[c(1L,.N)], by=.interval]$V1 ]
          events[,CUSTOM := last(CUSTOM), by = .interval]

          # convert to wide format with start and end date of intervals
          events[,EVENT := rep(c("DISP.START", "DISP.END"), nrow(events)/2)]

          all.events <- dcast(events, ID + CUSTOM + DURATION + SPECIAL.DURATION + CARRYOVER.DURATION + .interval ~ EVENT, value.var = "DATE")
          setorderv(all.events, cols = ".interval")

          # create event IDs
          all.events[,EVENT.ID := seq(from = event_id, length.out = nrow(all.events), by = 1)]

          # create all events table
          all.events <- cbind(all.events[,c("ID",
                                            "CUSTOM",
                                            "EVENT.ID",
                                            "DISP.START",
                                            "DURATION",
                                            "SPECIAL.DURATION",
                                            "CARRYOVER.DURATION"), with = FALSE],
                              data.table(DAILY.DOSE = as.numeric(presc.dose.i),
                                         episode.start = start.episode,
                                         episode.end = end.episode))

          return(all.events)
        }

        if(exists("debug.mode") && debug.mode==TRUE) print(paste("Event:", event));

         ## !Important: We assume that the prescribed dose can be accomodated with the dispensed medication

        #subset data to event
        curr_disp <- med_disp[event];
        orig.disp.date <- curr_disp[["DISP.DATE"]]

        # if current dispensing event is before first prescription date, don't calculate a duration
        if(orig.disp.date < first_presc[["PRESC.DATE"]])
        {
          med_event <- cbind(curr_disp[,c("ID", medication.class.colnames, "TOTAL.DOSE", "DISP.DATE"), with = FALSE],
                             DISP.START = orig.disp.date,
                             DURATION = 0,
                             DAILY.DOSE = NA,
                             SPECIAL.DURATION = NA);
          # if current dispensing event is after end of last prescription episode, don't calculate a duration (only when last prescription indicates termination)
        } else
        {
          #select prescription episodes ending after the original dispensing date
          episodes <- med_presc[orig.disp.date < episode.end | is.na(episode.end), which = TRUE];

          ## for each prescription episode, calculate the duration with the current dose
          total.dose.i <- curr_disp[["TOTAL.DOSE"]]; #dispensed dose
          presc.dose.i <- 0; # initialize prescibed dose as 0
          disp.start.date.i <- orig.disp.date; #start date of dispensing event

          ## check for carry-over status and adjust start date in case of carry-over from last event
          if( carryover == TRUE){
            if(length(last.disp.end.date) > 0 && !is.na(last.disp.end.date) && last.disp.end.date > disp.start.date.i ) {

              disp.start.date.i <- last.disp.end.date

              #select prescription episodes ending after the original dispensing date
              episodes <- med_presc[disp.start.date.i < episode.end | is.na(episode.end), which = TRUE];
            }
          }

          # if the current dispensing event is after the last prescription episode, don't calculate a duration
          if(length(episodes) == 0 | out.of.presc == TRUE)
          {
            med_event <- cbind(curr_disp[,c("ID", medication.class.colnames, "TOTAL.DOSE", "DISP.DATE"), with = FALSE],
                               DISP.START = orig.disp.date,
                               DURATION = NA,
                               DAILY.DOSE = NA,
                               SPECIAL.DURATION = NA);
          } else
          {
            #select prescription episodes ending after the original dispensing date and add the one immediately before
            curr_med_presc <- data.table(med_presc)

            # if supply should be finished with original dose, collapse consecutive episodes with dosage > 0
            if(split.on.dosage.change == FALSE){

              curr_med_presc[(orig.disp.date < episode.end | is.na(episode.end)) & DAILY.DOSE > 0,POS.DOSE := 1]
              curr_med_presc[,.episode := rleidv(.SD, cols = "POS.DOSE")]

              curr_med_presc[POS.DOSE == 1,episode.start := head(episode.start,1), by = .episode]; # first start date per episode
              curr_med_presc[POS.DOSE == 1,episode.end:= tail(episode.end,1), by = .episode]; # last end date per episode
              curr_med_presc[POS.DOSE == 1,DAILY.DOSE:= head(DAILY.DOSE,1), by = .episode]; # first dosage per episode

              curr_med_presc <- unique(curr_med_presc, by = c("episode.start", "episode.end"), fromLast = TRUE);
              curr_med_presc[,.episode := rleidv(curr_med_presc, cols = c("episode.start", "episode.end"))];

              #select prescription episodes ending after the original dispensing date
              episodes <- curr_med_presc[orig.disp.date < episode.end | is.na(episode.end), which = TRUE];

              }

            # rm.trt.episode <- FALSE; # will be set to TRUE in case of calculations during treatment interruptions

            stop <- 0;

            med_event <- NULL;
            event_id <- 0

            for(episode in episodes)
            {event_id <- event_id + 1

              presc.dose.i <- curr_med_presc[[episode,"DAILY.DOSE"]]; # prescribed daily dose
              start.episode <- curr_med_presc[episode,episode.start];
              end.episode <- curr_med_presc[episode,episode.end];

              if(presc.dose.i == 0) # if event happens during treatment interruption (prescribed dose = 0), check what to do
              {
                if(trt.interruption == "continue") # if trt.interruption is set to continue, continue with last prescribed dose
                {
                  presc.dose.i <- curr_med_presc[[episode-1,"DAILY.DOSE"]];

                  # adjust prescription episode to previous episode
                  start.episode <- curr_med_presc[episode-1,episode.start];
                  end.episode <- curr_med_presc[episode-1,episode.end];

                  stop <- 1;

                  # rm.trt.episode <- TRUE; # make sure that prescription start- and end date are removed at the end
                } else if(trt.interruption == "discard") # if trt.interruption is set to discard, don't calculate anything
                {
                  if(is.null(med_event))
                  {
                    med_event <- cbind(curr_disp[,c("ID", medication.class.colnames, "TOTAL.DOSE", "DISP.DATE"), with = FALSE],
                                       EVENT.ID = event_id,
                                       DISP.START = disp.start.date.i,
                                       DURATION = 0,
                                       DAILY.DOSE = NA,
                                       SPECIAL.DURATION = NA);
                  }

                  break
                } else
                {
                  episode <- episode + 1; # else skip to next episode
                  next;
                }
              }

              # if disp.start.date.i is after end.episode date, go to next episode.
              if( !is.na(curr_med_presc[episode,episode.end]) & disp.start.date.i >= curr_med_presc[episode,episode.end] ) {
                next;
              }

              # if it is not the first episode, adjust supply start date to prescription start date
              if(episode != episodes[1]) disp.start.date.i <- curr_med_presc[episode,episode.start];

              duration.i <- total.dose.i/presc.dose.i; # calculate duration

              disp.end.date.i <- disp.start.date.i + duration.i; # calculate end date of supply

              # add special durations during the supply period
              special.periods.duration.i <- 0;
              if(nrow(med_special.periods_events) != 0 & !is.na(duration.i))
              {
                # check for special durations within the episode
                med_special.periods_events_i <- med_special.periods_events[(DATE.IN <= end.episode|is.na(end.episode)) & DATE.OUT > start.episode];

                if(nrow(med_special.periods_events_i) > 0)
                {
                  all.events <- compute.special.intervals(med_special.periods_events_i);

                  event_id <- last(all.events$EVENT.ID)

                  sum.duration <- sum(all.events$DURATION, na.rm = TRUE)

                  # if last line is "discard", create med_event
                  if(!is.na(last(all.events$CUSTOM)) && last(all.events$CUSTOM) == "discard") {

                    med_event <-  rbind(med_event,
                                        cbind(curr_disp[,c("ID", medication.class.colnames, "TOTAL.DOSE", "DISP.DATE"), with = FALSE],
                                              all.events[,3:10]),
                                        fill = TRUE);

                    break;
                  } else if( duration.i == sum.duration ) # if supply is equal to the sum of durations
                  {
                    med_event <-  rbind(med_event,
                                        cbind(curr_disp[,c("ID", medication.class.colnames, "TOTAL.DOSE", "DISP.DATE"), with = FALSE],
                                              all.events[,3:10]),
                                        fill = TRUE);

                    break;

                  } else if(is.na(last(all.events$episode.end))) # if last event is not terminated
                  {
                    all.events[nrow(all.events), DURATION := DURATION + (duration.i-sum.duration)];

                    med_event <-  rbind(med_event,
                                        cbind(curr_disp[,c("ID", medication.class.colnames, "TOTAL.DOSE", "DISP.DATE"), with = FALSE],
                                              all.events[,3:10]),
                                        fill = TRUE);
                    break;
                  } else # if supply duration is longer than the sum of the durations
                  {
                    # calculate the carryover dose
                    oversupply <- duration.i-sum.duration; # calculate remaining days of oversupply
                    total.dose.i <- presc.dose.i*oversupply; # calculate remaining total dose

                    med_event <-  rbind(med_event,
                                        cbind(curr_disp[,c("ID", medication.class.colnames, "TOTAL.DOSE", "DISP.DATE"), with = FALSE],
                                              all.events[,3:10]),
                                        fill = TRUE);
                    next;
                  }
                }
              }

              # check various parameters to decide wheter to stop or continue

              # check if end of supply is before end of episode OR last row of prescription episodes is reached
              if( disp.end.date.i < curr_med_presc[episode,episode.end] | episode == last(episodes) )
              {
                stop <- 1;
              } else {
                episode <- episode + 1; # get next prescription episode
                next.presc.dose <- curr_med_presc[[episode,"DAILY.DOSE"]]; # get next episode's dosage

                # if there is a treatment interruption and trt.interruption is set to continue, stop
                if( next.presc.dose == 0 & trt.interruption == "continue" ) stop <- 1;

                # if there is no treatment interruption, but a dosage change and split.on.dosage.change is set FALSE, stop
                if( next.presc.dose != 0 & next.presc.dose != presc.dose.i & split.on.dosage.change == FALSE ) stop <- 1;
              }

              if( stop == 1 )
              {
                # if( rm.trt.episode == TRUE )
                # {
                #   start.episode <- as.Date(NA, format = date.format);
                #   end.episode <- as.Date(NA, format = date.format);
                # }

                med_event <- rbind(med_event,
                                   cbind(curr_disp[,c("ID", medication.class.colnames, "TOTAL.DOSE", "DISP.DATE"), with = FALSE],
                                         data.table(EVENT.ID = event_id,
                                                    DISP.START = disp.start.date.i,
                                                    DURATION = as.numeric(duration.i),
                                                    episode.start = start.episode,
                                                    episode.end = end.episode,
                                                    DAILY.DOSE = as.numeric(presc.dose.i),
                                                    SPECIAL.DURATION = as.numeric(special.periods.duration.i))),
                                   fill = TRUE);
                break;
              } else
              {
                duration.i <- end.episode - disp.start.date.i; # calculate duration until end of episode
                oversupply <- disp.end.date.i - end.episode; # calculate remaining days of oversupply
                total.dose.i <- presc.dose.i*oversupply; # calculate remaining total dose

                # if( rm.trt.episode == TRUE )
                # {
                #   start.episode <- as.Date(NA, format = date.format);
                #   end.episode <- as.Date(NA, format = date.format);
                # }

                #create medication event
                med_event <- rbind(med_event,
                                   cbind(curr_disp[,c("ID", medication.class.colnames, "TOTAL.DOSE", "DISP.DATE"), with = FALSE],
                                         data.table(EVENT.ID = event_id,
                                                    DISP.START = disp.start.date.i,
                                                    DURATION = as.numeric(duration.i),
                                                    episode.start = start.episode,
                                                    episode.end = end.episode,
                                                    DAILY.DOSE = as.numeric(presc.dose.i),
                                                    SPECIAL.DURATION = as.numeric(special.periods.duration.i))),
                                   fill = TRUE);
              }
            }
            med_event;
          }
        }
      }

      if(exists("debug.mode") && debug.mode==TRUE) print(paste("Medication:", med));

      ## subset data to medication

      setkeyv(pat_disp, medication.class.colnames);
      setkeyv(pat_presc, medication.class.colnames);

      med_disp <- pat_disp[list(disp_presc[med, medication.class.colnames, with = FALSE])];

      med_presc <- pat_presc[list(disp_presc[med, medication.class.colnames, with = FALSE])];

      setkeyv(med_disp, cols = "DISP.DATE");
      setkeyv(med_presc, cols = "PRESC.DATE");

      med_special.periods_events <- data.table(special.periods_events)
      if( !is.null(special.periods.data) )
      {
        special.colnames <- intersect(medication.class.colnames, colnames(special.periods.data.copy))

        if( length(special.colnames) > 0 ) {
          setkeyv(special.periods_events, special.colnames);
          med_special.periods_events <- special.periods_events[list(disp_presc[med, special.colnames, with = FALSE])];
        }

        setkeyv(med_special.periods_events, cols = "DATE.IN")
      }

      # determine date of initial prescription
      first_presc <- med_presc[1];

      # determine date of initial dispense
      first_disp <- med_disp[["DISP.DATE"]][1];

      #if force.presc.renew, trt.interruption, and split.on.dosage.change are not set globally, set for medication based on first dispensing event
      if( !is.logical(force.presc.renew) )
      {
        force.presc.renew <- as.logical(first_disp[[force.presc.renew]]);
      }
      if( !trt.interruption %in% c("continue", "discard", "carryover") )
      {
        trt.interruption <- as.logical(first_disp[[trt.interruption]]);
      }

      if( !is.logical(split.on.dosage.change) )
      {
        split.on.dosage.change <- as.logical(first_disp[[split.on.dosage.change]]);
      }

      ## calculate treatment interruptions and end of prescription date

      ## determine end of prescription and prescription interruptions if prescription reneval is enforced for each subsequent prescription event (requires the "visit" column)
      presc_interruptions <- data.table(NULL);
      if( force.presc.renew == TRUE )
      {
        presc_visit <- presc_events[[visit.colname]] %in% unique(med_presc[[visit.colname]]); # determine for each visit if medication was prescribed

        first_presc_event <- head(which(presc_visit),1); # extract first prescription event
        last_presc_event <- tail(which(presc_visit),1); # extract last prescription event

        presc_omit <- which(!presc_visit)[which(!presc_visit) > first_presc_event & which(!presc_visit) < last_presc_event]; # identify visits between first and last prescription with missing prescriptions

        interruption_dates <- presc_events[["PRESC.DATE"]][presc_omit]; # determine dates of treatment interruptions

        presc_interruptions <- med_presc[rep(1, length(presc_omit))]; # create table with one row for each interruption

        presc_interruptions[, c(visit.colname, "PRESC.DATE", "DAILY.DOSE", "episode.duration") :=
                              list(presc_events[[visit.colname]][presc_omit], interruption_dates, 0, NA)]; # adjust variables

        med_presc <- rbind(med_presc, presc_interruptions); # bind to existing prescriptions
        setkeyv(med_presc, cols = "PRESC.DATE"); # order by date

        med_presc[,.episode := rleidv(med_presc, cols = "DAILY.DOSE")]; # add counter for treatment episodes
      }

      setorder(med_presc);

      ## construct treatment episodes
      # create new .episode counter
      med_presc[,.episode := rleidv(med_presc, cols = c("DAILY.DOSE", "episode.duration"))];

      # if consecutive episodes with set end date, increase .episode counter

      if( nrow(med_presc) > 2 )
      {
        for( n in 2:(nrow(med_presc)))
        {
          if( !is.na(med_presc[n,"episode.duration", with = FALSE]) & !is.na(med_presc[n-1,"episode.duration", with = FALSE]) )
          {
            med_presc[n:nrow(med_presc), .episode := as.integer(.episode + 1)];
          }
        }
      } else if( nrow(med_presc) == 2 )
      {
        med_presc[!is.na(shift(episode.duration, type = "lag")) & !is.na(episode.duration), .episode := as.integer(.episode + 1)];
      }

      # add episodes with same dose but set end date to last episode
      .row <- med_presc[is.na(shift(episode.duration, type = "lag")) & shift(DAILY.DOSE, type = "lag") == DAILY.DOSE & !is.na(episode.duration), which = TRUE];
      if( length(.row)>0 )
      {
        med_presc[.row:nrow(med_presc),.episode := as.integer(.episode-1)];
      }

      ## set start and end of prescription dates per group
      med_presc[, `:=` (episode.start = PRESC.DATE, # set prescription date as start date
                        episode.end = PRESC.DATE)]; # set end date to prescription date ...

      med_presc[,episode.end := shift(episode.end, type = "lead")]; # ... and shift end dates up by one

      # adjust end date if prescription duration is provided and change start date of following prescriptions
      med_presc[!is.na(episode.duration) & ((PRESC.DATE + episode.duration) <= episode.end | is.na(episode.end)), episode.end := PRESC.DATE + episode.duration]; # only if prescription ends before the current end prescription date!
      end.limited.presc <- head(med_presc,-1)[!is.na(episode.duration) & ((PRESC.DATE + episode.duration) <= episode.end | is.na(episode.end))]$episode.end; #don't include last prescription episode
      med_presc[shift(!is.na(episode.duration), type = "lag") & shift((PRESC.DATE + episode.duration) <= episode.end, type = "lag"), episode.start := end.limited.presc];
      med_presc[PRESC.DATE>episode.start & DAILY.DOSE != 0,episode.start:=PRESC.DATE];

      # combine episodes with set durations with previous episodes of same dosage but unrestricted duration
      med_presc[shift(DAILY.DOSE,type="lag")==DAILY.DOSE & !is.na(shift(episode.duration,type="lag")) & shift(episode.end, type = "lag") == episode.start, .episode := as.integer(.episode-1)];

      # fill in start and end dates by group
      med_presc[,episode.start := head(episode.start,1), by = .episode]; # first start date per episode
      med_presc[,episode.end:= tail(episode.end,1), by = .episode]; # last end date per episode
      med_presc[,PRESC.DATE := min(PRESC.DATE), by = .episode]; # set PRESC.DATE to first start date

      # collapse episodes
      med_presc <- unique(med_presc, by = ".episode", fromLast = TRUE);
      med_presc[,.episode := rleidv(med_presc, cols = c("episode.start", "episode.end"))];

      # remove episodes where end date is before start date
      rm.episode <- med_presc[episode.end <= episode.start, which = TRUE];
      if( length(rm.episode) > 0 )
      {
        med_presc <- med_presc[-rm.episode];
      }
      med_presc[,.episode := rleidv(med_presc)];

      # collapse consecutive episodes where end date of the former is before start date of the latter
      med_presc[shift(episode.end,type = "lag") > episode.start & shift(DAILY.DOSE,type = "lag") == DAILY.DOSE,
                .episode := as.integer(.episode-1)];
      med_presc[,episode.start := head(episode.start,1), by = .episode]; # first start date per episode
      med_presc[,episode.end:= tail(episode.end,1), by = .episode]; # last end date per episode
      med_presc <- unique(med_presc, by = ".episode");
      med_presc[,.episode := rleidv(med_presc, cols = c("episode.start", "episode.end"))];

      # add treatment interruptions

      med_presc <- rbind(med_presc,med_presc[shift(episode.start,type = "lead")!=episode.end][,c("DAILY.DOSE", "episode.start", ".episode") := list(0, episode.end, 0)]);
      setorder(med_presc, episode.start, episode.end);
      end.trt.interruptions <- med_presc[shift(episode.end,type = "lag")!=episode.start]$episode.start;
      med_presc[.episode == 0, episode.end := end.trt.interruptions];

      if( force.init.presc == TRUE )
      {
        # if initial dispense is before initial prescription, adjust date of initial prescription to match initial dispense
        # but only if first prescription is unlimited
        if( first_disp < first(med_presc[["episode.start"]]) & is.na(head(med_presc[["episode.duration"]],1)) )
        {
          # adjust first prescription date
          first_presc[1, PRESC.DATE := first_disp];
          med_presc[1, episode.start := first_disp];
        }
      }

      ## calculate medication events for "simple" events not extending over multiple episodes or affected by special periods
      # add prescription events to dispensing events

      for( i in 1:nrow(med_presc) )
      {
        med_disp[DISP.DATE >= med_presc[i,episode.start] & (DISP.DATE < med_presc[i,episode.end] | is.na(med_presc[i,episode.end])),
                 c("episode.start", "episode.end", "DAILY.DOSE") := list(med_presc[i,episode.start], med_presc[i,episode.end],med_presc[i,DAILY.DOSE])];
      }

      med_disp[,DURATION := (TOTAL.DOSE)/(DAILY.DOSE)];
      med_disp[,`:=` (DISP.START = DISP.DATE,
                      DISP.END = DISP.DATE+DURATION)];

      med_disp[DISP.END > episode.end, .out := 1];

      # add special periods to dispensing events
      med_disp[,.special.periods := as.numeric(NA)];

      if( nrow(med_special.periods_events) != 0 ){

         for( i in 1:nrow(med_special.periods_events) )
        {
          med_disp[(DISP.END >= med_special.periods_events[i,DATE.IN] & DISP.START < med_special.periods_events[i,DATE.OUT])|(DISP.START >= med_special.periods_events[i,DATE.IN] & DISP.START < med_special.periods_events[i,DATE.OUT]),
                   .special.periods := 1];
        }
      }

      med_disp[DURATION == Inf | .out == 1 | .special.periods == 1, process.seq := 1]
      med_disp[,process.seq.num := rleidv(process.seq)]

      medication_events_rest <- NULL;

      out.of.presc <- FALSE # set flag for carryover processing

      if(carryover == TRUE){

        # compute carryover
        med_disp[,carryover.from.last := as.numeric(shift(DISP.START+DURATION, type = "lag")-DISP.START)]
        med_disp[1,carryover.from.last := 0]
        med_disp[,carryover.total := cumsum(carryover.from.last)]

        # get first row with carryover
        index <- suppressWarnings(min(which(med_disp$carryover.total > 0)))

        if(index <= nrow(med_disp)){
          med_disp[index:nrow(med_disp), process.seq := 1]
        }

        # create medication events before first carryover event
        medication_events <- med_disp[is.na(process.seq) & process.seq.num == 1,
                                      c("ID",
                                        medication.class.colnames,
                                        "TOTAL.DOSE",
                                        "DISP.DATE",
                                        "EVENT.ID",
                                        "DISP.START",
                                        "DURATION",
                                        "DAILY.DOSE",
                                        "episode.start",
                                        "episode.end"), with = FALSE];
        medication_events[,SPECIAL.DURATION := 0];

        # subset to events with carryover or special periods
        med_disp <- med_disp[process.seq == 1 | process.seq.num > 1];

        ## apply process_dispensing_events to each dispensing event
        last.disp.end.date <- last(medication_events[,DISP.START + DURATION])
        #carryover.total <- 0#ifelse(nrow(medication_events) > 0, last(medication_events$carryover.total), 0)

        if( nrow(med_disp) > 0 )
        {for(i in 1:nrow(med_disp)){

            medication_events_i <- process_dispensing_events(event = i)

            medication_events_rest <- rbind(medication_events_rest, medication_events_i, fill = TRUE)

            # if DURATION is NA, set flag for all future events
            if(is.na(last(medication_events_i[,DURATION]))) {
              out.of.presc <- TRUE
            } else {

              # cache last dispensing end date
              last.disp.end.date <- last(medication_events_i[, DISP.START + DURATION])
            }
          }
        }

      } else {
        medication_events <- med_disp[is.na(process.seq),
                                      c("ID",
                                        medication.class.colnames,
                                        "TOTAL.DOSE",
                                        "DISP.DATE",
                                        "EVENT.ID",
                                        "DISP.START",
                                        "DURATION",
                                        "DAILY.DOSE",
                                        "episode.start",
                                        "episode.end"), with = FALSE];
        medication_events[,SPECIAL.DURATION := 0];

        med_disp <- med_disp[process.seq == 1];

        if( nrow(med_disp) > 0 ) {

          medication_events_rest <- do.call(rbindlist,
                                            list(l = lapply(1:nrow(med_disp), FUN = function(i) process_dispensing_events(event = i)),
                                                 fill = TRUE));
        }
      }

      medication_events <- rbind(medication_events, medication_events_rest, fill = TRUE);

      setorderv(medication_events,cols=c("DISP.DATE", "DISP.START"));

      if( force.presc.renew == TRUE ) # add number of prescription interruptions
      {
        tot.presc.interruptions <- nrow(med_presc[DAILY.DOSE==0]);

        medication_events[,tot.presc.interruptions := tot.presc.interruptions];
      }

      if( split.on.dosage.change == TRUE ) # add number of dosage changes
      {
        tot.dosage.changes <- (nrow(med_presc) - 1 - 2*nrow(med_presc[DAILY.DOSE==0]));

        medication_events[,tot.dosage.changes := tot.dosage.changes];
      }

      # presc_episode_no_dispense <- med_presc[!medication_events[,c("DAILY.DOSE","episode.start","episode.end")],
      #                                         on = c("DAILY.DOSE","episode.start", "episode.end")];
      #
      # presc_episode_no_dispense[,c(".episode","VISIT", "episode.duration", "PRESC.DATE") := NULL];
      #
      # medication_events <- rbind(medication_events, presc_episode_no_dispense, fill = TRUE);

      # add episode number
      med_presc <- med_presc[DAILY.DOSE != 0, episode.ID := seq(.N)];

      # calculate duration
      med_presc[,episode.duration := as.numeric(episode.end-episode.start)];

      # compute prescription events
      prescription_events <- med_presc[DAILY.DOSE != 0,
                                       c("ID",
                                         medication.class.colnames,
                                         "DAILY.DOSE",
                                         "episode.ID",
                                         "episode.start",
                                         "episode.duration",
                                         "episode.end"), with = FALSE]

      return(list(DURATIONS = medication_events,
                  PRESCRIPTION_EPISODES = prescription_events));

################### end of process_medication ###################
    }

    if(exists("debug.mode") && debug.mode==TRUE) print(paste("Patient:",pat));

    # subset data to patient
    pat_disp <- disp.data.copy[ID == pat, c("ID",
                                       "DISP.DATE",
                                       "EVENT.ID",
                                       medication.class.colnames,
                                       "TOTAL.DOSE"), with = FALSE];

    pat_presc <- presc.data.copy[ID == pat, c("ID",
                                         "PRESC.DATE",
                                         medication.class.colnames,
                                         "DAILY.DOSE",
                                         "episode.duration"), with = FALSE];
    if(visit.colname %in% colnames(presc.data.copy)){
      pat_presc <- cbind(presc.data.copy[ID == pat, visit.colname, with = FALSE], pat_presc);
    };
    # sort by DCI
    setkeyv(pat_disp, cols = medication.class.colnames);
    setkeyv(pat_presc, cols = medication.class.colnames);

    # extract unique dispensed/prescribed DCIs
    disp_unique <- unique(pat_disp[,c(medication.class.colnames), with = FALSE]);
    presc_unique <- unique(pat_presc[,c(medication.class.colnames), with = FALSE]);

    # extract medications present in both dispensing and prescription database (by DCI, Unit, and Form)
    disp_presc <- merge(disp_unique, presc_unique, by = c(medication.class.colnames), all=FALSE);

    # extract unique dispensed/prescribed DCIs not present in both databases
    disp_no_presc <- disp_unique[!presc_unique];
    presc_no_disp <- presc_unique[!disp_unique];

    #create visits if not supplied
    if( !visit.colname %in% colnames(presc.data.copy) )
    {
      presc_events <- unique(pat_presc[,"PRESC.DATE"]);
      presc_events[,(visit.colname) := 0:(nrow(presc_events)-1)];
      pat_presc <- merge(pat_presc, presc_events, by = "PRESC.DATE");
      setorderv(pat_presc, medication.class.colnames);
    } else
    {
      presc_events <- unique(pat_presc[,c("PRESC.DATE", visit.colname), with = FALSE]); # extract prescription instances
    }

    # if duplicate visit numbers for different dates or vice versa, throw an error
    if( length(unique(presc_events[["PRESC.DATE"]])) != nrow(presc_events) )
    {
      {
        if( !suppress.warnings ) warning("Prescription dates and visit number don't match for patient Nr.", pat);
        return (NULL);
      }
    }

    # extract special periods
    if( !is.null(special.periods.data) )
    {
      special.periods_events <- special.periods.data.copy[ID == pat];
    } else
    {
      special.periods_events <- data.table(NULL);
    }

    setkeyv(presc_events, cols = "PRESC.DATE");

    # apply process_medication() function to each medication present in both databses
    patient_events <- NULL;
    if( nrow(disp_presc) != 0 )
    {
      patient_events <- lapply(1:nrow(disp_presc), FUN = function(i) process_medication(med = i));

      # patient_events <- do.call(rbindlist, list(l = lapply(1:nrow(disp_presc), FUN = function(i) process_medication(med = i)),
      #                                           fill = TRUE));
    }

    setkeyv(pat_disp, cols = medication.class.colnames);
    setkeyv(pat_presc, cols = medication.class.colnames);
    #
    patient_events[[1]][[1]] <- rbind(pat_disp[list(disp_no_presc[,medication.class.colnames, with = FALSE]), c("ID", "DISP.DATE", medication.class.colnames, "TOTAL.DOSE"), with = FALSE],
                                      pat_presc[list(presc_no_disp[,medication.class.colnames, with = FALSE]), c("ID", medication.class.colnames, "DAILY.DOSE"), with = FALSE],
                                      patient_events[[1]][[1]],
                                      fill = TRUE);

    # update progress bar
    if(progress.bar == TRUE) { setTxtProgressBar(pb, getTxtProgressBar(pb)+1) };

    patient_events;
  }

  # extract IDs of all patients present in dispensing and prescription database
  disp_presc_IDs <- sort(intersect(disp.data.copy[["ID"]], presc.data.copy[["ID"]]));

  # progress bar
  if(progress.bar == TRUE) {
    pb <- txtProgressBar(min = 0, max = length(disp_presc_IDs), style = 3);
  }

  # apply process_patient function
  setkeyv(disp.data.copy, cols = "ID");
  setkeyv(presc.data.copy, cols = "ID");

  events_output_list <- lapply(disp_presc_IDs, FUN = function(i) process_patient(pat = i));

  events_output_durations <- do.call(rbindlist, list(l = lapply(events_output_list, FUN = function(i) {
                             do.call(rbindlist, list(l = lapply(i, FUN = function(j) {
                                     j[[1]] }), fill = TRUE)) }), fill = TRUE));

  events_output_prescriptions <- do.call(rbindlist, list(l = lapply(events_output_list, FUN = function(i) {
                                 do.call(rbindlist, list(l = lapply(i, FUN = function(j) {
                                         j[[2]] }), fill = TRUE)) }), fill = TRUE));


  # events_output <- do.call(rbindlist, list(l = lapply(disp_presc_IDs, FUN = function(i) process_patient(pat = i)),
  #                                          fill = TRUE));

  # key by ID, medication class, and dispensing date
  setkeyv(events_output_durations, cols = c("ID", medication.class.colnames, "DISP.DATE"));
  setkeyv(events_output_prescriptions, cols = c("ID", medication.class.colnames));

  # convert column names
  setnames(events_output_durations,
           old = c("ID",
                   "DISP.DATE",
                   "DAILY.DOSE",
                   "TOTAL.DOSE"),
           new = c(ID.colname,
                   disp.date.colname,
                   presc.daily.dose.colname,
                   total.dose.colname)
  )

  setnames(events_output_prescriptions,
           old = c("ID",
                   "DAILY.DOSE"),
           new = c(ID.colname,
                   presc.daily.dose.colname)
  )

  # only return special periods for selected patients
  if(!is.null(special.periods.data)) {
    special.periods.data.copy <- special.periods.data.copy[ID %in% disp_presc_IDs]
  } else {special.periods.data.copy <- NULL}

if(!is.null(special.periods.data.copy)) {
  setnames(special.periods.data.copy,
           old = c("ID"),
           new = c(ID.colname))
}


if(progress.bar == TRUE) { close(pb) }

  attributes(events_output_durations)$carryover <- carryover

  if( !return.data.table )
  {
    events_output_durations <- as.data.frame(events_output_durations);
    events_output_prescriptions <- as.data.frame(events_output_prescriptions)

  }

  # set order of column names
  opt_cols <- c("SPECIAL.DURATION","tot.presc.interruptions","tot.dosage.changes","CARRYOVER.DURATION")
  opt_cols <- opt_cols[opt_cols %in% names(events_output_durations)]

  colorder <- c(ID.colname,
                medication.class.colnames,
                disp.date.colname,
                total.dose.colname,
                presc.daily.dose.colname,
                "EVENT.ID",
                "DISP.START",
                "DURATION",
                "episode.start",
                "episode.end",
                opt_cols)

  setcolorder(events_output_durations, colorder)

  summary <- "Event durations based on dispensing, prescription, and other data, which can be used with the CMA constructors in AdhereR."

  structure(list("event_durations" = events_output_durations,
                 "prescription_episodes" = events_output_prescriptions,
                 "special_periods" = special.periods.data,
                 "ID.colname" = ID.colname,
                 "medication.class.colnames" = medication.class.colnames,
                 "disp.date.colname" = disp.date.colname,
                 "total.dose.colname" = total.dose.colname,
                 "presc.date.colname" = presc.date.colname,
                 "presc.daily.dose.colname"  = presc.daily.dose.colname,
                 "presc.duration.colname" = presc.duration.colname,
                 "visit.colname"  = visit.colname,
                 "split.on.dosage.change" = split.on.dosage.change,
                 "force.init.presc" = force.init.presc,
                 "force.presc.renew" = force.presc.renew,
                 "trt.interruption" = trt.interruption,
                 "special.periods.method" = special.periods.method,
                 "date.format" = date.format),
            class = "event_durations");

}