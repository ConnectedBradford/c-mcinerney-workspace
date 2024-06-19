# Auxiliary function: call a given function sequentially or in parallel:
.compute.function <- function(fnc, # the function to compute
                              fnc.ret.vals=1, # how many distinct values (as elements in a list) does fnc return (really useful for binding results from multi-cpu processing)?
                              # Parallel processing:
                              parallel.backend=c("none","multicore","snow","snow(SOCK)","snow(MPI)","snow(NWS)")[1], # parallel backend to use
                              parallel.threads="auto", # specification (or number) of parallel threads
                              # The parameters with which to call the function:
                              # - NULL indicates that they are not used at all; the values, including defaults, must be fed by the caller
                              # - all consistency checks have been already been done in the caller (except for the parallel procssing params: these are checked here)
                              data=NULL, # this is a per-event *data.table* already keyed by patient ID and event date!
                              ID.colname=NULL, # the name of the column containing the unique patient ID
                              event.date.colname=NULL, # the start date of the event in the date.format format
                              event.duration.colname=NULL, # the event duration in days
                              event.daily.dose.colname=NULL, # the prescribed daily dose
                              medication.class.colname=NULL, # the classes/types/groups of medication
                              event.interval.colname=NULL, # contains number of days between the start of current event and the start of the next
                              gap.days.colname=NULL, # contains the number of days when medication was not available
                              carryover.within.obs.window=NULL, # if TRUE consider the carry-over within the observation window
                              carryover.into.obs.window=NULL, # if TRUE consider the carry-over from before the starting date of the observation window
                              carry.only.for.same.medication=NULL, # if TRUE the carry-over applies only across medication of same type
                              consider.dosage.change=NULL, # if TRUE carry-over is adjusted to reflect changes in dosage
                              followup.window.start=NULL, # if a number, date earliest event per participant + number of units, otherwise a date.format date or variable date
                              followup.window.start.unit=NULL, # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!)
                              followup.window.duration=NULL, # the duration of the follow-up window in the time units given below
                              followup.window.duration.unit=NULL, # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!)
                              observation.window.start=NULL, # the number of time units relative to followup.window.start, otherwise a date.format date or variable date
                              observation.window.start.unit=NULL, # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!)
                              observation.window.duration=NULL, # the duration of the observation window in time units
                              observation.window.duration.unit=NULL, # the time units; can be "days", "weeks", "months" or "years" (if months or years, using an actual calendar!)
                              date.format=NULL, # the format of the dates used in this function
                              suppress.warnings=NULL,
                              suppress.special.argument.checks=FALSE
)
{
  # Quick decision for sequential processing:
  if( parallel.backend == "none" || (is.numeric(parallel.threads) && parallel.threads == 1) )
  {
    # Single threaded: simply call the function with the given data:
    return (fnc(data=data,
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
                suppress.special.argument.checks=suppress.special.argument.checks
           ));
  }


  # Could be Multicore:
  if( parallel.backend == "multicore" )
  {
    if( .Platform$OS.type == "windows" )
    {
      # Can't do multicore on Windows!
      if( !suppress.warnings ) .report.ewms(paste0("Parallel processing backend \"multicore\" is not currently supported on Windows: will use SNOW instead.\n"), "warning", ".compute.function", "AdhereR");
      parallel.backend <- "snow"; # try to do SNOW...
    } else
    {
      # Okay, seems we're on some sort of *NIX, so can do multicore!

      # Pre-process parallel.threads:
      if( parallel.threads == "auto" )
      {
        parallel.threads <- getOption("mc.cores", 2L);
      } else if( is.na(parallel.threads) || !is.numeric(parallel.threads) || (parallel.threads < 1) || (parallel.threads %% 1 != 0) )
      {
        if( !suppress.warnings ) .report.ewms(paste0("Number of parallel processing threads \"",parallel.threads,"\" must be either \"auto\" or a positive integer; forcing \"auto\".\n"), "warning", ".compute.function", "AdhereR");
        parallel.threads <- getOption("mc.cores", 2L);
      }

      # Pre-split the participants into a number of chunks equal to the number of threads to reduce paying the overheads multiple times
      # and call the function for each thread in parallel:
      patids <- unique(data[,get(ID.colname)]); # data is a data.table, so careful with column selection!
      if( length(patids) < 1 ) return (NULL);
      tmp <- parallel::mclapply( lapply(parallel::splitIndices(length(patids), min(parallel.threads, length(patids))), function(i) patids[i]), # simulate snow::clusterSplit() to split patients by thread
                                 function(IDs) fnc(data=data[list(IDs),], # call the function sequentially for the patients in the current chunk
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
                                                   suppress.special.argument.checks=suppress.special.argument.checks),
                                 mc.cores=parallel.threads, # as many cores as threads
                                 mc.preschedule=FALSE # don't preschedule as we know that we have relatively few jobs to do
                               );

      # Combine the results (there may be multiple returned data.frames intermingled!)
      if( fnc.ret.vals == 1 )
      {
        # Easy: just one!
        ret.val <- data.table::rbindlist(tmp);
      } else
      {
        # Combine them in turn:
        ret.val <- lapply(1:fnc.ret.vals, function(i)
        {
          x <- data.table::rbindlist(lapply(tmp, function(x) x[[i]]));
        });
        if( length(tmp) > 0 ) names(ret.val) <- names(tmp[[1]]);
      }
      return (ret.val);
    }
  }

  # Could be SNOW:
  cluster.type <- switch(parallel.backend,
                         "snow"=, "snow(SOCK)"="SOCK", # socket cluster
                         "snow(MPI)"="MPI",            # MPI cluster (required Rmpi)
                         "snow(NWS)"="NWS",            # NWS cluster (requires nws)
                         NA                            # unknown type of cluster
                        );
  if( is.na(cluster.type) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("Unknown parallel processing backend \"",parallel.backend,"\": will force sequential (\"none\").\n"), "warning", ".compute.function", "AdhereR");
    # Force single threaded: simply call the function with the given data:
    return (fnc(data=data,
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
                suppress.special.argument.checks=suppress.special.argument.checks
           ));
  }

  patids <- unique(data[,get(ID.colname)]); # data is a data.table, so careful with column selection!
  if( length(patids) < 1 ) return (NULL);

  # Pre-process parallel.threads:
  if( length(parallel.threads) == 1 )
  {
    if( parallel.threads == "auto" )
    {
      parallel.threads <- 2L;
    } else if( is.na(parallel.threads) || (is.numeric(parallel.threads) && ((parallel.threads < 1) || (parallel.threads %% 1 != 0))) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Number of parallel processing threads \"",parallel.threads,"\", if numeric, must be either a positive integer; forcing \"auto\".\n"), "warning", ".compute.function", "AdhereR");
      parallel.threads <- 2L;
    }
    if( is.numeric(parallel.threads) ) parallel.threads <- min(parallel.threads, length(patids)); # make sure we're not starting more threads than patients
  }

  # Attempt to create the SNOW cluster:
  cluster <- parallel::makeCluster(parallel.threads, # process only "auto", otherwise trust makeCluster() to interpret the parameters
                                   type = cluster.type);
  if( is.null(cluster) )
  {
    if( !suppress.warnings ) .report.ewms(paste0("Failed to create the cluster \"",parallel.backend,"\" with parallel.threads \"",parallel.threads,"\": will force sequential (\"none\").\n"), "warning", ".compute.function", "AdhereR");
    # Force single threaded: simply call the function with the given data:
    return (fnc(data=data,
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
                suppress.special.argument.checks=suppress.special.argument.checks
           ));
  }

  # Pre-split the participants into a number of chunks equal to the number of created cluster nodes to reduce paying the overheads multiple times
  # and call the function for each cluster in parallel:
  tmp <- parallel::parLapply(cluster,
                             parallel::clusterSplit(cluster, patids),
                             function(IDs) fnc(data=data[list(IDs),], # call the function sequentially for the patients in the current chunk
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
                                               suppress.special.argument.checks=suppress.special.argument.checks));

  parallel::stopCluster(cluster); # stop the cluster

  # Combine the results (there may be multiple return data.frames intermingled!)
  if( fnc.ret.vals == 1 )
  {
    # Easy: just one!
    ret.val <- data.table::rbindlist(tmp);
  } else
  {
    # Combine them in turn:
    ret.val <- lapply(1:fnc.ret.vals, function(i)
    {
      x <- data.table::rbindlist(lapply(tmp, function(x) x[[i]]));
    });
    if( length(tmp) > 0 ) names(ret.val) <- names(tmp[[1]]);
  }
  return (ret.val);
}
