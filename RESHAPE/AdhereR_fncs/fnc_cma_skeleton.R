.cma.skeleton <- function(data=NULL,
                          ret.val=NULL,
                          cma.class.name=NA,
                          ID.colname=NA,
                          event.date.colname=NA,
                          event.duration.colname=NA,
                          event.daily.dose.colname=NA,
                          medication.class.colname=NA,
                          medication.groups.colname=NA,
                          flatten.medication.groups=NA,
                          followup.window.start.per.medication.group=NA,
                          event.interval.colname=NA,
                          gap.days.colname=NA,
                          carryover.within.obs.window=NA,
                          carryover.into.obs.window=NA,
                          carry.only.for.same.medication=NA,
                          consider.dosage.change=NA,
                          followup.window.start=NA,
                          followup.window.start.unit=NA,
                          followup.window.duration=NA,
                          followup.window.duration.unit=NA,
                          observation.window.start=NA,
                          observation.window.start.unit=NA,
                          observation.window.duration=NA,
                          observation.window.duration.unit=NA,
                          date.format=NA,
                          suppress.warnings=FALSE,
                          suppress.special.argument.checks=FALSE,
                          force.NA.CMA.for.failed.patients=TRUE,
                          parallel.backend="none",
                          parallel.threads=1,
                          .workhorse.function=NULL)
{
  # Convert to data.table, cache event date as Date objects keeping only the necessary columns, check for column name conflicts, and key by patient ID and event date:
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
  if( is.null(mg) )
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

    # Convert to data.frame and return:
    if( force.NA.CMA.for.failed.patients )
    {
      # Make sure patients with failed CMA estimations get an NA estimate!
      patids <- unique(data.copy[,get(ID.colname)]);
      if( length(patids) > nrow(tmp$CMA) )
      {
        setnames(tmp$CMA, 1, ".ID"); tmp$CMA <- merge(data.table(".ID"=patids, key=".ID"), tmp$CMA, all.x=TRUE);
      }
    }
    setnames(tmp$CMA, c(ID.colname,"CMA")); ret.val[["CMA"]] <- as.data.frame(tmp$CMA);
    ret.val[["event.info"]] <- as.data.frame(tmp$event.info);
    class(ret.val) <- c(cma.class.name, class(ret.val));
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
      .report.ewms(paste0("None of the medication classes (included __ALL_OTHERS__) selects any observation!\n"), "warning", cma.class.name, "AdhereR");
      return (NULL);
    }
    mb.obs <- mb.obs[,mg.to.eval]; # keep only the non-trivial ones

    # How is the FUW to be estimated?
    if( !followup.window.start.per.medication.group )
    {
      # The FUW and OW are estimated once per patient (i.e., all medication groups share the same FUW and OW):
      # Call the compute.event.int.gaps() function and use the results:
      event.info <- compute.event.int.gaps(data=as.data.frame(data.copy),
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
                                           return.data.table=FALSE);
      if( is.null(event.info) ) return (list("CMA"=NA, "event.info"=NULL));

      # Add the FUW and OW start dates to the data:
      data.copy <- merge(data.copy, unique(event.info[,c(ID.colname, ".FU.START.DATE", ".OBS.START.DATE")]), by=ID.colname, all.x=TRUE, all.y=FALSE);
      names(data.copy)[ names(data.copy) == ".FU.START.DATE" ]  <- ".FU.START.DATE.PER.PATIENT.ACROSS.MGS";
      names(data.copy)[ names(data.copy) == ".OBS.START.DATE" ] <- ".OBS.START.DATE.PER.PATIENT.ACROSS.MGS";

      # Adjust the corresponding params:
      actual.followup.window.start         <- ".FU.START.DATE.PER.PATIENT.ACROSS.MGS";
      actual.followup.window.start.unit    <- NA;
      actual.observation.window.start      <- ".OBS.START.DATE.PER.PATIENT.ACROSS.MGS";
      actual.observation.window.start.unit <- NA;
    } else
    {
      # The FUW and OW are estimated separately for each medication group -- nothing to do...
      actual.followup.window.start         <- followup.window.start;
      actual.followup.window.start.unit    <- followup.window.start.unit;
      actual.observation.window.start      <- observation.window.start;
      actual.observation.window.start.unit <- observation.window.start.unit;
    }

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
                               followup.window.start=actual.followup.window.start,
                               followup.window.start.unit=actual.followup.window.start.unit,
                               followup.window.duration=followup.window.duration,
                               followup.window.duration.unit=followup.window.duration.unit,
                               observation.window.start=actual.observation.window.start,
                               observation.window.start.unit=actual.observation.window.start.unit,
                               observation.window.duration=observation.window.duration,
                               observation.window.duration.unit=observation.window.duration.unit,
                               date.format=date.format,
                               suppress.warnings=suppress.warnings,
                               suppress.special.argument.checks=suppress.special.argument.checks);
      if( is.null(tmp) || (!is.null(tmp$CMA) && !inherits(tmp$CMA,"data.frame")) || is.null(tmp$event.info) ) return (NULL);

      if( !is.null(tmp$CMA) )
      {
        # Convert to data.frame and return:
        if( force.NA.CMA.for.failed.patients )
        {
          # Make sure patients with failed CMA estimations get an NA estimate!
          patids <- unique(data.copy[,get(ID.colname)]);
          if( length(patids) > nrow(tmp$CMA) )
          {
            setnames(tmp$CMA, 1, ".ID"); tmp$CMA <- merge(data.table(".ID"=patids, key=".ID"), tmp$CMA, all.x=TRUE);
          }
        }

        setnames(tmp$CMA, c(ID.colname,"CMA")); tmp$CMA <- as.data.frame(tmp$CMA);
      }
      if( !is.null(tmp$event.info) )
      {
        tmp$event.info <- as.data.frame(tmp$event.info);
      }
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
    }
    class(ret.val) <- c(cma.class.name, class(ret.val));
    return (ret.val);

  }
}