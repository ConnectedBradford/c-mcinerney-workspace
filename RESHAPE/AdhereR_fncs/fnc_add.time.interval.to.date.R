# Auxiliary function: add time units to date:
.add.time.interval.to.date <- function( start.date, # a Date object
                                        time.interval=0, # the number of "units" to add to the start.date (must be numeric and is rounded)
                                        unit="days", # can be "days", "weeks", "months", "years"
                                        suppress.warnings=FALSE
)
{
  # Checks
  if( !inherits(start.date,"Date") )
  {
    if( !suppress.warnings ) .report.ewms("Parameter start.date of .add.time.interval.to.date() must be a Date() object.\n", "error", ".add.time.interval.to.date", "AdhereR");
    return (NA);
  }
  if( !is.numeric(time.interval) )
  {
    if( inherits(time.interval, "difftime") ) # check if a difftime and, if so, attempt conversion to number of units
    {
      # Time difference:
      if( units(time.interval) == as.character(unit) )
      {
        time.interval <- as.numeric(time.interval);
      } else
      {
        # Try to convert it:
        time.interval <- as.numeric(time.interval, unit="days");
        time.interval <- switch( as.character(unit),
                                 "days"  = time.interval,
                                 "weeks" = time.interval/7,
                                 "months" = {if( !suppress.warnings ) .report.ewms("Converting to months assuming a 30-days month!", "warning", ".add.time.interval.to.date", "AdhereR"); time.interval/30;},
                                 "years"  = {if( !suppress.warnings ) .report.ewms("Converting to years assuming a 365-days year!", "warning", ".add.time.interval.to.date", "AdhereR"); time.interval/365;},
                                 {if( !suppress.warnings ) .report.ewms(paste0("Unknown unit '",unit,"' to '.add.time.interval.to.date'.\n"), "error", ".add.time.interval.to.date", "AdhereR"); return(NA);} # default
        );
      }
    } else
    {
      if( !suppress.warnings ) .report.ewms("Parameter start.date of .add.time.interval.to.date() must be a number.\n", "error", ".add.time.interval.to.date", "AdhereR");
      return (NA);
    }
  }
    
  # time.interval <- round(time.interval);
  # return (switch( as.character(unit),
  #                 "days"  = (start.date + time.interval),
  #                 "weeks" = (start.date + time.interval*7),
  #                 "months" = {tmp <- (start.date + months(time.interval));
  #                             i <- which(is.na(tmp));
  #                             if( length(i) > 0 ) tmp[i] <- start.date[i] + lubridate::days(1) + months(time.interval);
  #                             tmp;},
  #                 "years"  = {tmp <- (start.date + lubridate::years(time.interval));
  #                             i <- which(is.na(tmp));
  #                             if( length(i) > 0 ) tmp[i] <- start.date[i] + lubridate::days(1) + lubridate::years(time.interval);
  #                             tmp;},
  #                 {if( !suppress.warnings ) .report.ewms(paste0("Unknown unit '",unit,"' to '.add.time.interval.to.date'.\n"), "error", ".add.time.interval.to.date", "AdhereR"); NA;} # default
  # ));

  # Faster but assumes that the internal representation of "Date" object is in number of days since the begining of time (probably stably true):
  return (switch( as.character(unit),
                            "days"  = structure(unclass(start.date) + round(time.interval), class="Date"),
                            "weeks" = structure(unclass(start.date) + round(time.interval*7), class="Date"),
                            "months" = lubridate::add_with_rollback(start.date,
                                                                    lubridate::period(round(time.interval),"months"),
                                                                    roll_to_first=TRUE), # take care of cases such as 2001/01/29 + 1 month
                            "years"  = lubridate::add_with_rollback(start.date,
                                                                    lubridate::period(round(time.interval),"years"),
                                                                    roll_to_first=TRUE), # take care of cases such as 2000/02/29 + 1 year
                            {if( !suppress.warnings ) .report.ewms(paste0("Unknown unit '",unit,"' to '.add.time.interval.to.date'.\n"), "error", ".add.time.interval.to.date", "AdhereR"); NA;} # default
  ));
}
