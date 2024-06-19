# Auxiliary function: subtract two dates to obtain the number of days in between:
# WARNING! Faster than difftime() but makes the assumption that the internal representation of Date objects is the number of days since a given beginning of time
# (true as of R 3.4 and probably conserved in future versions)
.difftime.Dates.as.days <- function( start.dates, end.dates, suppress.warnings=FALSE )
{
  # Checks
  if( !inherits(start.dates,"Date") || !inherits(end.dates,"Date") )
  {
    if( !suppress.warnings ) .report.ewms("start.dates and end.dates to '.difftime.Dates.as.days' must be a Date() objects.\n", "error", ".difftime.Dates.as.days", "AdhereR");
    return (NA);
  }
  return (unclass(start.dates) - unclass(end.dates)); # the difference between the raw internal representations of Date objects is in days
}