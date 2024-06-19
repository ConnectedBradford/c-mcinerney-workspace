.apply.medication.groups <- function(medication.groups, # the medication groups
                                     data, # the data on which to apply them
                                     suppress.warnings=FALSE)
{
  if( is.null(medication.groups) || length(medication.groups) == 0 )
  {
    # by definition, there's nothing wrong with NULL
    return (list("groups"=NULL,
                 "obs_in_groups"=NULL,
                 "errors"=NULL));
  }

  # Basic checks:
  if( is.null(data) || !inherits(data, "data.frame") || nrow(data) == 0 )
  {
    if( !suppress.warnings ) .report.ewms("The data must be a non-emtpy data.frame (or derived) object!\n", "error", ".apply.medication.groups", "AdhereR");
    return (NULL);
  }
  data <- as.data.frame(data); # make sure we convert it to a data.frame

  if( !(is.character(medication.groups) || is.factor(medication.groups)) )
  {
    if( !suppress.warnings ) .report.ewms("The medication groups must be a vector of characters!\n", "error", ".apply.medication.groups", "AdhereR");
    return (NULL);
  }

  if( length(medication.groups) == 1 && (medication.groups %in% names(data)) )
  {
    # It is a column in the data: transform it into the corresponding explicit definitions:
    mg.vals <- unique(data[,medication.groups]); mg.vals <- mg.vals[!is.na(mg.vals)]; # the unique non-NA values
    if( is.null(mg.vals) || length(mg.vals) == 0 )
    {
      if( !suppress.warnings ) .report.ewms("The column '",medication.groups,"' in the data must contain at least one non-missing value!\n", "error", ".apply.medication.groups", "AdhereR");
      return (NULL);
    }
    mg.vals <- sort(mg.vals); # makes it easier to view if ordered
    medication.groups.defs <- paste0('(',medication.groups,' == "',mg.vals,'")'); names(medication.groups.defs) <- mg.vals;
    medication.groups <- medication.groups.defs;
  }

  if( length(names(medication.groups)) != length(medication.groups) || any(duplicated(names(medication.groups))) || ("" %in% names(medication.groups)) )
  {
    if( !suppress.warnings ) .report.ewms("The medication groups must be a named list with unique and non-empty names!\n", "error", ".apply.medication.groups", "AdhereR");
    return (NULL);
  }


  # The safe environment for evaluating the medication groups (no parent) containing only the needed variables and functions (as per https://stackoverflow.com/a/18391779):
  safe_env <- new.env(parent = emptyenv());
  # ... add the safe functions:
  for( .f in c(getGroupMembers("Math"),
               getGroupMembers("Arith"),
               getGroupMembers("Logic"),
               getGroupMembers("Compare"),
               "(", "[", "!") )
  {
    safe_env[[.f]] <- get(.f, "package:base");
  }
  # ... and variables:
  assign(".data", data, envir=safe_env); # the data
  assign(".res", matrix(NA, nrow=nrow(data), ncol=length(medication.groups)), envir=safe_env); # matrix of results from evaluating the medication classes
  assign(".evald", rep(FALSE, length(medication.groups)), envir=safe_env); # records which the medication groups were already evaluated (to speed up things)
  assign(".errs", NULL, envir=safe_env); # possible errors encountered during the evaluation
  # The safe eval function (use evalq to avoid searching in the current environment): evalq(expr, env = safe_env);

  # Check, transform, parse and add the group definitions as functions to the .mg_safe_env_original environment to be later used in the safe evaluation environment:
  mg <- data.frame("name"=names(medication.groups),
                   "def"=medication.groups,
                   "uid"=make.names(names(medication.groups), unique=TRUE, allow_=TRUE),
                   "fnc_name"=NA); # convert the names to valid identifiers
  mg$fnc_name <- paste0(".mg_fnc_", mg$uid);
  for( i in 1:nrow(mg) )
  {
    # Cache the definition:
    s <- mg$def[i];
    if( is.na(s) || is.null(s) || s == "" )
    {
      # Empty definition!
      if( !suppress.warnings ) .report.ewms(paste0("Error parsing the medication class definition '",mg$name[i],"': this definition is empty!\n"), "error", ".apply.medication.groups", "AdhereR");
      return (NULL);
    }

    # Search for medication group references of the form {name} :
    ss <- gregexpr("\\{[^(\\})]+\\}",s); #gregexpr("\\{[[:alpha:]]+\\}",s);
    if( ss[[1]][1] != (-1) )
    {
      # There's at least one group reference: replace it by the corresponding function call:
      ss_calls <- regmatches(s, ss)[[1]];
      ss_calls_replaced <- vapply(ss_calls, function(x)
      {
        if( length(ii <- which(mg$name == substr(x, 2, nchar(x)-1))) != 1 )
        {
          if( !suppress.warnings ) .report.ewms(paste0("Error parsing the medication class definition '",mg$name[i],"': there is a call to the undefined medication class '",substr(x, 2, nchar(x)-1),"'!\n"), "error", ".apply.medication.groups", "AdhereR");
          return (NA);
        } else
        {
          return (paste0("safe_env$",mg$fnc_name[ii],"()"));
        }
      }, character(1));
      if( anyNA(ss_calls_replaced) )
      {
        # Error parsing medication class definitions (the message should be already generated):
        return (NULL);
      }
      regmatches(s, ss)[[1]] <- ss_calls_replaced;
    }

    # Make it into a function definition:
    s_fnc <- paste0("function()
                    {
                      if( !safe_env$.evald[",i,"] )
                      {
                        # not already evaluated:
                        tmp <- try(with(safe_env$.data, ",s,"), silent=TRUE);
                        if( inherits(tmp, 'try-error') )
                        {
                          # fail gracefully and informatively:
                          safe_env$.errs <- rbind(safe_env$.errs,
                                                  data.frame('fnc'='",mg$fnc_name[i],"', 'error'=as.character(tmp)));
                          stop('Error in ",mg$fnc_name[i],"');
                        }
                        # sanity checks:
                        if( is.null(tmp) || !is.logical(tmp) || length(tmp) != nrow(safe_env$.data) )
                        {
                          # fail gracefully and informatively:
                          safe_env$.errs <- rbind(safe_env$.errs,
                                                  data.frame('fnc'='",mg$fnc_name[i],"', 'error'='Error: evaluation did not produce logical results'));
                          stop('Error in ",mg$fnc_name[i],"');
                        }
                        safe_env$.res[,",i,"] <- tmp; safe_env$.evald[",i,"] <- TRUE;
                      } # else, it should have already been evaluated!
                      # return the value
                      return (safe_env$.res[,",i,"]);
                    }");

    # Parse it:
    s_parsed <- NULL;
    try(s_parsed <- parse(text=s_fnc), silent=TRUE);
    if( is.null(s_parsed) )
    {
      if( !suppress.warnings ) .report.ewms(paste0("Error parsing the medication class definition '",mg$name[i],"'!\n"), "error", ".apply.medication.groups", "AdhereR");
      return (NULL);
    }

    # Add the function:
    safe_env[[mg$fnc_name[i]]] <- eval(s_parsed);
  }


  # Evaluate the medication groups on the data:
  for( i in 1:nrow(mg) )
  {
    # Call the corresponding function:
    res <- NULL;
    try(res <- eval(parse(text=paste0(mg$fnc_name[i],"()")), envir=safe_env), silent=TRUE); # the errors will be anyway recorded in the .errs variable in safe_env
    if( is.null(res) || inherits(res, 'try-error') || !safe_env$.evald[i] )
    {
      # Evaluation error:
      if( !is.null(safe_env$.errs) )
      {
        err_msgs <- unique(safe_env$.errs);
        err_msgs$error <- vapply(err_msgs$error, function(s) trimws(strsplit(s,":",fixed=TRUE)[[1]][2]), character(1));
        err_msgs$fnc <- vapply(err_msgs$fnc, function(s) mg$name[ mg$fnc_name == s ], character(1));
        ss <- gregexpr(".mg_fnc_",err_msgs$error); regmatches(err_msgs$error, ss) <- "";
        if( !suppress.warnings ) .report.ewms(paste0("Error(s) during the evaluation of medication class(es): ",paste0("'",err_msgs$error," (for ",err_msgs$fnc,")'",colapse=", "),"!\n"), "warning", ".apply.medication.groups", "AdhereR");
      }
      return (NULL);
    }
  }

  # Retrieve the results and compute __ALL_OTHERS__:
  groups_info <- cbind(rbind(mg[,c("name", "def")], c("__ALL_OTHERS__", "*all observations not included in the defined groups*")),
                       "evaluated"=c(safe_env$.evald, TRUE)); rownames(groups_info) <- NULL;
  obs_in_groups <- cbind(safe_env$.res, rowSums(!is.na(safe_env$.res) & safe_env$.res) == 0); colnames(obs_in_groups) <- c(mg$name, "__ALL_OTHERS__");

  # ... and return them:
  return (list("defs"=groups_info[,c("name", "def")], "obs"=obs_in_groups));
}