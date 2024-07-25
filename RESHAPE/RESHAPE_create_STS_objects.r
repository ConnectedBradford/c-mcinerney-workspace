# The purpose of this script is...
#

#####################
# Install packages. #
#####################
if( !"pacman" %in% installed.packages() ){ install.packages( "pacman" ) }
pacman::p_load(
    bigrquery
    ,tidyverse
    ,TraMineR
    ,TraMineRextras
    )


##########################
# Check for reequisites. #
##########################
if( !exists( "df_event_factor" ) )
    {
    run_cohort_generator <- readline("It looks like `RESHAPE_format_the_data.r has not been run\nDo you want to run `RESHAPE_format_the_data.r`? [Y, N]")
    if( run_cohort_generator == "Y"){ source( 'RESHAPE_format_the_data.r' ) }
}


##################################
# Create state-sequence objects. #
##################################
#
# Create state-sequence objects for TraMineR. In the previous iteration, the Clinical Review Board requested me to
# separate the test statuses from the prescriptions in the state distribution plot because they help to answer two
# distinct questions. To achieve this, I will create separate state-sequence objects in addition to the combined one.
#
# ## Including 'Unobserved' ## #
# Define list of events to focus on.
events_selection <-
    df_event_factor %>%
    dplyr::select( event_fct_order ) %>%
    droplevels() %>%
    dplyr::pull()

# Convert dataframe into a state sequence object for `TraMineR`.
df_seq <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::mutate_at( .vars = vars( event_value ), .funs = funs( as.factor ) ) %>%
    dplyr::filter( event_value %in% events_selection ) %>%
    dplyr::group_by( person_id ) %>%
    dplyr::mutate( rn = row_number() ) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    suppressWarnings()

# Create the transition matrix.
# In previous iterations, I used the `TraMineRextras::seqe2stm` function, but it requires a matrix with 2^nevents rows,
# where nevents = the count of events in `df_event_factor`. In this iterations, we have 48 events, which makes for an
# untenable matrix with 562,949,953,421,312 rows. Instead, I wrote the syntax below to create the transition matrix
# explicitly, via a loop rather than using matrix operations.
stm <- matrix(
        rep( events_selection, each = length( events_selection ) )
        ,ncol = length( events_selection )
        ) %>%
        `colnames<-`( events_selection ) %>%
        `rownames<-`( events_selection ) %>%
        rbind( data.frame( None = events_selection ) %>% t(), .)

# Convert the time series dataframe to a state sequence object.
# Similarly to the previous action, I have had to make some changes to the syntax from the previous iteration.
d <- TraMineRextras::TSE_to_STS(
        df_seq
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "event_value"
        ,stm = stm
        ,tmax = df_seq %>% dplyr::select( rn ) %>% dplyr::ungroup() %>% max()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine <-
    d %>%
    TraMineR::seqdef( labels = events_selection[ events_selection %in% statl] ) %>%
    suppressMessages()

# ## Excluding 'Unobserved' ## #
# Define list of events to focus on.
events_selection <-
    df_event_factor %>%
    dplyr::select( event_fct_order ) %>%
     dplyr::filter( event_fct_order != "Unobserved" ) %>%
    droplevels() %>%
    dplyr::pull()

# Convert dataframe into a state sequence object for `TraMineR`.
df_seq_excludingUnobserved <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::mutate_at( .vars = vars( event_value ), .funs = funs( as.factor ) ) %>%
    dplyr::filter( event_value %in% events_selection ) %>%
    dplyr::group_by( person_id ) %>%
    dplyr::mutate( rn = row_number() ) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    suppressWarnings()

# Create the transition matrix.
# In previous iterations, I used the `TraMineRextras::seqe2stm` function, but it requires a matrix with 2^nevents rows,
# where nevents = the count of events in `df_event_factor`. In this iterations, we have 48 events, which makes for an
# untenable matrix with 562,949,953,421,312 rows. Instead, I wrote the syntax below to create the transition matrix
# explicitly, via a loop rather than using matrix operations.
stm_excludingUnobserved <- matrix(
        rep( events_selection, each = length( events_selection ) )
        ,ncol = length( events_selection )
        ) %>%
        `colnames<-`( events_selection ) %>%
        `rownames<-`( events_selection ) %>%
        rbind( data.frame( None = events_selection ) %>% t(), .)

# Convert the time series dataframe to a state sequence object.
# Similarly to the previous action, I have had to make some changes to the syntax from the previous iteration.
d <- TraMineRextras::TSE_to_STS(
        df_seq_excludingUnobserved
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "event_value"
        ,stm = stm_excludingUnobserved
        ,tmax = df_seq_excludingUnobserved %>% dplyr::select( rn ) %>% dplyr::ungroup() %>% max()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_excludingUnobserved <-
    d %>%
    TraMineR::seqdef( labels = events_selection[ events_selection %in% statl] ) %>%
    suppressMessages()


# Make an test-only state-sequence object for some particular plots.
#
# I actually make two. One contains the "Unobserved" state, which is useful for tracking when sequences stop. The
# other exlcudes the "Unobserved" state, which is useful for plotting proportions of the remaining states without
# being distracted by the growing proportion of unobserved events.
#
# Define list of events to focus on.
events_selection <-
    df_event_factor %>%
    dplyr::select( event_fct_order ) %>%
    dplyr::filter( stringr::str_detect( event_fct_order, pattern = "(Test)" ) | event_fct_order == "Unobserved" ) %>%
    droplevels() %>%
    dplyr::pull()

# Convert dataframe into a state sequence object for `TraMineR`.
df_seq_test_only <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::mutate_at( .vars = vars( event_value ), .funs = funs( as.factor ) ) %>%
    dplyr::filter( event_value %in% events_selection ) %>%
    dplyr::group_by( person_id ) %>%
    dplyr::mutate( rn = row_number() ) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    suppressWarnings()

# Create the transition matrix.
# In previous iterations, I used the `TraMineRextras::seqe2stm` function, but it requires a matrix with 2^nevents rows,
# where nevents = the count of events in `df_event_factor`. In this iterations, we have 48 events, which makes for an
# untenable matrix with 562,949,953,421,312 rows. Instead, I wrote the syntax below to create the transition matrix
# explicitly, via a loop rather than using matrix operations.
stm_test_only <-
    matrix(
        rep( events_selection, each = length( events_selection ) )
        ,ncol = length( events_selection )
        ) %>%
    `colnames<-`( events_selection ) %>%
    `rownames<-`( events_selection ) %>%
    rbind( data.frame( None = events_selection ) %>% t(), .)

# Convert the time series dataframe to a state sequence object.
# Similarly to the previous action, I have had to make some changes to the syntax from the previous iteration.
d <- TraMineRextras::TSE_to_STS(
        df_seq_test_only
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "event_value"
        ,stm = stm_test_only
        ,tmax = df_seq_test_only %>% dplyr::select( rn ) %>% dplyr::ungroup() %>% max()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_test_only <-
    d %>%
    TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
    suppressMessages()

# Define list of events to focus on.
events_selection <-
    df_event_factor %>%
    dplyr::select( event_fct_order ) %>%
    dplyr::filter( stringr::str_detect( event_fct_order, pattern = "(Test)" ) ) %>%
    droplevels() %>%
    dplyr::pull()

# Convert dataframe into a state sequence object for `TraMineR`.
df_seq_test_only_excludingUnobserved <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::mutate_at( .vars = vars( event_value ), .funs = funs( as.factor ) ) %>%
    dplyr::filter( event_value %in% events_selection ) %>%
    dplyr::group_by( person_id ) %>%
    dplyr::mutate( rn = row_number() ) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    suppressWarnings()

# Create the transition matrix.
# In previous iterations, I used the `TraMineRextras::seqe2stm` function, but it requires a matrix with 2^nevents rows,
# where nevents = the count of events in `df_event_factor`. In this iterations, we have 48 events, which makes for an
# untenable matrix with 562,949,953,421,312 rows. Instead, I wrote the syntax below to create the transition matrix
# explicitly, via a loop rather than using matrix operations.
stm_test_only_excludingUnobserved <-
    matrix(
        rep( events_selection, each = length( events_selection ) )
        ,ncol = length( events_selection )
        ) %>%
    `colnames<-`( events_selection ) %>%
    `rownames<-`( events_selection ) %>%
    rbind( data.frame( None = events_selection ) %>% t(), .)

# Convert the time series dataframe to a state sequence object.
# Similarly to the previous action, I have had to make some changes to the syntax from the previous iteration.
d <- TraMineRextras::TSE_to_STS(
        df_seq_test_only_excludingUnobserved
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "event_value"
        ,stm = stm_test_only_excludingUnobserved
        ,tmax = df_seq_test_only_excludingUnobserved %>% dplyr::select( rn ) %>% dplyr::ungroup() %>% max()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_test_only_excludingUnobserved <-
    d %>%
    TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
    suppressMessages()


# Make an intervention-only state-sequence object for some particular plots.
#
# Like the test-only objects, I actually make two. One contains the "Unobserved" state, which is useful for tracking
# when sequences stop. The other exlcudes the "Unobserved" state, which is useful for plotting proportions of the remaining
# states without being distracted by the growing proportion of unobserved events.
#
# Define list of events to focus on.
events_selection <-
    df_event_factor %>%
    dplyr::select( event_fct_order ) %>%
    dplyr::filter( !stringr::str_detect( event_fct_order, pattern = "(Test)" ) | event_fct_order == "Unobserved" ) %>%
    droplevels() %>%
    dplyr::pull()

# Convert dataframe into a state sequence object for `TraMineR`.
df_seq_intervention_only <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::mutate_at( .vars = vars( event_value ), .funs = funs( as.factor ) ) %>%
    dplyr::filter( event_value %in% events_selection ) %>%
    dplyr::group_by( person_id ) %>%
    dplyr::mutate( rn = row_number() ) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    suppressWarnings()

# Create the transition matrix.
# In previous iterations, I used the `TraMineRextras::seqe2stm` function, but it requires a matrix with 2^nevents rows,
# where nevents = the count of events in `df_event_factor`. In this iterations, we have 48 events, which makes for an
# untenable matrix with 562,949,953,421,312 rows. Instead, I wrote the syntax below to create the transition matrix
# explicitly, via a loop rather than using matrix operations.
stm_intervention_only <-
    matrix(
        rep( events_selection, each = length( events_selection ) )
        ,ncol = length( events_selection )
        ) %>%
    `colnames<-`( events_selection ) %>%
    `rownames<-`( events_selection ) %>%
    rbind( data.frame( None = events_selection ) %>% t(), .)

# Convert the time series dataframe to a state sequence object.
# Similarly to the previous action, I have had to make some changes to the syntax from the previous iteration.
d <- TraMineRextras::TSE_to_STS(
        df_seq_intervention_only
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "event_value"
        ,stm = stm_intervention_only
        ,tmax = df_seq_intervention_only %>% dplyr::select( rn ) %>% dplyr::ungroup() %>% max()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_intervention_only <-
    d %>%
    TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
    suppressMessages()

# Define list of events to focus on.
events_selection <-
    df_event_factor %>%
    dplyr::select( event_fct_order ) %>%
    dplyr::filter( !stringr::str_detect( event_fct_order, pattern = "Test" ) & event_fct_order != "Unobserved" ) %>%
    droplevels() %>%
    dplyr::pull()

# Convert dataframe into a state sequence object for `TraMineR`.
df_seq_intervention_only_excludingUnobserved <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::mutate_at( .vars = vars( event_value ), .funs = funs( as.factor ) ) %>%
    dplyr::filter( event_value %in% events_selection ) %>%
    dplyr::group_by( person_id ) %>%
    dplyr::mutate( rn = row_number() ) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    suppressWarnings()

# Create the transition matrix.
# In previous iterations, I used the `TraMineRextras::seqe2stm` function, but it requires a matrix with 2^nevents rows,
# where nevents = the count of events in `df_event_factor`. In this iterations, we have 48 events, which makes for an
# untenable matrix with 562,949,953,421,312 rows. Instead, I wrote the syntax below to create the transition matrix
# explicitly, via a loop rather than using matrix operations.
stm_intervention_only_excludingUnobserved <-
    matrix(
        rep( events_selection, each = length( events_selection ) )
        ,ncol = length( events_selection )
        ) %>%
    `colnames<-`( events_selection ) %>%
    `rownames<-`( events_selection ) %>%
    rbind( data.frame( None = events_selection ) %>% t(), .)

# Convert the time series dataframe to a state sequence object.
# Similarly to the previous action, I have had to make some changes to the syntax from the previous iteration.
d <- TraMineRextras::TSE_to_STS(
        df_seq_intervention_only_excludingUnobserved
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "event_value"
        ,stm = stm_intervention_only_excludingUnobserved
        ,tmax = df_seq_intervention_only_excludingUnobserved %>% dplyr::select( rn ) %>% dplyr::ungroup() %>% max()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_intervention_only_excludingUnobserved <-
    d %>%
    TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
    suppressMessages()


# Make an test-and-intervention status state-sequence object for some particular plots.
#
# Like the test-only objects, I actually make two. One contains the "Unobserved" state, which is useful for tracking
# when sequences stop. The other exlcudes the "Unobserved" state, which is useful for plotting proportions of the
# remaining states without being distracted by the growing proportion of unobserved events.
#
# Define list of events to focus on.
events_selection <-
    df_TandI_factor %>%
    dplyr::select( TandI_fct_order ) %>%
    droplevels() %>%
    dplyr::pull()

# Create the data transition matrix.
stm_TandI <-
    matrix(
        rep( events_selection, each = length( events_selection ) )
        ,ncol = length( events_selection )
        ) %>%
    `colnames<-`( events_selection ) %>%
    `rownames<-`( events_selection ) %>%
    rbind( data.frame( None = events_selection ) %>% t(), .)

# Convert the time series dataframe to a state sequence object.
d <-
    TraMineRextras::TSE_to_STS(
        df_seq_test_only
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "TandI"
        ,stm = stm_TandI
        ,tmax = df_seq_test_only %>% dplyr::summarise( max( rn ) ) %>% dplyr::pull()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_TandI <-
   d %>%
   TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
   suppressMessages()

# Convert the time series dataframe to a state sequence object.
d <-
    TraMineRextras::TSE_to_STS(
        df_seq_test_only_excludingUnobserved
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "TandI"
        ,stm = stm_TandI
        ,tmax = df_seq_test_only %>% dplyr::summarise( max( rn ) ) %>% dplyr::pull()
    ) %>%
    dplyr::select( - a1 )
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_TandI_excludingUnobserved <-
   d %>%
   TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
   suppressMessages()


# Make an HMA state state-sequence object for some particular plots.
#
# Like the test-only objects, I actually make two. One contains the "Unobserved" state, which is useful
# for tracking when sequences stop. The other exlcudes the "Unobserved" state, which is useful for plotting
# proportions of the remaining states without being distracted by the growing proportion of unobserved events.
#
# Define list of events to focus on.
events_selection <-
    df_HMA_factor %>%
    dplyr::select( HMA_fct_order ) %>%
    droplevels() %>%
    dplyr::pull()

# Create the data transition matrix.
stm_HMA <-
    matrix(
        rep( events_selection, each = length( events_selection ) )
        ,ncol = length( events_selection )
        ) %>%
    `colnames<-`( events_selection ) %>%
    `rownames<-`( events_selection ) %>%
    rbind( data.frame( None = events_selection ) %>% t(), .)

# Convert the time series dataframe to a state sequence object.
seqdata_HMA <-
    df_seq_test_only %>%
    dplyr::filter( !is.na( HMA ) ) %>%
    dplyr::group_by( person_id ) %>%
    dplyr::mutate( rn = row_number() ) %>% # Need to reassign the row number.
    dplyr::ungroup() %>%
    as.data.frame()
d <-
    TraMineRextras::TSE_to_STS(
        seqdata_HMA
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "HMA"
        ,stm = stm_HMA
        ,tmax = seqdata_HMA %>% dplyr::summarise( max( rn ) ) %>% dplyr::pull()
    ) %>%
     dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_HMA <-
   d %>%
   TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
   suppressMessages()

# Convert the time series dataframe to a state sequence object.
seqdata_HMA_excludingUnobserved <-
    df_seq_test_only_excludingUnobserved %>%
    dplyr::filter( !is.na( HMA ) ) %>%
    dplyr::group_by( person_id ) %>%
    dplyr::mutate( rn = row_number() ) %>% # Need to reassign the row number.
    dplyr::ungroup() %>%
    as.data.frame()
d <-
    TraMineRextras::TSE_to_STS(
        seqdata_HMA_excludingUnobserved
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "HMA"
        ,stm = stm_HMA
        ,tmax = seqdata_HMA_excludingUnobserved %>% dplyr::summarise( max( rn ) ) %>% dplyr::pull()
    ) %>%
     dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_HMA_excludingUnobserved <-
   d %>%
   TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
   suppressMessages()


# Make an HMA-and-Test status state-sequence object for some particular plots.
#
# Like the test-only objects, I actually make two. One contains the "Unobserved" state, which is useful for tracking
# when sequences stop. The other exlcudes the "Unobserved" state, which is useful for plotting proportions of the
# remaining states without being distracted by the growing proportion of unobserved events.
#
# Define list of events to focus on.
events_selection <-
    df_HMAandTestStatus_factor %>%
    dplyr::select( HMAandTestStatus_fct_order ) %>%
    droplevels() %>%
    dplyr::pull()

# Create the data transition matrix.
stm_HMAandTestStatus <-
    matrix(
        rep( events_selection, each = length( events_selection ) )
        ,ncol = length( events_selection )
        ) %>%
    `colnames<-`( events_selection ) %>%
    `rownames<-`( events_selection ) %>%
    rbind( data.frame( None = events_selection ) %>% t(), .)

# Convert the time series dataframe to a state sequence object.    
d <-
    TraMineRextras::TSE_to_STS(
        seqdata_HMA
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "HMAandTestStatus"
        ,stm = stm_HMAandTestStatus
        ,tmax = seqdata_HMA %>% dplyr::summarise( max( rn ) ) %>% dplyr::pull()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_HMAandTestStatus <-
   d %>%
   TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
   suppressMessages()

# Convert the time series dataframe to a state sequence object.    
d <-
    TraMineRextras::TSE_to_STS(
        seqdata_HMA_excludingUnobserved
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "HMAandTestStatus"
        ,stm = stm_HMAandTestStatus
        ,tmax = seqdata_HMA_excludingUnobserved %>% dplyr::summarise( max( rn ) ) %>% dplyr::pull()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_HMAandTestStatus_excludingUnobserved <-
   d %>%
   TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
   suppressMessages()


# Make a test-and-multimorbidity status state-sequence object for some particular plots.
#
# Like the test-only objects, I actually make two. One contains the "Unobserved" state, which is useful for
# tracking when sequences stop. The other exlcudes the "Unobserved" state, which is useful for plotting
# proportions of the remaining states without being distracted by the growing proportion of unobserved events.
#
# Define list of events to focus on.
events_selection <-
    df_TandMultiMorb_factor %>%
    dplyr::select( TandMultiMorb_fct_order ) %>%
    droplevels() %>%
    dplyr::pull()

# Create the data transition matrix.
stm_TandMultiMorb <-
    matrix(
        rep( events_selection, each = length( events_selection ) )
        ,ncol = length( events_selection )
        ) %>%
    `colnames<-`( events_selection ) %>%
    `rownames<-`( events_selection ) %>%
    rbind( data.frame( None = events_selection ) %>% t(), .)

# Convert the time series dataframe to a state sequence object.
d <-
    TraMineRextras::TSE_to_STS(
        df_seq_test_only
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "TandMultiMorb"
        ,stm = stm_TandMultiMorb
        ,tmax = df_seq_test_only %>% dplyr::summarise( max( rn ) ) %>% dplyr::pull()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_TandMultiMorb <-
   d %>%
   TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
   suppressMessages()

# Convert the time series dataframe to a state sequence object.
d <-
    TraMineRextras::TSE_to_STS(
        df_seq_test_only_excludingUnobserved
        ,id = "person_id"
        ,timestamp = "rn"
        ,event = "TandMultiMorb"
        ,stm = stm_TandMultiMorb
        ,tmax = df_seq_test_only %>% dplyr::summarise( max( rn ) ) %>% dplyr::pull()
    ) %>%
    dplyr::select( - a1 ) 
statl <- TraMineR:::seqxtract(d, NULL, data.frame = TRUE) %>% seqstatl()
sts_seqmine_TandMultiMorb_excludingUnobserved <-
   d %>%
   TraMineR::seqdef( labels = events_selection[ events_selection %in% statl ] ) %>%
   suppressMessages()
