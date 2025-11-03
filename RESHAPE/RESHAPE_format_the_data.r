# The purpose of this script is...
#


#####################
# Install packages. #F
#####################
if( !"pacman" %in% installed.packages() ){ install.packages( "pacman" ) }
pacman::p_load(
    bigrquery
    ,tidyverse
    ,TraMineR
    ,TraMineRextras
    )


##########################
# Check for requisites. #
##########################
if( !exists( "df_log_PandT_longFormat" ) )
    {
    run_cohort_generator <- readline("It looks like `RESHAPE_cohort_generator.R has not been run\nDo you want to run `RESHAPE_cohort_generator.R`? [Y, N]")
    if( run_cohort_generator == "Y"){ source( 'RESHAPE_cohort_generator.r' ) }
}


####################
# Format the data. #
####################
# I create the simplified dataframe object from the first iteration.
df_log_PandT_longFormat_simplified <-
    df_log_PandT_longFormat %>%
    # Remove the numbers in the event name.
    dplyr::mutate_at(
        .vars = vars( event_name )
        ,.funs = ~ gsub( "_*[0-9]", "", . )
    )

# Convert to factor. # This becomes too unwieldy if I don't substantially amalgamate the drugs.
if( ( df_log_PandT_longFormat_simplified %>% dplyr::select( event_value ) %>% unique() %>% nrow() ) < 36 )
    {
    df_log_PandT_longFormat_simplified <-
        df_log_PandT_longFormat_simplified %>%
        dplyr::mutate(
            event_value = 
                    factor(
                        event_value
                        ,levels = df_event_factor %>% dplyr::select( event_fct_order ) %>% dplyr::pull()
                    )
            )
    } else { 
    message("\nNOTE: Prescriptions were not converted to factor type because there are too many.\n         Consider amalgamating them.")
    }

# Before I can begin analysing the data, I need to define the strata proposed by the Clinical Review Board.
# The stratifications were:
# - H.M.A.: Four strata defined by combinations of {'Expected', 'Shorter-than-expected'} testing intervals
#   and {'No observed change', 'Observed change'} in prescriptions. The H.M.A. acronym derives from the three
#   strata: (0,0)-Hold; (1,0)-Monitor; (0 or 1, 1)-Adjust.
# - Tests-and-Interventions: Twelve strata defined by combinations of the test statuses and whether the patient
#   is on one, two or three medications simulatneously.
# - Multimorbidity: Two strata defined by whether or not there is a record for at least one of the
#   multimorbidity diagnositic codes.
#
# The first thing I do is to add a variable that indicates the patient-specific test interval. This will be
# handy for bounding the variables I need to create.

# Append an age column.
# ## For some unknown reason, I have to reload `r_tbl_srpatient` or the `DateBirth` field isn't recognised in functions like `colnames()` and `dplyr::distinct`
# ## despite clearly showing up when `r_tbl_srpatient` is printed to screen. Figure that one out!
r_tbl_srpatient <- dplyr::tbl( con, sql( paste0( "SELECT * FROM ", project_id, ".cb_SRPatient" ) ) )
df_age <-
    r_tbl_srpatient %>%
    # Select only those patients in whom we are interested.
    dplyr::inner_join( qry_records_with_T2DM_diagnoses, by = join_by( person_id ) ) %>%
    dplyr::distinct( person_id, DateBirth ) %>%
    dplyr::collect() %>%
    tidyr::separate(
        col = DateBirth
        ,sep = 4
        ,into = c( "birth_year", "birth_month" )
    ) %>%
    dplyr::mutate(
        DOB = 
            dplyr::if_else(
                is.na( birth_year )
                ,NA
                ,ISOdate( year = birth_year, month = birth_month, day = 1, hour = 0 )
            )
    ) %>%
    # If there are multiple dates of death, chose the earliest one.
    dplyr::group_by( person_id ) %>%
    dplyr::mutate( DOB = min( DOB, na.rm = TRUE ) ) %>%
    dplyr::ungroup() %>% 
    tidyr::drop_na() %>%
    # Tidy up.
    dplyr::distinct( person_id, DOB ) %>%
    dplyr::arrange( person_id ) %>%
    suppressWarnings() 

df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified %>%
    dplyr::left_join(
        df_age
        ,by = join_by( person_id )
        ,relationship = 'many-to-one'
    ) %>%
    dplyr::mutate(
        Age = lubridate::interval( DOB, start_dttm ) / lubridate::years( 1 )
    )


# Include the date of diagnosis.
df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::group_by( person_id ) %>%
    dplyr::group_modify( ~ tibble::add_row( .x ), .by = person_id ) %>% 
    dplyr::left_join( qry_records_with_T2DM_diagnoses %>% collect(), by = join_by( person_id ) ) %>%
    dplyr::mutate( start_dttm = if_else( is.na( start_dttm ), date_diagnosis , start_dttm ) ) %>%
    dplyr::arrange( person_id, start_dttm, .group_by = TRUE ) %>%
    dplyr::mutate(
        end_dttm = dplyr::if_else( is.na( end_dttm ), lead( start_dttm ), end_dttm )
        ,event_name = dplyr::if_else( is.na( event_value ), "diagnosis", event_name )
    ) %>%
    dplyr::ungroup() 

                      
# Create a variable that indicates patient-specific test interval.
df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::group_by( person_id ) %>%
    dplyr::arrange( start_dttm, .by_group = TRUE ) %>%
    dplyr::mutate(
        idx_test_interval = 0 + cumsum( stringr::str_detect( event_value, pattern = "Test" ) | stringr::str_detect( event_name, pattern = "diagnosis" ) )
    ) %>% 
    dplyr::ungroup() %>%
    # Remove any calculated `idx_test_interval` for which the `event_value` is "Unobserved".
    dplyr::mutate(
        idx_test_interval =
            dplyr::if_else( event_value == "Unobserved", NA_integer_, idx_test_interval )
    ) %>%
    # Tidy up.
    dplyr::distinct() %>%
    dplyr::arrange( person_id, start_dttm )

#########################
# H.M.A. stratification #
#########################
#
# The two components of this stratification are the testing interval and the change in prescription. The testing
# interval requries me to create a variable indicating that the inter-test duration was between
# `val_testing_interval_LB` and `val_testing_interval_UB`. The change-in-prescription criterion refers to 
# whether a new medication-dose set was prescribed, where "new" means not seen in the previous three months.
#
# The first thing is to calculate the inter-test duration.
                      
df_inter_test_duration <-
    qry_log_test_longFormat %>%
    dplyr::collect() %>%
    # Include the date of diagnosis.
    dplyr::group_by( person_id ) %>%
    dplyr::group_modify( ~ tibble::add_row( .x ), .by = person_id ) %>% 
    dplyr::left_join( qry_records_with_T2DM_diagnoses %>% collect(), by = join_by( person_id ) ) %>%
    dplyr::mutate( start_dttm = if_else( is.na( start_dttm ), date_diagnosis , start_dttm ) ) %>%
    dplyr::arrange( person_id, start_dttm, .group_by = TRUE ) %>%
    dplyr::mutate(
        end_dttm = dplyr::if_else( is.na( end_dttm ), lead( start_dttm ), end_dttm )
        ,event_name = dplyr::if_else( is.na( event_value ), "diagnosis", event_name )
    ) %>%
    dplyr::ungroup() %>%
    # Calculate the duration between each test and the next.
    dplyr::mutate( inter_test_duration_cont = lubridate::interval( start_dttm, end_dttm  ) %>% `/`( months(1) ) ) %>% # So, `inter_test_duration` is the duration *to* the next rather than *since* the previous.
    # Create a discretised version of the inter-test duration variable.
    dplyr::mutate(
        inter_test_duration_discr =
            dplyr::if_else(
                inter_test_duration_cont < val_testing_interval_UB
                ,"Shorter than expected"
                ,"As expected"
            )
        ,inter_test_duration_discr = 
            dplyr::if_else(
                inter_test_duration_cont < val_testing_interval_LB
                ,"Not applicable"
                ,inter_test_duration_discr
            )
    ) %>% 
    # Make any same-day tests or anomalies = NA.
    dplyr::mutate(
        inter_test_duration_discr =
             dplyr::if_else(
                 is.na( inter_test_duration_cont )
                 ,"Not applicable"
                 ,inter_test_duration_discr
             )
    ) %>%
    # Tidy up.
    arrange( person_id, start_dttm ) %>%
    dplyr::select( c( person_id, start_dttm, inter_test_duration_cont, inter_test_duration_discr ) ) %>% 
    dplyr::distinct() %>%
    dplyr::arrange( person_id, start_dttm )


# Add the inter-test duration variables to the dataframe.
df_log_PandT_longFormat_simplified_StrataLabels <-
    # Join the data.frame with the variables indicating the inter-test duration.
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::left_join(
        df_inter_test_duration 
        ,by = join_by( person_id, start_dttm )
        ,relationship = "many-to-many"
        # The `relationship` argument is needed because a warning is raised for the situations
        # where multiple prescriptions are given on the same day.
    ) %>%
    # Fill in the values of the inter-test variables into subsequent rows until a new value is given.
    dplyr::group_by( person_id ) %>%
    dplyr::arrange( start_dttm, .group_by = TRUE ) %>%
    tidyr::fill( inter_test_duration_cont ) %>%
    tidyr::fill( inter_test_duration_discr ) %>%
    dplyr::mutate( inter_test_duration_cont = dplyr::if_else( inter_test_duration_discr == "Not applicable", NA_real_, inter_test_duration_cont ) ) %>%
    dplyr::ungroup()

# Next, I need to make the variable that indicates whether the prescription has changed `HMA_adjust_lookBack_window_in_wks` weeks previously.
if( is.null( HMA_adjust_lookBack_window ) )
{
    df_log_PandT_longFormat_simplified_StrataLabels <-
        df_log_PandT_longFormat_simplified_StrataLabels %>%
        # Create new columns that contain the test events and the prescription events.
        dplyr::select( person_id, start_dttm, event_value, idx_test_interval ) %>%
        dplyr::filter( !stringr::str_detect( event_value, pattern = "Test") ) %>%
        dplyr::filter( !stringr::str_detect( event_value, pattern = "Unobserved") ) %>%
        dplyr::mutate( meds_only = as.character( event_value ) ) %>%
        # Remove duplicated prescriptions occuring on the same day.
        dplyr::distinct( person_id, start_dttm, meds_only ) %>%
        dplyr::arrange( person_id, start_dttm ) %>%
        # Create an incremented index for prescriptions so that same-day drugs are marked as the same prescription.
        dplyr::group_by( person_id, start_dttm ) %>%
        dplyr::mutate( prescription_index = cur_group_id() ) %>%
        dplyr::ungroup() %>%
        dplyr::select( - start_dttm ) %>%
        # Create a nested list of drugs prescribed on a given day.
        tidyr::nest( .by = person_id, .key = "meds_as_nested_list" ) %>%
        dplyr::mutate(
                meds_in_current_prescription = 
                    lapply(
                        meds_as_nested_list
                        ,function(x) slider::slide_index(.x = x$meds_only, .i = x$prescription_index, .f = ~.x, .before = 0 ) %>% array()
                        )
        ) %>%
        tidyr::unnest( cols = c( meds_as_nested_list, meds_in_current_prescription ) ) %>%
        dplyr::mutate( meds_in_current_prescription = lapply( X = meds_in_current_prescription, FUN = \(x) sort( unique( unlist(x) ) ) ) ) %>%
        # Remove `meds_only` because it is no longer required.
        dplyr::select( -meds_only ) %>%
        dplyr::distinct() %>% 
        # Create a nested list of drugs prescribed in prescriptions within the look-back window. The `purrr:map2()` bit
        # removes the current prescription from the list.
        tidyr::nest( .by = person_id, .key = "meds_as_nested_list" ) %>%
        dplyr::mutate(
            meds_in_lookback_prescriptions = 
                lapply(
                    meds_as_nested_list
                    ,function(x)
                        purrr::map2(
                            slider::slide_index( .x = seq_along( x$meds_in_current_prescription ), .i = x$prescription_index, .f = ~.x, .before = HMA_adjust_lookBack_count )
                            ,seq_along( x$meds_in_current_prescription )
                            ,~x$meds_in_current_prescription[ setdiff( .x, .y ) ]
                        )
                    )
        ) %>%
        tidyr::unnest( cols = c( meds_as_nested_list, meds_in_lookback_prescriptions ) ) %>%
        # Unlist the drugs that are included in the look-back window, and keep and sort only unique values.
        dplyr::mutate(
            meds_in_lookback_prescriptions = lapply( X = meds_in_lookback_prescriptions, FUN = \(x) sort( unique( unlist(x) ) ) )
        ) %>%
        as.data.frame() %>% 
        # Alert if there is anything in the current prescription that was not in the look-back prescriptions.
        dplyr::mutate(
            HMA = purrr::map2( meds_in_current_prescription, meds_in_lookback_prescriptions, setdiff )
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate( HMA = dplyr::if_else( !all( HMA == "" ), 'Adjust', 'Not adjust' ) ) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    
    
} else {
    df_log_PandT_longFormat_simplified_StrataLabels <-
        df_log_PandT_longFormat_simplified_StrataLabels %>%
        # Create new columns that contain the test events and the prescription events.
        dplyr::select( person_id, start_dttm, event_value, idx_test_interval ) %>%
        dplyr::filter( !stringr::str_detect( event_value, pattern = "Test") ) %>%
        dplyr::filter( !stringr::str_detect( event_value, pattern = "Unobserved") ) %>%
        dplyr::mutate( meds_only = as.character( event_value ) ) %>%
        dplyr::select( -event_value ) %>%
        # Remove duplicated prescriptions occuring on the same day.
        dplyr::distinct( ) %>%
        # Create a list-column of the prescriptions given in the preceding months of a given prescription, `HMA_adjust_lookBack_window`.
        tidyr::nest( .by = person_id ) %>%
        dplyr::mutate(
            meds_set_window = 
                lapply(
                    data
                    ,function(x) slider::slide_index(.x = x$meds_only, .i = x$start_dttm, .f = ~.x, .before = HMA_adjust_lookBack_window ) %>% array()
                )
            ) %>%
        tidyr::unnest( cols = c( data, meds_set_window ) ) %>%
        # Remove the current prescription from `meds_set_window`.
        dplyr::mutate(
            meds_set_window = purrr::map2(
                meds_only
                ,meds_set_window
                ,function( med, med_set) {
                    idx <- which( med_set == med )
                    if ( length( idx ) > 0) {
                        med_set[ -idx[ 1 ] ]
                    } else {
                        med_set
                    }
                }
            )
        ) %>%
        # Join this data.frame to one that contains a list-column of the prescriptions given on a particular day.
        dplyr::left_join(
            x = dplyr::select( ., - meds_only )
            # Create a list-column of the prescriptions given on a particular day.
            ,y = stats::aggregate(
                        meds_only ~ person_id + start_dttm
                        ,data = .
                        ,FUN = list
                        ,na.action = NULL
                    ) %>% dplyr::arrange( person_id, start_dttm )
            ,by = join_by( person_id, start_dttm )
            ) %>%
        # Find the difference between the history of meds and the meds on the day.
        dplyr::mutate( new_meds = purrr::map2( meds_only, meds_set_window,  ~!all( .x %in% .y ) ) ) %>%
        as.data.frame() %>%
        dplyr::select( -c( meds_only, meds_set_window ) ) %>%
        # Label any additional meds as ' Adjust'.
        dplyr::rowwise() %>%
        dplyr::mutate( HMA = dplyr::if_else( new_meds, 'Adjust', 'Not adjust' ) ) %>%
        dplyr::ungroup() %>%
        as.data.frame() %>%
        dplyr::select( -new_meds ) %>%
        # Set the entire inter-test intervals to an 'Adjust' or 'Not adjust' value based
        # whether an adjust event occured during that window.
        dplyr::group_by( person_id, idx_test_interval ) %>%
        dplyr::mutate( HMA = dplyr::if_else( 'Adjust' %in% HMA, 'Adjust', 'Not adjust') ) %>%
        dplyr::ungroup() %>%
        as.data.frame() %>%
        # Rejoin to the original dataframe.
        dplyr::select( - start_dttm ) %>%
        dplyr::distinct( ) %>%
        dplyr::right_join(
                df_log_PandT_longFormat_simplified_StrataLabels
                ,by = join_by( person_id, idx_test_interval )
            ,relationship = "one-to-many"
            ) %>%
        # Impute the value for `idx_test_interval` for the event_value == `Unobserved` event.
        dplyr::group_by( person_id ) %>%
        dplyr::mutate(
            idx_test_interval =
                dplyr::if_else(
                    ( event_value == 'Unobserved' ) & ( !all(is.na(idx_test_interval)) )
                    ,max( idx_test_interval, na.rm = TRUE ) + 1
                    ,idx_test_interval
                )
        ) %>%
        dplyr::ungroup() %>%
        # Remove any record with `idx_test_interval` == 0. This exlcudes all prescriptions before
        # the first non-diagnostic test recorded. The syntax below will also exclude the diagnosis event.
        dplyr::filter( idx_test_interval != 0 ) %>%
        suppressWarnings()
}

# The final step is to combine the two components of the stratification's definition into a single variable called HMA.
# Take note that the 'Adjust' label overrules a 'Monitor' label.
df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::mutate(
        HMA =
            dplyr::case_when(
                event_value == "Unobserved" ~ "Unobserved"
                
                ,HMA == "Adjust" ~ "Adjust"
                ,inter_test_duration_discr == "As expected" ~ "Hold"
                ,inter_test_duration_discr == "Shorter than expected" ~ "Monitor"

                ,TRUE ~ NA_character_
            )
    ) %>%
    # Set as a factor datatype and order the factors.
    dplyr::mutate( HMA = factor( HMA, levels = df_HMA_factor %>% dplyr::select( HMA_fct_order ) %>% dplyr::pull() ) ) %>%
    dplyr::arrange( person_id, start_dttm )


#######################################
# HMA and Test Status stratification. #
#######################################
#
# First, I add a column indicating the most-recent test value.
test_status_rollover <- vector( mode = 'character', length = nrow( df_log_PandT_longFormat_simplified_StrataLabels ) )
for (i in 1:nrow( df_log_PandT_longFormat_simplified_StrataLabels ) )
    {
    # Focus on the value in the row of interest.
    val_i <- df_log_PandT_longFormat_simplified_StrataLabels$event_value[ i ]

    # Extract the test status or copy from previous.
    if ( stringr::str_detect( val_i, pattern = "(Test)" ) )
            {
            test_status_rollover[i] <- stringr::str_extract( val_i, '\\b\\w+$')
        } else {
            test_status_rollover[i] <-
                ifelse(
                    i == 1
                    ,stringr::str_extract( val_i, '\\b\\w+$')
                    ,ifelse(
                        ( df_log_PandT_longFormat_simplified_StrataLabels$person_id[ i ] !=
                             df_log_PandT_longFormat_simplified_StrataLabels$person_id[ i-1 ] )
                        ,stringr::str_extract( val_i, '\\b\\w+$')
                        ,test_status_rollover[ i-1 ] )
                )

        } # End IF
    } # End FOR

# Add column of test statuses.
df_log_PandT_longFormat_simplified_StrataLabels <- 
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    tibble::add_column( as.data.frame( test_status_rollover )
        )

# Next, I create the variable indicating the combination of HMA category and most-recent test status value.
df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    # The ordering of the my dplyr::case_when() arguments might seem odd at first.
    # This is because the dplyr::case_when() function applies its arguments
    # sequentially. I check the complicated scenarios for each test status before
    # checking the scenario where a test status is given but there is no value
    # for `n_meds_per_test_interval`.
    dplyr::mutate(
        HMAandTestStatus = dplyr::case_when(
            event_value == "Unobserved" ~ "Unobserved"

            ,( HMA == "Hold" ) & ( test_status_rollover == "Red" ) ~ "Hold Red"
            ,( HMA == "Hold" ) & ( test_status_rollover == "Amber" )  ~ "Hold Amber"
            ,( HMA == "Hold" ) & ( test_status_rollover == "Yellow" )  ~ "Hold Yellow"
            ,( HMA == "Hold" ) & ( test_status_rollover == "Green" )  ~ "Hold Green"

            ,( HMA == "Monitor" ) & ( test_status_rollover == "Red" ) ~ "Monitor Red"
            ,( HMA == "Monitor" ) & ( test_status_rollover == "Amber" )  ~ "Monitor Amber"
            ,( HMA == "Monitor" ) & ( test_status_rollover == "Yellow" )  ~ "Monitor Yellow"
            ,( HMA == "Monitor" ) & ( test_status_rollover == "Green" )  ~ "Monitor Green"

            ,( HMA == "Adjust" ) & ( test_status_rollover == "Red" ) ~ "Adjust Red"
            ,( HMA == "Adjust" ) & ( test_status_rollover == "Amber" )  ~ "Adjust Amber"
            ,( HMA == "Adjust" ) & ( test_status_rollover == "Yellow" )  ~ "Adjust Yellow"
            ,( HMA == "Adjust" ) & ( test_status_rollover == "Green" )  ~ "Adjust Green"

            ,TRUE ~ NA
        )
    ) %>%
    dplyr::mutate(
        HMAandTestStatus = 
            factor(
                HMAandTestStatus
                ,levels = df_HMAandTestStatus_factor %>% dplyr::select( HMAandTestStatus_fct_order ) %>% dplyr::pull()
            )
    )

# Finally, I create a reference table indicating the HMA-and-Test Status stratification                
HMAandTestStatus_display_table <-
    data.frame(
        Value = df_HMAandTestStatus_factor %>% dplyr::pull(HMAandTestStatus_fct_order)
        ,`HMA.component` = c( rep( c( 'Hold', 'Monitor' ,'Adjust' ), each = 4 ), 'Unobserved' )
        ,`Test.component` = c( rep( c( 'Test Status = Red', 'Test Status = Amber'
                                      , 'Test Status = Yellow', 'Test Status = Green')
                                   , times = 3 )
                              , 'Unobserved' )
    )
                      

###########################
# T-and-I stratification. #
###########################
#
# The two components of this stratification are the test status, T, and the degree of intervention in the
# previous inter-test interval, I. The test status is already encoded in the test_status_rollover variable.
# The degree of intervention will, in this diabetes case study, require me to count the unique medication
# names in the event_value column that exists between testing events.
#
# Firstly, I add a column that indicates the count of unique medications prescribed in inter-test intervals.
df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::left_join(
        df_log_PandT_longFormat_simplified_StrataLabels %>%
        dplyr::filter( !stringr::str_detect( event_value, pattern = "(Test)" ) ) %>%
        dplyr::group_by( person_id, idx_test_interval ) %>%
        dplyr::summarise( n_meds_per_test_interval = n_distinct( event_value ) ) %>%
        dplyr::ungroup() %>%
        dplyr::filter( idx_test_interval != 0 )
        ,by = join_by( person_id, idx_test_interval )
    ) %>%
    dplyr::mutate(
        n_meds_per_test_interval =
            dplyr::if_else(
                stringr::str_detect( event_value, pattern = "(Test)" )
                , NA_integer_
                , n_meds_per_test_interval
            )
    ) %>%
    # Fill the test-status row with the values for `n_meds_per_test_interval` that
    # was calculated for the respective intervention-event row.
    dplyr::group_by( person_id, idx_test_interval ) %>%
    tidyr::fill( n_meds_per_test_interval, .direction = "up" ) %>%
    dplyr::ungroup()
                                  
df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    # The ordering of the my dplyr::case_when() arguments might seem odd at first.
    # This is because the dplyr::case_when() function applies its arguments
    # sequentially. I check the complicated scenarios for each test status before
    # checking the scenario where a test status is given but there is no value
    # for `n_meds_per_test_interval`.
    dplyr::mutate(
        TandI = dplyr::case_when(
            event_value == "Unobserved" ~ "Unobserved"

            ,( test_status_rollover == "Red" ) & ( n_meds_per_test_interval == 1 ) ~ "Red One Rx"
            ,( test_status_rollover == "Red" ) & ( n_meds_per_test_interval == 2 ) ~ "Red Two Rx"
            ,( test_status_rollover == "Red" ) & ( n_meds_per_test_interval > 2 ) ~ "Red More Rx"
            ,( test_status_rollover == "Red" ) ~ "Red Zero Rx"

            ,( test_status_rollover == "Amber" ) & ( n_meds_per_test_interval == 1 ) ~ "Amber One Rx"
            ,( test_status_rollover == "Amber" ) & ( n_meds_per_test_interval == 2 ) ~ "Amber Two Rx"
            ,( test_status_rollover == "Amber" ) & ( n_meds_per_test_interval > 2 ) ~ "Amber More Rx"
            ,( test_status_rollover == "Amber" ) ~ "Amber Zero Rx"

            ,( test_status_rollover == "Yellow" ) & ( n_meds_per_test_interval == 1 ) ~ "Yellow One Rx"
            ,( test_status_rollover == "Yellow" ) & ( n_meds_per_test_interval == 2 ) ~ "Yellow Two Rx"
            ,( test_status_rollover == "Yellow" ) & ( n_meds_per_test_interval > 2 ) ~ "Yellow More Rx"
            ,( test_status_rollover == "Yellow" ) ~ "Yellow Zero Rx"

            ,( test_status_rollover == "Green" ) & ( n_meds_per_test_interval == 1 ) ~ "Green One Rx"
            ,( test_status_rollover == "Green" ) & ( n_meds_per_test_interval == 2 ) ~ "Green Two Rx"
            ,( test_status_rollover == "Green" ) & ( n_meds_per_test_interval > 2 ) ~ "Green More Rx"
            ,( test_status_rollover == "Green" ) ~ "Green Zero Rx"

            ,TRUE ~ NA
        )
    ) %>%
    # Set as a factor datatype and order the factors.
    dplyr::mutate(
        TandI =
            factor(
                TandI
                ,levels = df_TandI_factor %>% dplyr::select( TandI_fct_order ) %>% dplyr::pull()
            )
    )


# Finally, I create the variable indicating the Tests-and-Interventions stratification by combining the
# test status with the variable indicating the count of unique medications prescribed in inter-test
# intervals (i.e. with n_meds_per_test_interval).
TandI_display_table <-
    data.frame(
        Value = df_TandI_factor %>% dplyr::pull(TandI_fct_order)
        ,`Test.component` = c( rep( c( 'Test Status = Red', 'Test Status = Amber'
                                      , 'Test Status = Yellow', 'Test Status = Green')
                                   , each = 4 )
                              , 'Unobserved' )
        ,`Intervention.component` = c( rep( c( 'Zero', 'One', 'Two', '> Two' ), times = 4 ), 'Unobserved' )
    )
                      

##################################                      
# Multimorbidity stratification. #
##################################
#
# The strata are defined by whether or not there is a record for at least one of the multimorbidity diagnositic
# codes. Patients with a record for at least one of the multimorbidity diagnositic codes have already been
# identified in qry_log_multimorb_longFormat. This BigQuery query result needs 'collecting' and joining to
# df_log_PandT_longFormat_simplified_StrataLabels. Then, I need to create the indicator variable with values of
# Multimorbid = FALSE before the date_multimorb and Multimorbid = TRUE on or after date_multimorb.
#
# I also include a stratification called TandMultiMorb that, like TandI, combines the test status values with
# the multimorbidity values. This stratification has eight levels with an extra for errors.

df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    # Create the MultiMorb stratification.
    dplyr::left_join( qry_log_multimorb_longFormat %>% dplyr::collect(), by = join_by( person_id ) ) %>%
    dplyr::mutate( MultiMorb = dplyr::if_else( date_multimorb <= start_dttm, TRUE, FALSE ) ) %>%
    dplyr::mutate( MultiMorb = dplyr::if_else( is.na( MultiMorb ), FALSE, MultiMorb ) ) %>%
    dplyr::select( - date_multimorb ) %>%
    # Create the TandMultiMorb stratification.
    dplyr::mutate(
        TandMultiMorb = dplyr::case_when(
            event_value == "Unobserved" ~ "Unobserved"

            ,( test_status_rollover == "Red" ) & ( MultiMorb == TRUE ) ~ "Red Multimorbid"
            ,( test_status_rollover == "Amber" ) & ( MultiMorb == TRUE ) ~ "Amber Multimorbid"
            ,( test_status_rollover == "Yellow" ) & ( MultiMorb == TRUE ) ~ "Yellow Multimorbid"
            ,( test_status_rollover == "Green" ) & ( MultiMorb == TRUE ) ~ "Green Multimorbid"

            ,( test_status_rollover == "Red" ) & ( MultiMorb == FALSE ) ~ "Red Not multimorbid"
            ,( test_status_rollover == "Amber" ) & ( MultiMorb == FALSE ) ~ "Amber Not multimorbid"
            ,( test_status_rollover == "Yellow" ) & ( MultiMorb == FALSE ) ~ "Yellow Not multimorbid"
            ,( test_status_rollover == "Green" ) & ( MultiMorb == FALSE ) ~ "Green Not multimorbid"

            ,TRUE ~ NA
        )
    ) %>%
    # Set as a factor datatype and order the factors.
    dplyr::mutate(
        TandMultiMorb =
              factor(
                  TandMultiMorb
                  ,levels = df_TandMultiMorb_factor %>% dplyr::select( TandMultiMorb_fct_order ) %>% dplyr::pull()
              )
    )


# Finally, I create the variable indicating the Tests-and-Multimorbidity stratification by combining the test
# status with the variable indicating multimorbidity.
TandMultiMorb_display_table <-
    data.frame(
        Value = df_TandMultiMorb_factor %>% dplyr::pull(TandMultiMorb_fct_order)
        ,`Test.component` = c( rep( c( 'Test Status = Red', 'Test Status = Amber'
                                      , 'Test Status = Yellow', 'Test Status = Green'), times = 2 ), 'Unobserved' )
        ,`Multimorbidity.component` = c( rep( c( 'Multimorbid', 'Not multimorbid' ), each = 4 ), 'Unobserved' )
    )