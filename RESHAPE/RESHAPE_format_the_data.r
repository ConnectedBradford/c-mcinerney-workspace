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
                      
# Create a variable that indicates patient-specific test interval.
df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified %>%
    dplyr::group_by( person_id ) %>%
    dplyr::arrange( start_dttm ) %>%
    dplyr::mutate(
        idx_test_interval = 0 + cumsum( stringr::str_detect( event_value, pattern = "Test" ) )
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
    df_log_PandT_longFormat_simplified_StrataLabels %>%
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
    dplyr::arrange( person_id, start_dttm ) %>%
    dplyr::mutate(
        end_dttm = dplyr::if_else( is.na( end_dttm ), lead( start_dttm ), end_dttm )
        ,event_name = dplyr::if_else( is.na( event_value ), "diagnosis", event_name )
    ) %>%
    dplyr::ungroup() 


#########################
# H.M.A. stratification #
#########################
#
# The two components of this stratification are the testing interval and the change in prescription. The testing
# interval requries me to create a variable indicating that the inter-test duration was between
# `val_testing_interval_LB` and `val_testing_interval_UB`. The change-in-prescription criterion refers to 
# whether a new medication-dose set was prescribed, where "new" means not seen in the previous three months.
#
# The first thing is to calculate the inter-test duration. In this iteration, I will have to ignore the interval
# between the diagnostic test and the subsequent test because I am deliberately looking at records ten years
# after the diagnostic test.# ****************************************************
                            # **************************************************** This isn't always true. For the paper I might be doing straight after diagnosis.
                            # ****************************************************
                      
df_inter_test_duration <-
    qry_log_test_longFormat %>%
    dplyr::collect() %>%
    # Include the date of diagnosis.
    dplyr::group_by( person_id ) %>%
    dplyr::group_modify( ~ tibble::add_row( .x ), .by = person_id ) %>% 
    dplyr::left_join( qry_records_with_T2DM_diagnoses %>% collect(), by = join_by( person_id ) ) %>%
    dplyr::mutate( start_dttm = if_else( is.na( start_dttm ), date_diagnosis , start_dttm ) ) %>%
    dplyr::arrange( person_id, start_dttm ) %>%
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
    dplyr::distinct() %>% dplyr::arrange(person_id, start_dttm)

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
    dplyr::arrange( start_dttm ) %>%
    tidyr::fill( inter_test_duration_cont ) %>%
    tidyr::fill( inter_test_duration_discr ) %>%
    dplyr::mutate( inter_test_duration_cont = dplyr::if_else( inter_test_duration_discr == "Not applicable", NA_real_, inter_test_duration_cont ) ) %>%
    dplyr::ungroup()
  
                      
# Next, I need to make the variable that indicates whether the prescription has changed *since*??????
df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    # Create new columns that contain the test events and the prescription events.
    dplyr::select( person_id, start_dttm, event_value, idx_test_interval ) %>%
    dplyr::mutate(
        tempCol_tests = dplyr::if_else( stringr::str_detect( event_value, pattern = "Test"), event_value, NA_character_ )
    ) %>%
    dplyr::mutate(
        tempCol_prescrps = dplyr::if_else( !stringr::str_detect( event_value, pattern = "Test"), event_value, NA_character_ )
    ) %>%
    # Fill the test column until a new value is given.
    dplyr::arrange( person_id, start_dttm, tempCol_prescrps ) %>%
    tidyr::fill( tempCol_tests ) %>%
    # Tidy away unecessary columns.
    #dplyr::select( - c( event_value ) ) %>%
    # Remove repeated prescriptions within an inter-test period.
    dplyr::distinct() %>%
    # Pivot wider so that each row represents an inter-test period. Order within the grouping doesn't matter because
    # I'm just collecting the set of prescriptions rather than the sequence of them.
    dplyr::group_by( person_id, tempCol_tests ) %>%
    dplyr::mutate( rn = row_number() ) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(
        names_from = rn
        ,values_from = tempCol_prescrps
        ) %>%
    # List-concatonate the repeated prescriptions into one string.
    tidyr::unite(
        col = "meds_set"
        ,dplyr::select(
            .
            ,- c('person_id', 'tempCol_tests', 'start_dttm', 'idx_test_interval')
        ) %>% colnames()
        ,sep = " , "
        ,na.rm = TRUE
    ) %>%
    # Create a list of the `meds_set`s that have occurred in the previous 
    # `HMA_adjust_lookBack_window_in_wks`, inclusive of the given row's `med_set`.
    dplyr::mutate(
        meds_set = strsplit( meds_set, split = " , ", fixed = TRUE )
    ) %>% 
    dplyr::group_by( person_id ) %>%
    dplyr::mutate(
        meds_set_window = slider::slide_index(.x = meds_set, .i = start_dttm, .f = ~.x, .before = HMA_adjust_lookBack_window )
    ) %>%
    dplyr::ungroup() %>%
    # Remove the current `meds_set` from the `meds_set_window`.
    dplyr::rowwise() %>%
    dplyr::mutate(
        meds_set_window = 
            meds_set_window %>%
            head( -1 ) %>%
            unlist() %>%
            unique() %>%
            list()
    ) %>%
    dplyr::ungroup() %>%
    # Check if the current `meds_set` is not an element of the `meds_set_window`. Label as
    # 'HMA_adjust` == TRUE if so.
    dplyr::group_by( person_id ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        HMA_adjust =
            ifelse(
                ( length( meds_set ) > 0 )
                ,!meds_set %in% meds_set_window
                ,NA
                )
        ,HMA_adjust = dplyr::if_else( all( meds_set == "Unobserved" ), NA, HMA_adjust)
    ) %>%
    dplyr::ungroup() %>%
    # Set entire inter-test intervals to an 'Adjust' or 'Not adjust' value based
    # whether an adjust event occured during that window.
    dplyr::group_by( person_id, idx_test_interval ) %>%
    dplyr::mutate( adjust_interval = dplyr::if_else( 'TRUE' %in% HMA_adjust, 'Adjust', 'Not adjust') ) %>%
    dplyr::ungroup() %>%
    # Join the `new_meds` column to the original dataframe.
    dplyr::distinct( person_id, idx_test_interval, adjust_interval ) %>%
    dplyr::right_join(
        df_log_PandT_longFormat_simplified_StrataLabels
        ,by = join_by( person_id, idx_test_interval )
    ) %>%
    # Impute the value for `idx_test_interval` for the event_value == `Unobserved` event.
    dplyr::group_by( person_id ) %>%
    dplyr::mutate(
        idx_test_interval = dplyr::if_else( ( event_value == 'Unobserved' ) & ( !all(is.na(idx_test_interval)) ), max( idx_test_interval, na.rm = TRUE ) + 1, idx_test_interval )
    ) %>%
    dplyr::ungroup() %>%
    # Remove any record with `idx_test_interval` == 0. This exlcudes all prescriptions before
    # the first non-diagnostic test recorded. The syntax below will also exclude the diagnosis event.
    dplyr::filter( idx_test_interval != 0 ) %>%
    dplyr::arrange( person_id, start_dttm ) %>%
    suppressWarnings()
                     
# The final step is to combine the two components of the stratification's definition into a single variable called HMA.
# Take note that the 'Adjust' label overrules a 'Monitor' label.
df_log_PandT_longFormat_simplified_StrataLabels <-
    df_log_PandT_longFormat_simplified_StrataLabels %>%
    dplyr::mutate(
        HMA =
            dplyr::case_when(
                event_value == "Unobserved" ~ "Unobserved"
                
                ,adjust_interval == "Adjust" ~ "Adjust"
                , inter_test_duration_discr == "As expected" ~ "Hold"
                , inter_test_duration_discr == "Shorter than expected" ~ "Monitor"

                ,TRUE ~ NA_character_
            )
    ) %>%
    dplyr::select( -adjust_interval ) %>%
    # Set as a factor datatype and order the factors.
    dplyr::mutate( HMA = factor( HMA, levels = df_HMA_factor %>% dplyr::select( HMA_fct_order ) %>% dplyr::pull() ) )


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