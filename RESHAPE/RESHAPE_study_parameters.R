# RESHAPE_study_parameters.R
#
# The purpose of this script is to specify the parameters that will be used by
# the other scripts in the project.

# Study dates
# ## The date before which a patient must have had their diagnosis.
date_diagnosis_threshold <- lubridate::ymd('2000-01-01')

# ## The date after which test and intervention records will be studied.
followup_delay_in_years <- 0
date_followup_start <- date_diagnosis_threshold + lubridate::years( followup_delay_in_years )

# ## The date before which test and intervention records will be studied.
followup_duration_in_years <- 10
date_followup_end <- date_followup_start + lubridate::years( followup_duration_in_years )


# Set the duration of the window back in time to review prescriptions when identifying
# the 'Adjust' state.
# - Use `HMA_adjust_lookBack_window` if you want to look back in chronological time.
# - Use `HMA_adjust_lookBack_count` if you want to looks back through prescriptions.
HMA_adjust_lookBack_window <- NULL #lubridate::weeks( 16 )
HMA_adjust_lookBack_count <- 5 # NULL
if( !is.null( HMA_adjust_lookBack_window ) & !is.null( HMA_adjust_lookBack_count ) )  
{
  stop( "Only one of `HMA_adjust_lookBack_window` or `HMA_adjust_lookBack_count` can be specified as non NULL." )
}


# Set upper and lower thresholds for acceptable values of the test.
test_value_cutoff_lower <- 20
test_value_cutoff_upper <- 200


# Threshold for the expected interval between subsequent tests, in months
val_testing_interval_LB <- 2
val_testing_interval_UB <- 5


# Set values for meaningful changes in the values of the test.
val_meaningful_test_improvement <- -10
val_meaningful_test_disimprovement <- 10


# Set window within which to search for repeated (but not repeat) prescriptions.
window_repeated_prescription_months <- 3


# Set number of tests, treatments, or iterations after diagnosis that should be tracked.
n_iterations <- followup_duration_in_years*2


# Set the window within which mutimorbidity diagnoses and the index diagnosis must fit in, in months.
multimorb_inclusion_window_months <- 60


# Set the window outwith which at least two mutimorbidity diagnoses must be of each other, in months.
multimorb_gap_window_months <- 1