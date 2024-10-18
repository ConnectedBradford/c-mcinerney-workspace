# The purpose of this script is to create the sample cohort that satisfies the following characteristics:
#
#
#
#####################
# Install pacakges. #
#####################
if( !"pacman" %in% installed.packages() ){ install.packages( "pacman" ) }
pacman::p_load(
    bigrquery
    ,paletteer
    ,tidyverse
    ,slider
    )

##############################
## Set and load requisites. ##
##############################

# Setup connection to GCP.
project_id = "yhcr-prd-bradfor-bia-core"
con <- DBI::dbConnect( drv = bigquery(), project = project_id ) %>% suppressWarnings()

# Define R tibbles from GCP tables.
r_tbl_srcode <- dplyr::tbl( con, sql( "SELECT * FROM CB_FDM_PrimaryCare.tbl_srcode" ) )
r_tbl_srpatient <- dplyr::tbl( con, sql( "SELECT * FROM CB_FDM_PrimaryCare.tbl_srpatient" ) )
r_tbl_BNF_DMD_SNOMED_lkp <- dplyr::tbl( con, sql( "SELECT * FROM ref_Lookups.lkp_BNF_DMD_SNOMED_lkp" ) )
r_tbl_srprimarycaremedication <- dplyr::tbl( con, sql( "SELECT * FROM CB_FDM_PrimaryCare.tbl_srprimarycaremedication" ) )
r_tbl_srappointment <- dplyr::tbl( con, sql( "SELECT * FROM CB_FDM_PrimaryCare.tbl_srappointment" ) )

# Clinical code lists (BNF, SNOMED-CT, etc).
# ## Medication codes.
codes_SNOMED_diagnoses_of_interest <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-dmtype2_cod-20200812.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_test_of_interest <-
    readr::read_csv(file = paste0( 'codelists/', 'opensafely-glycated-haemoglobin-hba1c-tests-3e5b1269.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_BNF_metformin <-
    readr::read_csv(file = paste0( 'codelists/', 'ciaranmci-metformin-bnf-0601022b0-and-child-bnf-codes-only-43e7d87e.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
names_metformin_meds <-
    r_tbl_BNF_DMD_SNOMED_lkp %>%
    dplyr::filter( BNF_Code %in% codes_BNF_metformin ) %>%
    dplyr::select( DMplusD_ProductDescription )
BNF_meds_of_interest <-
    readr::read_csv(file = paste0( 'codelists/', 'ciaranmci-bnf-section-61-drugs-for-diabetes-207573b7.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )
codes_BNF_meds_of_interest <- BNF_meds_of_interest$code
meds_names_simplified <-
    readr::read_csv( file = paste0( 'codelists/', 'meds_names_simplified.csv' ),
                    col_types = cols( term = col_character(), term_simplified = col_character() ) )
# ## ## Amalgamate medicine preparations.
names_meds_of_interest <-
    r_tbl_BNF_DMD_SNOMED_lkp %>%
    dplyr::filter( BNF_Code %in% codes_BNF_meds_of_interest ) %>%
    dplyr::select( BNF_Code, DMplusD_ProductDescription )
shortNames_meds_of_interest <-
    names_meds_of_interest %>%
    # Create a shortened version of the BNF code to use as a linking key.
    dplyr::mutate(
        BNF_Code_short =
            dplyr::if_else(
                stringr::str_sub( BNF_Code, end = 6 ) == "060101"
                ,"060101"
                ,stringr::str_sub( BNF_Code, end = 9 )
            )
        ) %>%
    # Join with tables of medications of interest.
    dplyr::collect() %>%
    dplyr::left_join(
         BNF_meds_of_interest %>% filter( nchar( code ) == 9 )
        ,by = join_by( BNF_Code_short == code )
    ) %>%
    dplyr::left_join(
        meds_names_simplified
        ,by = join_by( DMplusD_ProductDescription == term )
    ) %>%
    # Fill in the NAs in `term_simplified` for insulin and metformin hydrochloride.
    dplyr::mutate(
        term_simplified = dplyr::case_when(
            !is.na( term_simplified ) ~ term_simplified
            
            ,BNF_Code_short == "060101"  ~ "Insulin"
            ,term == "Metformin hydrochloride" ~ "Metformin hydrochloride"
            
            ,TRUE ~ NA
        )
    ) %>%
    # Remove glucose preparations.
    dplyr::filter( term_simplified != "Glucose" ) %>%
    # Tidy up.
    dplyr::rename(
        drug_name = DMplusD_ProductDescription
        ,drug_name_short = term_simplified
    ) %>%
    dplyr::select( c( drug_name, drug_name_short ) ) %>%
    dplyr::distinct()

# ## Multimorbidity codes.
codes_SNOMED_depression <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-depr_cod-20210127.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_anxiety <-
    readr::read_csv(file = paste0( 'codelists/', 'opensafely-anxiety-disorders-6aef605a.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_substanceAbuse <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-illsub_cod-20210127.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_dementia <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-dem_cod-20210127.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_SMI <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-mh_cod-20210127.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_arthritis <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-c19rarth_cod-20210127.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_lungCancer <-
    readr::read_csv(file = paste0( 'codelists/', 'opensafely-lung-cancer-snomed-2020-04-15.csv' ),
                    col_types = cols( id = col_character(), name = col_character() ) )$id
codes_SNOMED_otherCancer <-
    readr::read_csv(file = paste0( 'codelists/', 'opensafely-cancer-excluding-lung-and-haematological-snomed-2020-04-15.csv' ),
                    col_types = cols( id = col_character(), name = col_character() ) )$id
codes_SNOMED_haemCancer <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-c19haemcan_cod-20210127.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_hypertension <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-hyp_cod-20210127.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_hyperlipid <-
    readr::read_csv(file = paste0( 'codelists/', 'ciaranmci-hyperlipidaemia-68d4d269.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_cardiacArrhythmia <-
    readr::read_csv(file = paste0( 'codelists/', 'ciaranmci-cardiac-arrhythmias-4fc55f34.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_asthma <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-ast_cod-20210127.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_CAD <-
    readr::read_csv(file = paste0( 'codelists/', 'ciaranmci-coronary-artery-disease-36b5ebff.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_COPD <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-copd_cod-20210127.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_osteo <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-osteo_cod-20200812.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_CKD1and2 <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-ckd1and2atrisk1_cod-20200812.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_CKD3to5 <-
    readr::read_csv(file = paste0( 'codelists/', 'nhsd-primary-care-domain-refsets-ckdatrisk2_cod-20210127.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_CHF <-
    readr::read_csv(file = paste0( 'codelists/', 'ciaranmci-congestive-heart-failure-1da678ca.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_hepatitis <-
    readr::read_csv(file = paste0( 'codelists/', 'ciaranmci-hepatitis-04970595.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code
codes_SNOMED_all_multimorbidity_diagnoses <-
    c(
        codes_SNOMED_anxiety
        ,codes_SNOMED_arthritis
        ,codes_SNOMED_asthma
        ,codes_SNOMED_CAD
        ,codes_SNOMED_cardiacArrhythmia
        ,codes_SNOMED_CHF
        ,codes_SNOMED_CKD1and2
        ,codes_SNOMED_CKD3to5
        ,codes_SNOMED_COPD
        ,codes_SNOMED_depression
        ,codes_SNOMED_dementia
        ,codes_SNOMED_haemCancer
        ,codes_SNOMED_hepatitis
        ,codes_SNOMED_hyperlipid
        ,codes_SNOMED_hypertension
        ,codes_SNOMED_lungCancer
        ,codes_SNOMED_osteo
        ,codes_SNOMED_otherCancer
        ,codes_SNOMED_SMI
        ,codes_SNOMED_substanceAbuse
    )
# ## Did-not-attend.
codes_SNOMED_didNotAttend <-
    readr::read_csv(file = paste0( 'codelists/', 'ciaranmci-did-not-attend-098119da.csv' ),
                    col_types = cols( code = col_character(), term = col_character() ) )$code

# Load the CSV file containing that matches all medications in `codes_SNOMED_diagnoses_of_interest`
# with an expected duration of prescription.
##meds_expected_prescription_duration <-
 ##   readr::read_csv(file = 'meds_expected_prescription_duration.csv',
  ##                  col_types = cols( term = col_character(), duration = col_numeric() ) )

# Study dates
# ## The date before which a patient must have had their diagnosis.
date_diagnosis_threshold <- lubridate::ymd('2000-01-01') # NULL
# ## The date after which test and intervention records will be studied.
followup_delay_in_years <- 0
date_followup_start <- date_diagnosis_threshold + lubridate::years( followup_delay_in_years )
# ## The date before which test and intervention records will be studied.
followup_duration_in_years <- 10
date_followup_end <- date_followup_start + lubridate::years( followup_duration_in_years )

# Set the duration of the window back in time to review prescriptions when identifying
# the HMA status.
HMA_adjust_lookBack_window <- lubridate::weeks( 16 )

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
n_iterations <- followup_duration_in_years*4

# Set the window within which mutimorbidity diagnoses and the index diagnosis must fit in, in months.
multimorb_inclusion_window_months <- 60

# Set the window outwith which at least two mutimorbidity diagnoses must be of each other, in months.
multimorb_gap_window_months <- 1

# Set ordering of factors labels for the event variable.
temp_meds_names <-
    shortNames_meds_of_interest %>%
    dplyr::distinct( drug_name_short ) %>%
    dplyr::arrange( drug_name_short ) %>%
    dplyr::pull( drug_name_short )

df_event_factor <-
    data.frame(
        event_fct_order =
            factor(
                c( "Test Status = Red", "Test Status = Amber"
                  ,"Test Status = Yellow", "Test Status = Green"
                  ,temp_meds_names
                  ,"Unobserved"
                 )
                ,levels = c( "Test Status = Red", "Test Status = Amber"
                            ,"Test Status = Yellow", "Test Status = Green"
                            ,temp_meds_names
                            ,"Unobserved" 
                           )
            )
        ,event_colours_order = c( "firebrick1", "darkorange", "gold", "limegreen"
                                 ,paletteer_c( "grDevices::Plasma" ,temp_meds_names %>% length() )
                                 ,"grey"
                                 )
    )

# Set ordering of H.M.A. stratification.
df_HMA_factor <-
    data.frame(
        HMA_fct_order =
            factor(
                c("Hold", "Monitor", "Adjust", "Unobserved")
                ,levels = c("Hold", "Monitor","Adjust", "Unobserved")
            )
        ,HMA_colours_order = c( "plum1", "lightcoral", "cornflowerblue", "grey" )
    )


# Set ordering of HMAandTestStatus stratification.
df_HMAandTestStatus_factor <-
    data.frame(
        HMAandTestStatus_fct_order =
            factor(
                c(
                  "Hold Red" ,"Hold Amber" ,"Hold Yellow" ,"Hold Green"
                  ,"Monitor Red","Monitor Amber", "Monitor Yellow", "Monitor Green"
                  ,"Adjust Red", "Adjust Amber", "Adjust Yellow", "Adjust Green"
                  ,"Unobserved"
                 )
                ,levels = 
                    c(
                      "Hold Red" ,"Hold Amber" ,"Hold Yellow" ,"Hold Green"
                      ,"Monitor Red","Monitor Amber", "Monitor Yellow", "Monitor Green"
                      ,"Adjust Red", "Adjust Amber", "Adjust Yellow", "Adjust Green"
                      ,"Unobserved"
                     )
            )
        ,HMAandTestStatus_colours_order =
                c( rep( c("firebrick1", "darkorange", "gold", "limegreen"), times = 3), "grey" )
    )

# Set ordering of TandI stratification.
df_TandI_factor <-
    data.frame(
        TandI_fct_order =
            factor(
                c(
                    'Red Zero Rx', 'Red One Rx', 'Red Two Rx', 'Red More Rx',
                    'Amber Zero Rx', 'Amber One Rx', 'Amber Two Rx', 'Amber More Rx',
                    'Yellow Zero Rx', 'Yellow One Rx', 'Yellow Two Rx', 'Yellow More Rx',
                    'Green Zero Rx', 'Green One Rx', 'Green Two Rx', 'Green More Rx',
                    'Unobserved'
                )
                ,levels =
                c(
                    'Red Zero Rx', 'Red One Rx', 'Red Two Rx', 'Red More Rx',
                    'Amber Zero Rx', 'Amber One Rx', 'Amber Two Rx', 'Amber More Rx',
                    'Yellow Zero Rx', 'Yellow One Rx', 'Yellow Two Rx', 'Yellow More Rx',
                    'Green Zero Rx', 'Green One Rx', 'Green Two Rx', 'Green More Rx',
                    'Unobserved'
                )
            )
        ,TandI_colours_order =
                c(
                    "red" # 1: Test status == Red & Zero prescriptions
                    ,"firebrick3" # 2: Test status == Red & One prescriptions
                    ,"firebrick4" # 3: Test status == Red & Two prescriptions
                    ,"darkred" # 4: Test status == Red & Three prescriptions

                    ,"orange" # 5: Test status == Amber & Zero prescriptions
                    ,"darkorange2" # 6: Test status == Amber & One prescriptions
                    ,"darkorange3" # 7: Test status == Amber & Two prescriptions
                    ,"darkorange4" # 8: Test status == Amber & Three prescriptions

                    ,"yellow" # 9: Test status == Yellow & Zero prescriptions
                    ,"gold2" # 10: Test status == Yellow & One prescriptions
                    ,"gold3" # 11: Test status == Yellow & Two prescriptions
                    ,"gold4" # 12: Test status == Yellow & Three prescriptions

                    ,"green" # 13: Test status == Green & Zero prescriptions
                    ,"chartreuse2" # 14: Test status == Green & One prescriptions
                    ,"chartreuse3" # 15: Test status == Green & Two prescriptions
                    ,"chartreuse4" # 16: Test status == Green & Three prescriptions

                    ,"grey" # 0 = Unobserved
                )
    )

# Set ordering of TandMultiMorb stratification.
df_TandMultiMorb_factor <-
    data.frame(
        TandMultiMorb_fct_order =
            factor(
                c(
                    'Red Multimorbid', 'Amber Multimorbid', 'Yellow Multimorbid', 'Green Multimorbid',
                    'Red Not multimorbid', 'Amber Not multimorbid', 'Yellow Not multimorbid', 'Green Not multimorbid',
                    'Unobserved'
                )
                ,levels =
                c(
                    'Red Multimorbid', 'Amber Multimorbid', 'Yellow Multimorbid', 'Green Multimorbid',
                    'Red Not multimorbid', 'Amber Not multimorbid', 'Yellow Not multimorbid', 'Green Not multimorbid',
                    'Unobserved'
                )
            )
        ,TandMultiMorb_colours_order =
                c(
                    "red4" # 1: Test status == Red & Multimorbid
                    ,"orangered3" # 2: Test status == Amber & Multimorbid
                    ,"gold2" # 3: Test status == Yellow & Multimorbid
                    ,"darkgreen" # 4: Test status == Green & Multimorbid
                    
                    ,"red" # 5: Test status == Red & Not multimorbid
                    ,"orange" # 6: Test status == Amber & Not multimorbid
                    ,"yellow" # 7: Test status == Yellow & Not multimorbid
                    ,"green" # 8: Test status == Green & Not multimorbid

                    ,"grey" # 0 = Unobserved
                )
    )


########################
# Define study cohort. #
########################
# This was trivial in previous iterations because I relied on clinicall-coded diagnoses, only. This time I
# will also permit diagnosis to be indicated by abnormal HbA1c (i.e. >48 mmol/mol). The justification is that
# concurrent work by CMI and colleagues shows that clinicall-coded dates of diabetes diagnosis disagree with
# HbA1c values by more than 10 years. The HbA1c values are surperior indicators but the validity of HbA1c
# values decreases as we go further back in the record because they weren't used diagnostically until a while
# after introduction. Therefore, in the early days, we will have to trust the clinically-coded date.
#
# In an email sent 26th July 2024, CB suggested a three-option algorthin for identifying the date of diagnosis
# assuming a `date_diagnosis_threshold` == "2000-01-01":
#
# 1. If the clinically-coded diagnosis is before April 2004 AND there are raised HbA1c values after April 2003
#    (note the difference in years), then use the clinically-coded date.
#    - The actual statement from CB was: if the clinically-coded diagnosis is before April 2004 AND some "recent"
#      raised HbA1c, then use the clinically-coded date, on the assumption that the date is not miscoded. This
#      requries a threshold for "recent" (or does he mean "concurrent"?).
# 2. If the clinically-coded diagnosis is after April 2004 AND there are raised HbA1c values before April 2003
#    (note the difference in years), then use the earliest raised HbA1c date.
# 3. If the clinically-coded diagnosis is after April 2004 AND the first raised HbA1c value is after April 2003
#    (note the difference in years), then use the clinically-coded date.
#
# The years 2003 and 2004 are relative to the assumed baseline of 2000, so reference to them will be calculated
# as lubridate::year( date_diagnosis_threshold + lubridate::years( x ) ), where x = {3, 4}, respectively.
# CB's definition does not specify what happens when both the clinically-coded and serological-informed dates
# are before 2003. In these cases, I will choose the serologically-coded date.
#
# My approach to identifying the date of diagnosis will be to start with the original method of using clinical
# codes, and then only change the date if it satisfies option #2, above.
#
# First, I need to convert historic A1c% values to up-to-date mmol/mol values
# (Formala taken from https://ebmcalc.com/GlycemicAssessment.htm). Then I will define the cohort of relevant records
#
r_tbl_srcode_tests <-
    r_tbl_srcode %>%
    dplyr::select( person_id, dateevent, numericvalue, snomedcode ) %>%
    # Convert A1c% values to mmol/mol values. The rationale for the following syntax is
    # in `RESHAPE_check_HbA1c_SNOMED_codes.ipynb`.
    # ## Apply transformations, where appropriate
    dplyr::mutate( numericvalue = numericvalue %>% as.numeric() ) %>%
    dplyr::filter(
        dplyr::case_when(
            snomedcode %in% codes_SNOMED_test_of_interest ~ numericvalue > 2.15
            ,.default = numericvalue == numericvalue
        )
    ) %>%         
    dplyr::mutate(
        numericvalue = 
            dplyr::if_else(
                ( snomedcode %in% codes_SNOMED_test_of_interest ) & ( numericvalue < 20 )
                ,( numericvalue %>% `-`( 2.15 ) ) %>% `*`( 10.929 )
                ,numericvalue
            )
        ,testType = 
                dplyr::case_when(
                   snomedcode == '1003671000000109'  ~ "Median 58"
                    ,snomedcode =='999791000000106'  ~ "Median 42"
                    ,snomedcode == '1049301000000100' ~ "Median 43"
                    ,TRUE ~ NA_character_
                )
        ) %>%
    dplyr::filter( ( numericvalue >= 20 ) | is.na( numericvalue ) ) %>%
    # ## Apply bandpass filter.
    dplyr::filter(
        dplyr::case_when(
            snomedcode %in% codes_SNOMED_test_of_interest ~ dplyr::between( numericvalue, test_value_cutoff_lower, test_value_cutoff_upper )
            ,.default = numericvalue == numericvalue
        )
    ) %>%
    # Recast `numericvalue` to string.
    dplyr::mutate( numericvalue = numericvalue %>% as.character() )

r_tbl_srcode_diagnoses <-
    r_tbl_srcode %>%
    dplyr::select( person_id, dateevent, numericvalue, snomedcode ) %>%
    dplyr::filter( snomedcode %in% codes_SNOMED_diagnoses_of_interest ) 

r_tbl_srcode_other <-
    r_tbl_srcode %>%
    dplyr::select( person_id, dateevent, numericvalue, snomedcode ) %>%
    dplyr::filter( ( !snomedcode %in% codes_SNOMED_test_of_interest ) & ( !snomedcode %in% codes_SNOMED_diagnoses_of_interest ) )

r_tbl_srcode <-
    r_tbl_srcode_tests %>%
    dplyr::union_all( r_tbl_srcode_diagnoses ) %>%
    dplyr::union_all( r_tbl_srcode_other )

rm( r_tbl_srcode_tests, r_tbl_srcode_diagnoses, r_tbl_srcode_other )

date_coded_diag_threshold <- date_diagnosis_threshold + lubridate::years( 4 )
qry_records_with_T2DM_diagnoses_coded <-
    r_tbl_srcode %>%
    # Identify records with a clinical code for Type 2 Diabetes Mellitus.
    dplyr::filter( snomedcode %in% codes_SNOMED_diagnoses_of_interest ) %>%
    dplyr::group_by( person_id ) %>%
    dplyr::summarise( date_coded_diagnosis = min( dateevent, na.rm = TRUE ) ) %>%
    dplyr::ungroup() %>%
    dplyr::select( person_id, date_coded_diagnosis ) %>%
    # Identify records whose clinically-coded date of diagnosis is after April 2004.
    dplyr::mutate(
        coded_diagnosis_before_April2004 =
            dplyr::if_else(
                ( date_coded_diagnosis < date_coded_diag_threshold )
                ,TRUE
                ,FALSE
            )
    )

date_HbA1c_diag_threshold <- date_diagnosis_threshold + lubridate::years( 3 )  
qry_records_with_T2DM_diagnoses_HbA1c <-
    r_tbl_srcode %>%
    dplyr::filter( snomedcode %in% codes_SNOMED_test_of_interest) %>%
    dplyr::select( person_id, dateevent, numericvalue, testType ) %>%
    # Filter for raised test scores.
    dplyr::filter( as.numeric( numericvalue ) > 48 ) %>%
    # Identify records with raised HbA1c values before April 2003.
    dplyr::mutate(
        raised_HbA1c_before_April2003 =
            dplyr::if_else(
                ( dateevent  < date_HbA1c_diag_threshold )
                ,TRUE
                ,FALSE
            )
    ) %>%
    # Filter for the earliest test that has satisfied the criteria so far.
    dplyr::group_by( person_id ) %>%
    dplyr::filter( dateevent == min( dateevent ) ) %>%
    # Filter for the largest score because, sometimes, a patient has two tests on the same date, one of which seems spuriously low.
    dplyr::filter( as.numeric( numericvalue ) == max( as.numeric( numericvalue ) ) ) %>%
    dplyr::ungroup() %>%
    # Tidy up.
    dplyr::rename( date_HbA1c_diagnosis = dateevent )

qry_records_with_T2DM_diagnoses <-
    # Join the two query results.
    qry_records_with_T2DM_diagnoses_coded %>%
    dplyr::left_join( qry_records_with_T2DM_diagnoses_HbA1c, by = join_by( person_id ) ) %>%
    # Choose the date based on CB's definitions.
    dplyr::mutate(
        date_diagnosis = 
            dplyr::case_when(
                is.na( date_HbA1c_diagnosis ) ~ date_coded_diagnosis
                ,is.na( date_coded_diagnosis ) ~ date_HbA1c_diagnosis
                # Option 1.
                ,( coded_diagnosis_before_April2004 == TRUE ) & ( raised_HbA1c_before_April2003 == FALSE ) ~ date_coded_diagnosis
                # Option 2.
                ,( coded_diagnosis_before_April2004 == FALSE ) & ( raised_HbA1c_before_April2003 == TRUE ) ~ date_HbA1c_diagnosis
                # Option 3.
                ,( coded_diagnosis_before_April2004 == FALSE ) & ( raised_HbA1c_before_April2003 == FALSE ) ~ date_coded_diagnosis
                
                ,( coded_diagnosis_before_April2004 == TRUE ) & ( raised_HbA1c_before_April2003 == TRUE ) ~ date_coded_diagnosis

                ,TRUE ~ NA_character_
            )
    ) %>%
    # Tidy up.
    dplyr::select( person_id, date_diagnosis, date_coded_diagnosis )
# Filter for records who diagnosis date is before the threshold.
if( !is.null( date_diagnosis_threshold ) )
    {
    qry_records_with_T2DM_diagnoses <- dplyr::filter( qry_records_with_T2DM_diagnoses, date_coded_diagnosis <= date_diagnosis_threshold )
} 
qry_records_with_T2DM_diagnoses <- dplyr::select( qry_records_with_T2DM_diagnoses, person_id, date_diagnosis ) %>% dplyr::distinct()


# Retrieve dates of prescriptions in the follow-up period.
adjusted_date_followup_start <-date_followup_start - HMA_adjust_lookBack_window
qry_log_prescription_longFormat <-
    qry_records_with_T2DM_diagnoses %>%
    dplyr::left_join( r_tbl_srprimarycaremedication, by = join_by( person_id ) ) %>% 
    # Select every record that has a prescription for any diabetes medication.
    dplyr::inner_join( names_meds_of_interest, by = join_by( nameofmedication == DMplusD_ProductDescription ) ) %>%
    # Filter for records within the follow-up period. I subtract `HMA_adjust_loockBack_window` I will need to 
    # know patients' prescriptions before follow-up when calculating the H.M.A. state.
    dplyr::filter( dplyr::between( dateevent, adjusted_date_followup_start, date_followup_end ) ) %>%
    # Create a new variable that is the concatenation of `nameofmedication` and `medicationdosage`.
    dplyr::mutate( meds_nameAndDose = paste0( nameofmedication, " DOSE: ", medicationdosage ) ) %>% 
    # Select the variables of interest.
    dplyr::select( person_id, date_diagnosis, dateevent, nameofmedication, meds_nameAndDose ) %>%
    dplyr::distinct() %>%
    # Filter for however many subsequent prescriptions were specified in `n_iteraions`, then number them.
    dplyr::group_by( person_id ) %>%
    dbplyr::window_order( person_id, dateevent ) %>%
    dplyr::mutate( new_prescp_day = if_else( dateevent != lag( dateevent), 1, 0 ) ) %>%
    tidyr::replace_na( list( new_prescp_day = 1 ) ) %>%
    dplyr::mutate( i_prescrp = cumsum( new_prescp_day ) ) %>%
    dplyr::filter( i_prescrp <= n_iterations + 1 ) %>% # The "+ 1" is added so that we always include the final "Unobserved" state.
    dplyr::mutate( event_name = paste0( "prescription_", i_prescrp ) ) %>% 
    dplyr::ungroup() %>%
    # Rename columns.
    dplyr::rename(
        start_dttm = dateevent
        ,event_value = nameofmedication
    ) %>%
    # Set the end time of each prescription as the start time of the previous one.
    dplyr::group_by( person_id ) %>%
    dbplyr::window_order( person_id, start_dttm ) %>%
    dplyr::mutate( end_dttm = dplyr::if_else( person_id != lag( person_id ), start_dttm, lead( start_dttm ) ) ) %>%
    dplyr::mutate( end_dttm = dplyr::if_else( is.na( end_dttm ), start_dttm, end_dttm ) ) %>%
    dplyr::ungroup() %>%
    # Tidy up.
    dplyr::select( person_id, start_dttm, event_name, event_value, end_dttm )


# Retrieve dates of tests in the follow-up period.
qry_log_test_longFormat <- 
    qry_records_with_T2DM_diagnoses %>%
    dplyr::left_join( r_tbl_srcode, by = join_by( person_id ) ) %>% 
    # Filter records for only those that refer to the test of interest.
    dplyr::filter( snomedcode %in% codes_SNOMED_test_of_interest ) %>% 
    # Extract the required fields.
    dplyr::select( person_id, dateevent, numericvalue, testType ) %>%
    dplyr::distinct() %>%
    # Filter out test scores of 0, which are expected to be data entry anomalies, given my selection of columns.
    dplyr::filter( numericvalue != "0" ) %>%
     # Filter for records within the follow-up period.
    dplyr::filter( dplyr::between( dateevent, date_followup_start, date_followup_end ) ) %>%
    # Rename columns.
    dplyr::rename(
        start_dttm = dateevent
        ,event_value_numeric = numericvalue
    ) %>% 
    # Filter for however many subsequent tests were specified in `n_iteraions`, then number them.
    dplyr::group_by( person_id ) %>% 
    dbplyr::window_order( start_dttm ) %>% 
    dplyr::filter( row_number() <= n_iterations + 1 ) %>% # The "+ 1" is added so that we always include the final "Unobserved" state.
    dplyr::mutate( event_name = paste0( "test_", row_number() - 1 ) ) %>% 
    dplyr::ungroup() %>% 
    # Set the end time each test as the start time of the previous one.
    dplyr::group_by( person_id ) %>%
    dbplyr::window_order( start_dttm ) %>%
    dplyr::mutate( end_dttm = lead( start_dttm ) ) %>% 
    dplyr::ungroup() %>%
    # Convert the numeric test values into labels.
    dplyr::mutate( event_value_numeric = as.numeric( event_value_numeric ) ) %>%
    dplyr::mutate(
        event_value = dplyr::case_when(
            is.na( event_value_numeric ) ~ NA_character_
            ,event_value_numeric > 70 ~ "Test Status = Red"
            ,event_value_numeric > 58  ~ "Test Status = Amber"
            ,event_value_numeric > 48 ~ "Test Status = Yellow"
            ,TRUE ~ "Test Status = Green"
            )
    ) %>%
    dplyr::rename( HbA1c = event_value_numeric ) %>%
    # Tidy up.
    dplyr::select( person_id, start_dttm, event_name, event_value, end_dttm, HbA1c, testType ) %>%
    dplyr::distinct()


# Create a variable to indicate the date from which we no longer observe the patient in the record.
patients_of_interest_with_death_dates <-
    r_tbl_srpatient %>%
    # ## Select only those patients in whom we are interested.
    dplyr::inner_join( qry_records_with_T2DM_diagnoses, by = join_by( person_id ) )
if( !is.na( dplyr::count(patients_of_interest_with_death_dates) %>% dplyr::pull() ) )
{
    df_date_unobserved <-
        patients_of_interest_with_death_dates %>%
        dplyr::distinct( person_id, datedeath ) %>% 
        dplyr::collect() %>%
        tidyr::separate(
            col = datedeath
            ,sep = 4
            ,into = c( "death_year", "death_month" )
        ) %>%
        dplyr::mutate(
            date_unobserved = 
                dplyr::if_else(
                    is.na( death_year )
                    ,NA
                    ,ISOdate( year = death_year, month = as.numeric( death_month) %>% `+`( 1 ) %>% as.character(), day = 1, hour = 0 )-1
                )
        ) %>%
        # If there are multiple dates of death, chose the earliest one.
        dplyr::group_by( person_id ) %>%
        dplyr::mutate(
            date_unobserved = 
                dplyr::if_else(
                    is.na( date_unobserved )
                    ,min( date_unobserved, na.rm = TRUE )
                    ,date_unobserved
                    )
        ) %>%
        dplyr::ungroup() %>%
        tidyr::drop_na() %>%
        # Tidy up.
        dplyr::select( person_id, date_unobserved ) %>%
        dplyr::arrange( person_id ) %>%
        suppressWarnings()
} else {
    df_date_unobserved <-
        data.frame(
            person_id = integer(0)
            ,date_unobserved = Date(0)
        )
}

# Append the prescriptions and tests dataframe logs.
df_log_PandT_longFormat <-
    qry_log_prescription_longFormat %>%
    dplyr::mutate( HbA1c = NA_integer_) %>%
    dplyr::union_all( qry_log_test_longFormat ) %>%
    dplyr::ungroup() %>% 
    dplyr::distinct() %>%
    dplyr::arrange( person_id, start_dttm ) %>%
    # Include a date from which we stop observing the patient in the record.
    # This will be either the date of death or one second after the `end_dttm`
    # of the last event in a patient's record.
    dplyr::collect() %>%
    dplyr::group_by( person_id ) %>%
    dplyr::group_modify( ~ tibble::add_row( .x ), .by = person_id ) %>% 
    dplyr::left_join( df_date_unobserved, by = join_by( person_id ) ) %>%
    dplyr::mutate( start_dttm = if_else( is.na( start_dttm ), date_unobserved , start_dttm ) ) %>%
    dplyr::ungroup() %>%
    # Replace the name of the medication in `event_value` with the short name created earlier.
    dplyr::left_join( shortNames_meds_of_interest, by = join_by( event_value == drug_name ) ) %>%
    dplyr::mutate( event_value = dplyr::if_else( is.na( drug_name_short ), event_value, drug_name_short ) ) %>%
    dplyr::select( -drug_name_short ) %>%
    # Update `event_name` so that it shows "Unobserved", where appropriate.
    dplyr::mutate( event_value = if_else( is.na( event_value ), "Unobserved", event_value ) ) %>%
    # Update `start_dttm` so that it shows one second after the last
    # `start_dttm` whenever `date_unobserved` is NA. 
    dplyr::mutate( start_dttm = if_else( is.na( start_dttm ), lag( start_dttm ) + lubridate::seconds(1), start_dttm ) ) %>%
    # Update `end_dttm` for the final observed event. It should no longer read
    # as NA because the "Unobserved" follows it.
    dplyr::mutate( end_dttm = if_else( lead( event_value ) == "Unobserved", lead( start_dttm ), end_dttm ) ) %>%
    suppressWarnings() %>%
    dplyr::select( - date_unobserved ) %>%
    dplyr::arrange( person_id, start_dttm )

# Include the date of diagnosis for every patient, and merge it with a test date if it is on that date.

                          

# Retrieve indication of diagnoses used in the calculation of comorbidity.
#
# In this iteration, I will use the Mayo Clinic's definition* of two or more of 19 specified conditions, where diabetes is one,
# separated by 30 days ([Rocca et al. (2014)](http://dx.doi.org/10.1016/j.mayocp.2014.07.010)). This means that I only need to
# identify patients that have at least one of the diagnostic codes in `codes_SNOMED_all_multimorbidity_diagnoses` more than
# 30 days before or after their diagnosis for Type 2 Diabetes Mellitus.
#
# *Rocca et al. (2014) used ICD-10 codes but I use SNOMED-CT codes.                          
qry_log_multimorb_longFormat <-
    r_tbl_srcode %>%
    # Filter for only those patients that we already identified in our prescription table.
    dplyr::inner_join( qry_log_prescription_longFormat %>% dplyr::select( person_id ), by = join_by( person_id ) ) %>%
    # Filter records for only those that refer to one of the diagnoses that constitute the multimorbidity variable.
    dplyr::filter( snomedcode %in% codes_SNOMED_all_multimorbidity_diagnoses ) %>% 
    # Extract the required fields.
    dplyr::select( person_id, dateevent, snomedcode ) %>%
    # Rename columns.
    dplyr::rename( date_multimorb = dateevent ) %>%
    # Include each person's date of diagnosis.
    dplyr::left_join( qry_records_with_T2DM_diagnoses, by = join_by( person_id ) ) %>%
    # Filter for multimorbidity diagnoses recorded outwith the multimorbidity gap window of the date of diagnosis of the main diagnosis of interest.
    dplyr::filter(
        ( sql( paste0( 'DATE_DIFF( date_multimorb, date_diagnosis, MONTH)' ) ) %>% abs() ) >= multimorb_gap_window_months
    ) %>%
    dplyr::select( - c( date_diagnosis, snomedcode ) ) %>%
    dplyr::distinct() %>% 
    # Filter for the earliest date of the multimorbdity diagnosis.
    dplyr::group_by( person_id ) %>%
    dplyr::filter(
        date_multimorb == min( date_multimorb )
    ) %>%
    dplyr::ungroup() %>%
    # Tidy up.
    dplyr::arrange( person_id, date_multimorb )          
              