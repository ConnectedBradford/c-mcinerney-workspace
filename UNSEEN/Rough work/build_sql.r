cmhd_condition <- function() {

  ##############################################################################
  # Complex Mental Health Difficulties Condition Table
  #
  # Create table that tallies the records of the following complex mental health
  # difficulties:
  # 1. Bipolar disorder2509	1058567	0.237019
  # 2. Developmental academic disorder	1026	1058567	0.096923
  # 3. Schizophrenia
  # 4. Borderline personality disorder	496	1058567	0.046856
  # 5. Chronic depression	1815	1058567	0.171458
  # 6. Chronic post-traumatic stress disorder	116	1058567	0.010958
  # 7. Dysthymia	859	1058567	0.081147
  # 9. Personality disorder	1907	1058567	0.180149
  #
  ##############################################################################

  sql <- str_c(
    "CREATE TABLE `",project,".",target_dataset,".",target_table_prefix,"CMHD_Condition`
    AS
    SELECT DISTINCT
      p.person_id,
      v.condition_occurrence_id,
      v.condition_concept_id,
      c.concept_name AS Condition_concept_name,
      v.condition_source_value,
      c.concept_code,
      v.condition_start_date,
      p.AgeJan19,
      cv.care_site_id,
      cv.care_site_name,
      1 AS Delirium
    FROM
      `",project,".",cdm_source_dataset,".condition_occurrence` AS v
    Left JOIN
      `",project,".",cdm_source_dataset,".concept` AS c
      On v.condition_concept_id = c.concept_id
    Left JOIN
      `", project,".",target_dataset,".",target_table_prefix,"person_table` AS p
      ON v.person_id=p.person_id
    LEFT JOIN
      `", project,".",target_dataset,".",target_table_prefix,"Vist_Info` AS cv
      ON cv.visit_occurrence_id=v.visit_occurrence_id
    WHERE
      c.concept_id IN (", paste0(selected_snomed_codes_delirium, collapse = ", "),")
      AND
      v.condition_start_date BETWEEN DATE(\"",project_start_date,"\") AND DATE(\"",project_end_date,"\")
      AND
      cv.place_of_service_concept_id=4263714;"

  )
  return(sql)
}