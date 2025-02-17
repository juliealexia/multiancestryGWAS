library(tidyverse)
library(bigrquery)

# This query represents dataset "bcdata" for domain "person" and was generated for All of Us Controlled Tier Dataset v7
dataset_12212783_person_sql <- paste("
    SELECT
        person.person_id,
        person.gender_concept_id,
        p_gender_concept.concept_name as gender,
        person.birth_datetime as date_of_birth,
        person.race_concept_id,
        p_race_concept.concept_name as race,
        person.ethnicity_concept_id,
        p_ethnicity_concept.concept_name as ethnicity,
        person.sex_at_birth_concept_id,
        p_sex_at_birth_concept.concept_name as sex_at_birth 
    FROM
        `person` person 
    LEFT JOIN
        `concept` p_gender_concept 
            ON person.gender_concept_id = p_gender_concept.concept_id 
    LEFT JOIN
        `concept` p_race_concept 
            ON person.race_concept_id = p_race_concept.concept_id 
    LEFT JOIN
        `concept` p_ethnicity_concept 
            ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id 
    LEFT JOIN
        `concept` p_sex_at_birth_concept 
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id  
    WHERE
        person.PERSON_ID IN (SELECT
            distinct person_id  
        FROM
            `cb_search_person` cb_search_person  
        WHERE
            cb_search_person.person_id IN (SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `cb_search_all_events` 
                WHERE
                    (concept_id IN(SELECT
                        DISTINCT c.concept_id 
                    FROM
                        `cb_criteria` c 
                    JOIN
                        (SELECT
                            CAST(cr.id as string) AS id       
                        FROM
                            `cb_criteria` cr       
                        WHERE
                            concept_id IN (4112853)       
                            AND full_text LIKE '%_rank1]%'      ) a 
                            ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                            OR c.path LIKE CONCAT('%.', a.id) 
                            OR c.path LIKE CONCAT(a.id, '.%') 
                            OR c.path = a.id) 
                    WHERE
                        is_standard = 1 
                        AND is_selectable = 1) 
                    AND is_standard = 1 )) criteria ) )", sep="")

# Formulate a Cloud Storage destination path for the data exported from BigQuery.
# NOTE: By default data exported multiple times on the same day will overwrite older copies.
#       But data exported on a different days will write to a new location so that historical
#       copies can be kept as the dataset definition is changed.
person_12212783_path <- file.path(
  Sys.getenv("WORKSPACE_BUCKET"),
  "bq_exports",
  Sys.getenv("OWNER_EMAIL"),
  strftime(lubridate::now(), "%Y%m%d"),  # Comment out this line if you want the export to always overwrite.
  "person_12212783",
  "person_12212783_*.csv")
message(str_glue('The data will be written to {person_12212783_path}. Use this path when reading ',
                 'the data into your notebooks in the future.'))

# Perform the query and export the dataset to Cloud Storage as CSV files.
# NOTE: You only need to run `bq_table_save` once. After that, you can
#       just read data from the CSVs in Cloud Storage.
bq_table_save(
  bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), dataset_12212783_person_sql, billing = Sys.getenv("GOOGLE_PROJECT")),
  person_12212783_path,
  destination_format = "CSV")

# Read the data directly from Cloud Storage into memory.
# NOTE: Alternatively you can `gsutil -m cp {person_12212783_path}` to copy these files
#       to the Jupyter disk.
read_bq_export_from_workspace_bucket <- function(export_path) {
  col_types <- cols(gender = col_character(), race = col_character(), ethnicity = col_character(), sex_at_birth = col_character())
  bind_rows(
    map(system2('gsutil', args = c('ls', export_path), stdout = TRUE, stderr = TRUE),
        function(csv) {
          message(str_glue('Loading {csv}.'))
          chunk <- read_csv(pipe(str_glue('gsutil cat {csv}')), col_types = col_types, show_col_types = FALSE)
          if (is.null(col_types)) {
            col_types <- spec(chunk)
          }
          chunk
        }))
}
dataset_12212783_person_df <- read_bq_export_from_workspace_bucket(person_12212783_path)

bc_pheno <- phenos
bc_pheno$BC <- ifelse(bc_pheno$IID %in% dataset_12212783_person_df$person_id,1,0)
bc_pheno <- bc_pheno[ALL_cov$sex==1,]
print(paste0("ALL no.cases: ",sum(bc_pheno$BC)," no. total: ",nrow(bc_pheno)))
fwrite(bc_pheno,"ALL_BC.txt",sep="\t")
AFR_bc <- bc_pheno[bc_pheno$IID %in% AFR_cov$IID,]
print(paste0("AFR no.cases: ",sum(AFR_bc$BC)," no. total: ",nrow(AFR_bc)))
fwrite(AFR_bc,"AFR_BC.txt",sep="\t")
AMR_bc <- bc_pheno[bc_pheno$IID %in% AMR_cov$IID,]
print(paste0("AMR no.cases: ",sum(AMR_bc$BC)," no. total: ",nrow(AMR_bc)))
fwrite(AMR_bc,"AMR_BC.txt",sep="\t")
EAS_bc <- bc_pheno[bc_pheno$IID %in% EAS_cov$IID,]
print(paste0("EAS no.cases: ",sum(EAS_bc$BC)," no. total: ",nrow(EAS_bc)))
fwrite(EAS_bc,"EAS_BC.txt",sep="\t")
EUR_bc <- bc_pheno[bc_pheno$IID %in% EUR_cov$IID,]
print(paste0("EUR no.cases: ",sum(EUR_bc$BC)," no. total: ",nrow(EUR_bc)))
fwrite(EUR_bc,"EUR_BC.txt",sep="\t")
SAS_bc <- bc_pheno[bc_pheno$IID %in% SAS_cov$IID,]
print(paste0("SAS no.cases: ",sum(SAS_bc$BC)," no. total: ",nrow(SAS_bc)))
fwrite(SAS_bc,"SAS_BC.txt",sep="\t")
MID_bc <- bc_pheno[bc_pheno$IID %in% MID_cov$IID,]
print(paste0("MID no.cases: ",sum(MID_bc$BC)," no. total: ",nrow(MID_bc)))
fwrite(bc_pheno,"MID_BC.txt",sep="\t")
bc_cov <- covs[covs$sex==1]
bc_cov$sex <- NULL
fwrite(bc_cov,"ALL_BC_cov.txt",sep="\t")
AFR_bc_cov <- bc_cov[bc_cov$IID %in% AFR_cov$IID,]
fwrite(AFR_bc_cov,"AFR_BC_cov.txt",sep="\t")
AMR_bc_cov <- bc_cov[bc_cov$IID %in% AMR_cov$IID,]
fwrite(AMR_bc_cov,"AMR_BC_cov.txt",sep="\t")
EAS_bc_cov <- bc_cov[bc_cov$IID %in% EAS_cov$IID,]
fwrite(EAS_bc_cov,"EAS_BC_cov.txt",sep="\t")
EUR_bc_cov <- bc_cov[bc_cov$IID %in% EUR_cov$IID,]
fwrite(EUR_bc_cov,"EUR_BC_cov.txt",sep="\t")
SAS_bc_cov <- bc_cov[bc_cov$IID %in% SAS_cov$IID,]
fwrite(SAS_bc_cov,"SAS_BC_cov.txt",sep="\t")
MID_bc_cov <- bc_cov[bc_cov$IID %in% MID_cov$IID,]
fwrite(MID_bc_cov,"MID_BC_cov.txt",sep="\t")