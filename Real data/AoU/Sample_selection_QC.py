##Run only once
#lines = ["AFR_MRMEGA.txt", "AMR_MRMEGA.txt", "EUR_MRMEGA.txt", "EAS_MRMEGA.txt","SAS_MRMEGA.txt","MID_MRMEGA.txt"]

#with open("mrmega.in", "w") as file:
#    for line in lines:
#        file.write(line + "\n")

!gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv .
!gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv .

ancestry_pred = pd.read_csv("ancestry_preds.tsv", delimiter="\t")
related_samples = pd.read_csv("relatedness_flagged_samples.tsv", delimiter="\t")

# Convert the list_column into separate columns
df_split = ancestry_pred['pca_features'].str.split(',', expand=True)

# Rename the new columns appropriately
df_split.columns = [f'PC{i+1}' for i in range(df_split.shape[1])]
df_split['PC1'] = df_split['PC1'].str.slice(1)
df_split['PC16'] = df_split['PC16'].str[:-1]
ancestry_pred = pd.concat([ancestry_pred, df_split], axis=1)
ancestry_pred.drop(columns=['pca_features'], inplace=True)
ancestry_pred.to_csv('ancestry_tab.txt', sep='\t', index=False)

#Get sex_at_birth info
dataset_84367389_person_sql = """
    SELECT
        person.person_id,
        person.birth_datetime as date_of_birth,
        person.sex_at_birth_concept_id,
        p_sex_at_birth_concept.concept_name as sex_at_birth, 
    FROM
        `""" + os.environ["WORKSPACE_CDR"] + """.person` person 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_sex_at_birth_concept 
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id  
    WHERE
        person.PERSON_ID IN (SELECT
            distinct person_id  
        FROM
            `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` cb_search_person  
        WHERE
            cb_search_person.person_id IN (SELECT
                person_id 
            FROM
                `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p 
            WHERE
                has_array_data = 1 ) )"""

dataset_sex_info = pd.read_gbq(
    dataset_84367389_person_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")

dataset_sex_info.head(5)

##################
#SAMPLE SELECTION#
##################

eur_df = ancestry_pred[ancestry_pred['ancestry_pred_other'] == 'eur'].copy()
eur_df['FID']=0
eur = eur_df[['FID','research_id']]
eur.columns = ['FID', 'IID']
#123072 EUR
mask = ~eur['IID'].isin(related_samples['sample_id'])
eur_filtered = eur[mask]
#118149 EUR after removing related samples
eur_filtered.to_csv('EUR_IDs_sub.txt', sep='\t', index=False)
eur_sexfiltered = pd.merge(eur_filtered, dataset_sex_info, left_on='IID', right_on='person_id', how='left')
eur_sexfiltered = eur_sexfiltered[eur_sexfiltered['sex_at_birth'].isin(['Male','Female'])]
eur_sexfiltered = eur_sexfiltered.drop(columns=['person_id','date_of_birth','sex_at_birth_concept_id','sex_at_birth'])
#115708 EUR after removing non Male or Female samples
eur_sexfiltered.to_csv('EUR_IDs_filt.txt', sep='\t', index=False)

afr_df['FID']=0
afr = afr_df[['FID','research_id']]
afr.columns = ['FID', 'IID']
#53944 AFR
mask = ~afr['IID'].isin(related_samples['sample_id'])
afr_filtered = afr[mask]
#48419 AFR after removing related samples
afr_filtered.to_csv('AFR_IDs_sub.txt', sep='\t', index=False)
afr_sexfiltered = pd.merge(afr_filtered, dataset_sex_info, left_on='IID', right_on='person_id', how='left')
afr_sexfiltered = afr_sexfiltered[afr_sexfiltered['sex_at_birth'].isin(['Male','Female'])]
afr_sexfiltered = afr_sexfiltered.drop(columns=['person_id','date_of_birth','sex_at_birth_concept_id','sex_at_birth'])
#47210 AFR after removing non Male or Female samples
afr_sexfiltered.to_csv('AFR_IDs_filt.txt', sep='\t', index=False)

eas_df = ancestry_pred[ancestry_pred['ancestry_pred_other'] == 'eas'].copy()
eas_df['FID']=0
eas = eas_df[['FID','research_id']]
eas.columns = ['FID', 'IID']
#5381 EAS
mask = ~eas['IID'].isin(related_samples['sample_id'])
eas_filtered = eas[mask]
#5208 EAS after removing related samples
eas_filtered.to_csv('EAS_IDs_sub.txt', sep='\t', index=False)
eas_sexfiltered = pd.merge(eas_filtered, dataset_sex_info, left_on='IID', right_on='person_id', how='left')
eas_sexfiltered = eas_sexfiltered[eas_sexfiltered['sex_at_birth'].isin(['Male','Female'])]
eas_sexfiltered = eas_sexfiltered.drop(columns=['person_id','date_of_birth','sex_at_birth_concept_id','sex_at_birth'])
#5153 EAS after removing non Male or Female samples
eas_sexfiltered.to_csv('EAS_IDs_filt.txt', sep='\t', index=False)

amr_df = ancestry_pred[ancestry_pred['ancestry_pred_other'] == 'amr'].copy()
amr_df['FID']=0
amr = amr_df[['FID','research_id']]
amr.columns = ['FID', 'IID']
#40838 AMR
mask = ~amr['IID'].isin(related_samples['sample_id'])
amr_filtered = amr[mask]
#37125 AMR after removing related samples
amr_filtered.to_csv('AMR_IDs_sub.txt', sep='\t', index=False)
amr_sexfiltered = pd.merge(amr_filtered, dataset_sex_info, left_on='IID', right_on='person_id', how='left')
amr_sexfiltered = amr_sexfiltered[amr_sexfiltered['sex_at_birth'].isin(['Male','Female'])]
amr_sexfiltered = amr_sexfiltered.drop(columns=['person_id','date_of_birth','sex_at_birth_concept_id','sex_at_birth'])
#36502 AMR after removing non Male or Female samples
amr_sexfiltered.to_csv('AMR_IDs_filt.txt', sep='\t', index=False)

mid_df = ancestry_pred[ancestry_pred['ancestry_pred_other'] == 'mid'].copy()
mid_df['FID']=0
mid = mid_df[['FID','research_id']]
mid.columns = ['FID', 'IID']
#528 MID
mask = ~mid['IID'].isin(related_samples['sample_id'])
mid_filtered = mid[mask]
#513 MID after removing related samples
mid_filtered.to_csv('MID_IDs_sub.txt', sep='\t', index=False)
mid_sexfiltered = pd.merge(mid_filtered, dataset_sex_info, left_on='IID', right_on='person_id', how='left')
mid_sexfiltered = mid_sexfiltered[mid_sexfiltered['sex_at_birth'].isin(['Male','Female'])]
mid_sexfiltered = mid_sexfiltered.drop(columns=['person_id','date_of_birth','sex_at_birth_concept_id','sex_at_birth'])
#497 MID after removing non Male or Female samples
mid_sexfiltered.to_csv('MID_IDs_filt.txt', sep='\t', index=False)

sas_df = ancestry_pred[ancestry_pred['ancestry_pred_other'] == 'sas'].copy()
sas_df['FID']=0
sas = sas_df[['FID','research_id']]
sas.columns = ['FID', 'IID']
#2342 SAS
mask = ~sas['IID'].isin(related_samples['sample_id'])
sas_filtered = sas[mask]
#2274 SAS after removing related samples
sas_filtered.to_csv('SAS_IDs_sub.txt', sep='\t', index=False)
sas_sexfiltered = pd.merge(sas_filtered, dataset_sex_info, left_on='IID', right_on='person_id', how='left')
sas_sexfiltered = sas_sexfiltered[sas_sexfiltered['sex_at_birth'].isin(['Male','Female'])]
sas_sexfiltered = sas_sexfiltered.drop(columns=['person_id','date_of_birth','sex_at_birth_concept_id','sex_at_birth'])
#2247 SAS after removing non Male or Female samples
sas_sexfiltered.to_csv('SAS_IDs_filt.txt', sep='\t', index=False)

ids_tot = pd.concat([afr_sexfiltered, amr_sexfiltered, eas_sexfiltered, eur_sexfiltered, mid_sexfiltered, sas_sexfiltered], ignore_index=True)
ids_tot.to_csv('ALL_IDs_filt.txt', sep='\t', index=False)

!gsutil -u $GOOGLE_PROJECT -m cp gs://fc-aou-datasets-controlled/v7/microarray/plink_v7.1/arrays* plink/

%%writefile ID_extract.sh
#!/bin/bash

set -o errexit
set -o nounset

plink2 \
--bed "${BED}" \
--bim "${BIM}" \
--fam "${FAM}" \
--keep "${IDS}" \
--maf 0.01 \
--make-bed \
--out "${OUT_PATH}/${ETH}_filt"


%%writefile ID_extract_ALL.sh
#!/bin/bash

set -o errexit
set -o nounset

plink2 \
--bed "${BED}" \
--bim "${BIM}" \
--fam "${FAM}" \
--keep "${IDS}" \
--make-bed \
--out "${OUT_PATH}/${ETH}_filt"

%%writefile ID_extract_task.R

bucket = 'gs://fc-secure-35e2cef0-fb25-4753-a639-8e21bee4018b/data/'
plink_dir = 'gs://fc-aou-datasets-controlled/v7/microarray/plink_v7.1/'

tasks = data.frame(check.names = FALSE)

for (eth in c("AFR","AMR","EUR","EAS","SAS","MID")) {
    tasks = rbind(tasks, data.frame(
        '--input BED'=paste0(plink_dir,'arrays.bed'),
        '--input BIM'=paste0(plink_dir,'arrays.bim'),
        '--input FAM'=paste0(plink_dir,'arrays.fam'),
        '--input IDS'=paste0(bucket,eth,'_IDs_filt.txt'),
        '--env ETH'=eth,
        '--input-recursive INPUT_PATH'=paste0(bucket),
        '--output-recursive OUT_PATH'=paste0(bucket),
        check.names = FALSE
    ))
}

write.table(tasks, 
            file="ID_extract_task.txt", 
            row.names=F, col.names=T, 
            sep='\t', quote=F)

%%writefile ID_extract_task_ALL.R

bucket = 'gs://fc-secure-35e2cef0-fb25-4753-a639-8e21bee4018b/data/'
plink_dir = 'gs://fc-aou-datasets-controlled/v7/microarray/plink_v7.1/'

tasks = data.frame(check.names = FALSE)

for (eth in c("ALL")) {
    tasks = rbind(tasks, data.frame(
        '--input BED'=paste0(plink_dir,'arrays.bed'),
        '--input BIM'=paste0(plink_dir,'arrays.bim'),
        '--input FAM'=paste0(plink_dir,'arrays.fam'),
        '--input IDS'=paste0(bucket,eth,'_IDs_filt.txt'),
        '--env ETH'=eth,
        '--input-recursive INPUT_PATH'=paste0(bucket),
        '--output-recursive OUT_PATH'=paste0(bucket),
        check.names = FALSE
    ))
}

write.table(tasks, 
            file="ID_extract_task_ALL.txt", 
            row.names=F, col.names=T, 
            sep='\t', quote=F)


!Rscript ID_extract_task.R
!Rscript ID_extract_task_ALL.R

%%writefile QC_set.sh
#!/bin/bash

set -o errexit
set -o nounset

plink2 \
--bed "${BED}" \
--bim "${BIM}" \
--fam "${FAM}" \
--autosome \
--snps-only just-acgt \
--make-bed \
--out "${OUT_PATH}/${ETH}_QC"

%%writefile QC_task.R

bucket = 'gs://fc-secure-35e2cef0-fb25-4753-a639-8e21bee4018b/data/'
plink_dir = 'gs://fc-aou-datasets-controlled/v7/microarray/plink_v7.1/'

tasks = data.frame(check.names = FALSE)

for (eth in c("EUR")) {
    tasks = rbind(tasks, data.frame(
        '--input BED'=paste0(bucket,eth,'_af.bed'),
        '--input BIM'=paste0(bucket,eth,'_af.bim'),
        '--input FAM'=paste0(bucket,eth,'_af.fam'),
        '--env ETH'=eth,
        '--input-recursive INPUT_PATH'=paste0(bucket),
        '--output-recursive OUT_PATH'=paste0(bucket),
        check.names = FALSE
    ))
}

write.table(tasks, 
            file="QC_task.txt", 
            row.names=F, col.names=T, 
            sep='\t', quote=F)

!Rscript QC_task.R

