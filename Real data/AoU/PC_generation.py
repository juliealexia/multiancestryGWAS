%%writefile LD_prune.sh
#!/bin/bash

set -o errexit
set -o nounset

plink2 \
--bed "${BED}" \
--bim "${BIM}" \
--fam "${FAM}" \
--indep-pairwise 500 kb 1 0.1 \
--out "${OUT_PATH}/prunedsnps"

%%writefile LD_prune_task.R

bucket = 'gs://fc-secure-35e2cef0-fb25-4753-a639-8e21bee4018b/data/'
plink_dir = 'gs://fc-aou-datasets-controlled/v7/microarray/plink_v7.1/'

tasks = data.frame(check.names = FALSE)

for (eth in c("EUR")) {
    tasks = rbind(tasks, data.frame(
        '--input BED'=paste0(bucket,'EUR_QC.bed'),
        '--input BIM'=paste0(bucket,'EUR_QC.bim'),
        '--input FAM'=paste0(bucket,'EUR_QC.fam'),
        '--input-recursive INPUT_PATH'=paste0(bucket),
        '--output-recursive OUT_PATH'=paste0(bucket),
        check.names = FALSE
    ))
}

write.table(tasks, 
            file="LD_prune_task.txt", 
            row.names=F, col.names=T, 
            sep='\t', quote=F)

!Rscript LD_prune_task.R

%%writefile VCF_prep.sh
#!/bin/bash

set -o errexit
set -o nounset

plink2 \
--bed "${BED}" \
--bim "${BIM}" \
--fam "${FAM}" \
--extract "${SNPS}" \
--export vcf \
--out "${OUT_PATH}/${ETH}_subset"

%%writefile VCF_prep.R

bucket = 'gs://fc-secure-35e2cef0-fb25-4753-a639-8e21bee4018b/data/'

tasks = data.frame(check.names = FALSE)

for (eth in c("ALL","MID","SAS","EAS","AFR","AMR","EUR")) {
    tasks = rbind(tasks, data.frame(
        '--input BED'=paste0(bucket,eth,'_filt.bed'),
        '--input BIM'=paste0(bucket,eth,'_filt.bim'),
        '--input FAM'=paste0(bucket,eth,'_filt.fam'),
        '--input SNPS'=paste0(bucket,'prunedsnps.prune.in'),
        '--env ETH'=eth,
        '--input-recursive INPUT_PATH'=paste0(bucket),
        '--output-recursive OUT_PATH'=paste0(bucket),
        check.names = FALSE
    ))
}

write.table(tasks, 
            file="VCF_prep.txt", 
            row.names=F, col.names=T, 
            sep='\t', quote=F)

!Rscript VCF_prep.R

%%writefile pca.sh
#!/bin/bash

set -o errexit
set -o nounset

plink2 \
--vcf "${VCF}" \
--make-bed \
--pca approx \
--out "${OUT_PATH}/${ETH}_pca"

%%writefile pca_task.R

bucket = 'gs://fc-secure-35e2cef0-fb25-4753-a639-8e21bee4018b/data/'

tasks = data.frame(check.names = FALSE)

for (eth in c("ALL","MID","SAS","EAS","AFR","AMR","EUR")) {
    tasks = rbind(tasks, data.frame(
        '--input VCF'=paste0(bucket,eth,'_subset.vcf'),
        '--env ETH'=eth,
        '--input-recursive INPUT_PATH'=paste0(bucket),
        '--output-recursive OUT_PATH'=paste0(bucket),
        check.names = FALSE
    ))
}

write.table(tasks, 
            file="pca_task.txt", 
            row.names=F, col.names=T, 
            sep='\t', quote=F)

!Rscript pca_task.R

#Reformat and add FID column
for eth in ['AFR','AMR','EAS','EUR','MID','SAS','ALL']:
    file_name = f'{eth}_pca.eigenvec'
    df = pd.read_csv(file_name, sep='\t')
    df.rename(columns={'#IID': 'IID'}, inplace=True)
    df.insert(0, 'FID', 0)
    df.to_csv(file_name, sep='\t', index=False)

#Look up sex and age to covariate file
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

reference_date = datetime(2024, 1, 1)

#Add sex and age info to make covariate file. Age is age on January 1st 2024 for null phenotype Type I error simulations
for eth in ['AFR','AMR','EAS','EUR','MID','SAS','ALL']:
    file_name = f'{eth}_pca.eigenvec'
    df = pd.read_csv(file_name, sep='\t')
    df = pd.merge(df,dataset_sex_info[['person_id','date_of_birth','sex_at_birth']], left_on='IID',right_on='person_id')
    df.drop(columns=['person_id'], inplace=True)
    df['date_of_birth'] = pd.to_datetime(df['date_of_birth'])
    df['age'] = ((pd.to_datetime(reference_date) - df['date_of_birth'].dt.tz_localize(None)).dt.days / 365.25).round()
    df.drop(columns=['date_of_birth'], inplace=True)
    df['sex'] = df['sex_at_birth'].apply(lambda x: 1 if x == 'Female' else 0)
    df.drop(columns=['sex_at_birth'], inplace=True)
    out_file_name = f'{eth}_cov_Null.txt'
    df.to_csv(out_file_name, sep='\t', index=False)

