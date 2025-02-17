#Height

dataset_84367389_person_sql = """
    SELECT
        measurement.PERSON_ID,
        measurement.MEASUREMENT_DATETIME,
        m_ext.src_id as DATA_SOURCE,
        measurement.MEASUREMENT_CONCEPT_ID, 
        m_standard_concept.concept_name as STANDARD_CONCEPT_NAME, 
        m_standard_concept.vocabulary_id as STANDARD_VOCABULARY, 
        measurement.UNIT_CONCEPT_ID, 
        m_unit.concept_name as UNIT_CONCEPT_NAME, 
        measurement.MEASUREMENT_TYPE_CONCEPT_ID, 
        m_type.concept_name as MEASUREMENT_TYPE_CONCEPT_NAME, 
        measurement.MEASUREMENT_SOURCE_CONCEPT_ID, 
        m_source_concept.concept_name as SOURCE_CONCEPT_NAME, 
        m_source_concept.vocabulary_id as SOURCE_VOCABULARY,
        measurement.VALUE_AS_NUMBER,
        measurement.VALUE_AS_CONCEPT_ID, 
        m_value.concept_name as VALUE_AS_CONCEPT_NAME, 
        m_operator.concept_name as OPERATOR_CONCEPT_NAME
    FROM
        `""" + os.environ["WORKSPACE_CDR"] + """.measurement` measurement 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.measurement_ext` m_ext 
            ON measurement.measurement_id = m_ext.measurement_id  
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_standard_concept 
            ON measurement.measurement_concept_id = m_standard_concept.concept_id  
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_unit 
            ON measurement.unit_concept_id = m_unit.concept_id  
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_type 
            ON measurement.measurement_type_concept_id = m_type.concept_id  
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_source_concept 
            ON measurement.measurement_source_concept_id = m_source_concept.concept_id  
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_value 
            ON measurement.value_as_concept_id = m_value.concept_id  
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_operator 
            ON measurement.operator_concept_id = m_operator.concept_id  
    WHERE
        (measurement_concept_id IN (3036277, 3023540, 3019171) OR 
    measurement_source_concept_id IN (903133))"""

dataset_84367390_person_df = pd.read_gbq(
    dataset_84367389_person_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")

dataset_84367390_person_df['MEASUREMENT_DATETIME'] = pd.to_datetime(dataset_84367390_person_df['MEASUREMENT_DATETIME'])
height_df = dataset_84367390_person_df.sort_values(by=['PERSON_ID', 'MEASUREMENT_DATETIME'])
most_recent_height = height_df.drop_duplicates(subset='PERSON_ID', keep='last')
array_height_info = pd.merge(dataset_sex_info, most_recent_height, left_on='person_id', right_on='PERSON_ID', how='left')
array_height_info['date_of_birth'] = pd.to_datetime(array_height_info['date_of_birth'])
array_height_info['age'] = ((array_height_info['MEASUREMENT_DATETIME'] - array_height_info['date_of_birth']).dt.days / 365.25).round()
#keep only one type of measurement metric for cohesion
array_height_info = array_height_info[array_height_info['UNIT_CONCEPT_NAME'].isin(['centimeter'])]
ID = pd.read_csv('ALL_cov_Null.txt',sep='\t')
ID = ID.iloc[:, :2]
height = pd.merge(ID,array_height_info[['person_id','age','VALUE_AS_NUMBER']], left_on='IID',right_on='person_id')
height = height.rename(columns={'VALUE_AS_NUMBER': 'height'})
height = height.dropna(subset=['height'])
height = height[height['height'] <= 240]
height = height[height['height'] > 0]
height = height.drop(columns=['person_id'])

ALLcov = pd.read_csv('ALL_cov_Null.txt', sep='\t')
AFRcov = pd.read_csv('AFR_cov_Null.txt', sep='\t')
AMRcov = pd.read_csv('AMR_cov_Null.txt', sep='\t')
EURcov = pd.read_csv('EUR_cov_Null.txt', sep='\t')
EAScov = pd.read_csv('EAS_cov_Null.txt', sep='\t')
SAScov = pd.read_csv('SAS_cov_Null.txt', sep='\t')
MIDcov = pd.read_csv('MID_cov_Null.txt', sep='\t')

#Remove old age column - as it is trait dependent
ALLcov = ALLcov.drop(columns=['age'])
AFRcov = AFRcov.drop(columns=['age'])
AMRcov = AMRcov.drop(columns=['age'])
EURcov = EURcov.drop(columns=['age'])
EAScov = EAScov.drop(columns=['age'])
SAScov = SAScov.drop(columns=['age'])
MIDcov = MIDcov.drop(columns=['age'])

#update for HEIGHT age
ALLcov_height = pd.merge(ALLcov,height[['IID','age']], on='IID')
AFRcov_height = pd.merge(AFRcov,height[['IID','age']], on='IID')
AMRcov_height = pd.merge(AMRcov,height[['IID','age']], on='IID')
EURcov_height = pd.merge(EURcov,height[['IID','age']], on='IID')
EAScov_height = pd.merge(EAScov,height[['IID','age']], on='IID')
SAScov_height = pd.merge(SAScov,height[['IID','age']], on='IID')
MIDcov_height = pd.merge(MIDcov,height[['IID','age']], on='IID')

ALLcov_height.to_csv('ALL_Height_cov.txt', sep='\t', index=False)
AFRcov_height.to_csv('AFR_Height_cov.txt', sep='\t', index=False)
AMRcov_height.to_csv('AMR_Height_cov.txt', sep='\t', index=False)
EURcov_height.to_csv('EUR_Height_cov.txt', sep='\t', index=False)
EAScov_height.to_csv('EAS_Height_cov.txt', sep='\t', index=False)
SAScov_height.to_csv('SAS_Height_cov.txt', sep='\t', index=False)
MIDcov_height.to_csv('MID_Height_cov.txt', sep='\t', index=False)

!gsutil -u $GOOGLE_PROJECT cp *_Height_cov.txt gs://fc-secure-35e2cef0-fb25-4753-a639-8e21bee4018b/data

heightpheno = height.drop(columns=['age'])
heightpheno.to_csv('ALL_Height.txt', sep='\t', index=False)
heightpheno_AFR = heightpheno[heightpheno['IID'].isin(AFRcov_height['IID'])]
heightpheno_AFR.to_csv('AFR_Height.txt', sep='\t', index=False)
heightpheno_AMR = heightpheno[heightpheno['IID'].isin(AMRcov_height['IID'])]
heightpheno_AMR.to_csv('AMR_Height.txt', sep='\t', index=False)
heightpheno_EUR = heightpheno[heightpheno['IID'].isin(EURcov_height['IID'])]
heightpheno_EUR.to_csv('EUR_Height.txt', sep='\t', index=False)
heightpheno_EAS = heightpheno[heightpheno['IID'].isin(EAScov_height['IID'])]
heightpheno_EAS.to_csv('EAS_Height.txt', sep='\t', index=False)
heightpheno_SAS = heightpheno[heightpheno['IID'].isin(SAScov_height['IID'])]
heightpheno_SAS.to_csv('SAS_Height.txt', sep='\t', index=False)
heightpheno_MID = heightpheno[heightpheno['IID'].isin(MIDcov_height['IID'])]
heightpheno_MID.to_csv('MID_Height.txt', sep='\t', index=False)

!gsutil -u $GOOGLE_PROJECT cp *_Height.txt gs://fc-secure-35e2cef0-fb25-4753-a639-8e21bee4018b/data


